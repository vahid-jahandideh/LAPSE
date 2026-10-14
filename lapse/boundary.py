"""Boundary relation algebra and conditional exposure probability.

Consequences obtained only by adding public correctness to a one-sided relation
do not affect additional shares. See ``docs/beta-definition.md`` for the precise
definition of the exposure event.
"""

from dataclasses import asdict, dataclass
from typing import Tuple

from lapse import core as rpt


def _binary(field):
    return ((isinstance(field, rpt.PrimeField) and field.p == 2)
            or (isinstance(field, rpt.GF2m) and field.m == 1))


def _pack(row):
    return sum((int(value) & 1) << i for i, value in enumerate(row))


def _binary_basis(rows):
    basis = {}
    for value in rows:
        while value:
            pivot = value.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = value
                break
            value ^= basis[pivot]
    return basis


def _binary_reduce(value, basis):
    for pivot in sorted(basis, reverse=True):
        if value >> pivot & 1:
            value ^= basis[pivot]
    return value


def _binary_dependencies(rows):
    basis, dependencies = {}, []
    for index, value in enumerate(rows):
        coefficients = 1 << index
        while value:
            pivot = value.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = (value, coefficients)
                break
            previous, combination = basis[pivot]
            value ^= previous
            coefficients ^= combination
        if value == 0:
            dependencies.append(coefficients)
    return dependencies


def _reduce(row, basis, field):
    row = list(row)
    for pivot in sorted(basis):
        coefficient = row[pivot]
        if not field.eq0(coefficient):
            row = [field.sub(a, field.mul(coefficient, b))
                   for a, b in zip(row, basis[pivot])]
    return row


def _basis(rows, field):
    basis = {}
    for row in rows:
        row = _reduce(row, basis, field)
        pivot = next((i for i, value in enumerate(row) if not field.eq0(value)), None)
        if pivot is not None:
            inverse = field.inv(row[pivot])
            basis[pivot] = [field.mul(inverse, value) for value in row]
    return basis


def _dependencies(rows, field):
    count = len(rows)
    basis, dependencies = {}, []
    for index, original in enumerate(rows):
        row = list(original)
        coefficients = [field.zero()] * count
        coefficients[index] = field.one()
        for pivot in sorted(basis):
            coefficient = row[pivot]
            if not field.eq0(coefficient):
                previous, combination = basis[pivot]
                row = [field.sub(a, field.mul(coefficient, b))
                       for a, b in zip(row, previous)]
                coefficients = [field.sub(a, field.mul(coefficient, b))
                                for a, b in zip(coefficients, combination)]
        pivot = next((i for i, value in enumerate(row) if not field.eq0(value)), None)
        if pivot is None:
            dependencies.append(coefficients)
        else:
            inverse = field.inv(row[pivot])
            basis[pivot] = ([field.mul(inverse, value) for value in row],
                            [field.mul(inverse, value) for value in coefficients])
    return dependencies


@dataclass(frozen=True)
class BoundaryRelations:
    """Bases use the coordinate order in ``boundary_shares``; RHSs are omitted."""

    boundary_shares: Tuple[str, ...]
    affected_shares: Tuple[str, ...]
    baseline_dimension: int
    relation_dimension: int
    input_only_dimension: int
    output_only_dimension: int
    coupled_dimension: int
    native_recovered: bool
    relation_basis: Tuple[Tuple[int, ...], ...]
    input_only_basis: Tuple[Tuple[int, ...], ...]
    output_only_basis: Tuple[Tuple[int, ...], ...]


@dataclass(frozen=True)
class BetaEstimate(rpt.RecoveryEstimate):
    """Exposure estimate conditioned on the selected boundary site not leaking."""

    boundary_share: str = ""
    conditional: bool = True
    p_effective: float = 0.0

    @property
    def beta(self):
        return self.estimate


class _BoundaryModel:
    def __init__(self, ir):
        self.ir = ir
        self.field = ir.field
        self.inputs, self.outputs = rpt.candidates_inputs_outputs(ir)
        self.names = tuple(self.inputs + self.outputs)
        self.n = len(self.inputs)
        if self.n == 0 or len(self.outputs) != self.n:
            raise ValueError("Boundary analysis requires equally sized nonempty input and output sharings")
        declared_outputs = sum(1 if symbol.size is None else symbol.size
                               for symbol in ir.symbols.values() if symbol.role == 'output')
        if declared_outputs != self.n:
            raise ValueError("Every declared output share must be assigned before boundary analysis")
        self.width = len(ir.sources)
        self.rows = [list(ir.vectors[name]) + [self.field.zero()] * (self.width - len(ir.vectors[name]))
                     for name in self.names]
        self.native = [self.field.zero()] * self.width
        output_sum = [self.field.zero()] * self.width
        for row in self.rows[:self.n]:
            self.native = [self.field.add(a, b) for a, b in zip(self.native, row)]
        for row in self.rows[self.n:]:
            output_sum = [self.field.add(a, b) for a, b in zip(output_sum, row)]
        if self.native != output_sum:
            raise ValueError("Boundary correctness fails: input and output share sums differ")
        if rpt.rank_over_field(self.rows, self.field) != 2 * self.n - 1:
            raise ValueError("Boundary analysis requires a fresh refresh: correctness must be its only leak-free relation")
        if rpt.rank_over_field(self.rows[:self.n], self.field) != self.n:
            raise ValueError("Input shares must be independent source coordinates")
        self.binary = _binary(self.field)
        if self.binary:
            self.packed = [_pack(row) for row in self.rows]
            self.packed_native = _pack(self.native)
            self.packed_wires = {name: _pack(row) for name, row in ir.vectors.items()}

    def resolve(self, name):
        if name in self.names:
            return name
        matches = [candidate for candidate in self.names if candidate.split('#', 1)[0] == name]
        if len(matches) != 1:
            raise ValueError(f"Unknown or ambiguous boundary share {name!r}; choose one of {', '.join(self.names)}")
        return matches[0]

    def analyze(self, leaked, *, details=False, target=None):
        if self.binary:
            basis = _binary_basis([self.packed_wires[name] for name in leaked])
            residuals = [_binary_reduce(row, basis) for row in self.packed]
            relations = _binary_dependencies(residuals)
            inputs = _binary_dependencies(residuals[:self.n])
            outputs = _binary_dependencies(residuals[self.n:])
            known = _binary_reduce(self.packed_native, basis) == 0
            input_support, output_support = 0, 0
            for row in inputs:
                input_support |= row
            for row in outputs:
                output_support |= row
            support = input_support | (output_support << self.n)
        else:
            rows = [list(self.ir.vectors[name]) + [self.field.zero()] *
                    (self.width - len(self.ir.vectors[name])) for name in leaked]
            basis = _basis(rows, self.field)
            residuals = [_reduce(row, basis, self.field) for row in self.rows]
            relations = _dependencies(residuals, self.field)
            inputs = _dependencies(residuals[:self.n], self.field)
            outputs = _dependencies(residuals[self.n:], self.field)
            known = all(self.field.eq0(value) for value in _reduce(self.native, basis, self.field))
            support = 0
            for offset, side in ((0, inputs), (self.n, outputs)):
                for row in side:
                    for index, value in enumerate(row):
                        if not self.field.eq0(value):
                            support |= 1 << (offset + index)
        coupled = len(relations) - len(inputs) - len(outputs) - (0 if known else 1)
        if coupled < 0:
            raise ArithmeticError("Inconsistent boundary quotient dimensions")
        if coupled > 0:
            support = (1 << len(self.names)) - 1
        if not details:
            return bool(support >> target & 1)
        if self.binary:
            unpack = lambda value, width: tuple((value >> i) & 1 for i in range(width))
            relations = [unpack(row, 2 * self.n) for row in relations]
            inputs = [unpack(row, self.n) for row in inputs]
            outputs = [unpack(row, self.n) for row in outputs]
        zeros = (self.field.zero(),) * self.n
        return BoundaryRelations(
            self.names,
            tuple(name for i, name in enumerate(self.names) if support >> i & 1),
            1, len(relations), len(inputs), len(outputs), coupled, known,
            tuple(tuple(row) for row in relations),
            tuple(tuple(row) + zeros for row in inputs),
            tuple(zeros + tuple(row) for row in outputs),
        )


def boundary_shares(ir):
    """List original input and final output SSA wires of a validated refresh."""
    return _BoundaryModel(ir).names


def inspect_boundary_relations(ir, leaked_ssa):
    """Inspect exact relation spaces for an explicit subset of SSA leakage sites."""
    model = _BoundaryModel(ir)
    leaked_ssa = list(leaked_ssa)
    allowed = set(rpt.filter_leakables(ir.created, ir.secret_ssa, None))
    unknown = set(leaked_ssa) - allowed
    if unknown:
        raise ValueError(f"Unknown or nonleakable SSA sites: {', '.join(sorted(unknown))}")
    return model.analyze(leaked_ssa, details=True)


def estimate_beta(ir, p, boundary_share, *, allow_roles=None,
                  method="auto", samples=10000, max_enum=20, q=None,
                  seed=None, rng=None, pilot_samples=1000):
    """Estimate affection conditioned on the selected SSA site not leaking.

    Exactly that site is removed from independent Bernoulli sampling; copies
    with the same value remain separate sites. ``p=1`` is undefined because the
    conditioning event has probability zero. The exposure event need not be
    monotone; both MC and ordinary likelihood-ratio IS still apply.
    """
    p = rpt._validate_probability(p, 'p')
    if p == 1.0:
        raise ValueError("beta is undefined at p=1: the conditioning event has probability zero")
    model = _BoundaryModel(ir)
    selected = model.resolve(boundary_share)
    candidates = rpt.filter_leakables(ir.created, ir.secret_ssa, allow_roles)
    if selected not in candidates:
        raise ValueError("The selected boundary share must be included in the allowed leakage roles")
    candidates = [name for name in candidates if name != selected]
    target = model.names.index(selected)

    def event(indices):
        return model.analyze([candidates[index] for index in indices], target=target)

    result = rpt.estimate_event(event, len(candidates), p, method=method, samples=samples,
                                max_enum=max_enum, q=q, seed=seed, rng=rng,
                                pilot_samples=pilot_samples)
    return BetaEstimate(**asdict(result), boundary_share=selected,
                        conditional=True, p_effective=p + (1.0 - p) * result.estimate)
