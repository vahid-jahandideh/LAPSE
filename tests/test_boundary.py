"""Analytic oracles and field-independent checks for explicit beta policies."""

import itertools
import unittest

from lapse import core as rpt
from lapse.boundary import boundary_shares, estimate_beta, inspect_boundary_relations


def refresh(field="bit", n=2):
    lines = ["circuit refresh();", f"input {field} a[{n}];",
             f"rand {field} r[{n - 1}];", f"output {field} b[{n}];"]
    for i in range(n - 1):
        lines.append(f"b[{i}] = a[{i}] + r[{i}];")
    expression = f"a[{n - 1}]" + "".join(f" - r[{i}]" for i in range(n - 1))
    lines.extend((f"b[{n - 1}] = {expression};",
                  f"const {field} S = sum_reduce(a, 0, {n});"))
    return rpt.build_ir(rpt.parse_program("\n".join(lines)), {})


class BoundaryTests(unittest.TestCase):
    def test_no_leakage_has_only_correctness(self):
        for field in ("bit", "gf(257)", "gf2m(2, 7)"):
            ir = refresh(field, 3)
            for policy in ("literal", "separated"):
                result = inspect_boundary_relations(ir, [], policy=policy)
                self.assertEqual(result.baseline_dimension, 1)
                self.assertEqual(result.relation_dimension, 1)
                self.assertEqual(result.coupled_dimension, 0)
                self.assertEqual(result.affected_shares, ())
                self.assertEqual(estimate_beta(ir, 0, 'a[0]', policy=policy).estimate, 0)

    def test_one_input_leak_does_not_affect_other_shares(self):
        ir = refresh(n=3)
        result = inspect_boundary_relations(ir, ['a[0]#0'], policy='separated')
        self.assertEqual(result.affected_shares, ('a[0]#0',))
        self.assertEqual(result.input_only_dimension, 1)
        self.assertEqual(result.output_only_dimension, 0)
        self.assertEqual(result.coupled_dimension, 0)
        literal = inspect_boundary_relations(ir, ['a[0]#0'], policy='literal')
        self.assertEqual(literal.affected_shares, boundary_shares(ir))

    def test_mask_leak_creates_coupled_residual(self):
        for field in ("bit", "gf(257)", "gf2m(2, 7)"):
            ir = refresh(field, 3)
            result = inspect_boundary_relations(ir, ['r[0]#0'], policy='separated')
            self.assertEqual(result.coupled_dimension, 1)
            self.assertEqual(result.input_only_dimension, 0)
            self.assertEqual(result.output_only_dimension, 0)
            self.assertEqual(result.affected_shares, boundary_shares(ir))

    def test_independent_input_and_output_leaks_are_separable(self):
        for field in ("bit", "gf(257)", "gf2m(2, 7)"):
            ir = refresh(field, 3)
            result = inspect_boundary_relations(ir, ['a[0]#0', 'b[1]#0'], policy='separated')
            self.assertEqual(result.coupled_dimension, 0)
            self.assertEqual(result.affected_shares, ('a[0]#0', 'b[1]#0'))

    def test_exact_two_share_analytic_probability(self):
        for field in ("bit", "gf(257)", "gf2m(2, 7)"):
            ir = refresh(field)
            for p in (0.01, 0.2, 0.5, 0.9):
                with self.subTest(field=field, p=p):
                    separated = estimate_beta(ir, p, 'a[0]', policy='separated', method='exact')
                    self.assertAlmostEqual(separated.estimate, p - p*p + 2*p**3 - p**4)
                    self.assertAlmostEqual(separated.p_effective, p + (1-p)*separated.estimate)
                    self.assertEqual(separated.wire_count, 4)
                    self.assertEqual(separated.boundary_share, 'a[0]#0')
                    self.assertTrue(separated.conditional)
                    literal = estimate_beta(ir, p, 'a[0]', policy='literal', method='exact')
                    self.assertAlmostEqual(literal.estimate, 1 - (1-p)**4)

    def test_identical_boundary_copy_remains_a_distinct_leakage_site(self):
        ir = rpt.build_ir(rpt.parse_program(
            "circuit copy(); input bit a; output bit b = a; const bit S = a;"), {})
        for target in ('a', 'b'):
            result = estimate_beta(ir, 0.2, target, policy='separated', method='exact')
            self.assertAlmostEqual(result.estimate, 0.2)
            self.assertEqual(result.wire_count, 1)

    def test_separated_policy_is_explicitly_nonmonotone(self):
        ir = refresh()
        less = inspect_boundary_relations(ir, ['r[0]#0'], policy='separated')
        more = inspect_boundary_relations(ir, ['r[0]#0', 'a[0]#0'], policy='separated')
        self.assertIn('a[1]#0', less.affected_shares)
        self.assertNotIn('a[1]#0', more.affected_shares)
        self.assertEqual(more.affected_shares, ('a[0]#0', 'b[0]#0'))

    def test_relation_bases_match_independently_enumerated_gf2_space(self):
        ir = refresh()
        names = boundary_shares(ir)
        wires = rpt.filter_leakables(ir.created, ir.secret_ssa, None)
        for mask in range(1 << len(wires)):
            leaked = [wires[i] for i in range(len(wires)) if mask >> i & 1]
            result = inspect_boundary_relations(ir, leaked, policy='separated')
            leaked_rows = [ir.vectors[name] for name in leaked]
            brute = []
            for coefficients in itertools.product((0, 1), repeat=len(names)):
                vector = [0] * len(ir.sources)
                for coefficient, name in zip(coefficients, names):
                    if coefficient:
                        vector = [a ^ b for a, b in zip(vector, ir.vectors[name])]
                if rpt.in_span(leaked_rows, vector, ir.field):
                    brute.append(coefficients)
            self.assertEqual(len(brute), 1 << result.relation_dimension)
            self.assertTrue(all(row in brute for row in result.relation_basis))
            self.assertEqual(rpt.rank_over_field(list(map(list, result.relation_basis)), ir.field),
                             result.relation_dimension)

    def test_sampling_is_reproducible_and_matches_exact(self):
        ir = refresh()
        exact = estimate_beta(ir, 0.2, 'a[0]', policy='separated', method='exact')
        for method in ('mc', 'is'):
            kwargs = dict(policy='separated', method=method, samples=10000, seed=71, q=0.3)
            result = estimate_beta(ir, 0.2, 'a[0]', **kwargs)
            self.assertEqual(result, estimate_beta(ir, 0.2, 'a[0]', **kwargs))
            self.assertLess(abs(result.estimate-exact.estimate), 6*result.standard_error)

    def test_reject_invalid_refresh_and_conditioning(self):
        ir = refresh()
        for kwargs in ({'p': 1, 'boundary_share': 'a[0]', 'policy': 'separated'},
                       {'p': .2, 'boundary_share': 'missing', 'policy': 'separated'},
                       {'p': .2, 'boundary_share': 'a[0]', 'policy': 'ambiguous'},
                       {'p': .2, 'boundary_share': 'a[0]', 'policy': 'separated',
                        'allow_roles': {'rand', 'output'}}):
            with self.assertRaises(ValueError):
                estimate_beta(ir, **kwargs)
        for source in (
            "input bit a[2]; output bit b[2]; b[0]=a[0]; b[1]=a[1];",
            "input bit a[2]; output bit b[2]; b[0]=a[0]; b[1]=0;",
            "input bit a[2]; output bit b[2]; b[0]=a[0];",
        ):
            invalid = rpt.build_ir(rpt.parse_program('circuit invalid(); '+source), {})
            with self.assertRaises(ValueError):
                boundary_shares(invalid)


if __name__ == '__main__':
    unittest.main()
