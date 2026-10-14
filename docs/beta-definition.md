# Boundary exposure β

This page explains the boundary exposure rate β: what it measures, how LAPSE
decides which boundary shares a leakage instance affects, and how to compute β
with the tool.

## What β measures

A refresh gadget takes an input sharing $`\mathbf X=(X_1,\ldots,X_n)`$ and
returns a fresh output sharing $`\mathbf Y=(Y_1,\ldots,Y_n)`$ of the same
native, so that $`\sum_i X_i=\sum_i Y_i`$. When wires inside the refresh leak,
the attacker may learn linear relations among the boundary shares. A boundary
share that takes part in such a relation is *affected*.

For a boundary share $`U`$ and probing rate $`p`$, LAPSE computes

```math
\beta_U(p)=\Pr\bigl[\,U\in\operatorname{Aff}(\mathcal L)\ \big|\ U\notin\mathcal L\,\bigr],
```

the probability that $`U`$ is affected although $`U`$ itself did not leak.
Adding direct leakage gives the effective rate at which $`U`$ is exposed:

```math
p' = p + (1-p)\,\beta_U(p).
```

These are the exposure rate $`\beta(n,p)`$ and the inflated rate $`p'`$ of the
refresh-boundary exposure reduction in Section 6.2 of the paper. LAPSE reports
both.

## Quick start

```bash
lapse --file two_share.circ --n 2 --beta-p 0.1 --boundary-share 'a[0]' --method exact
```

```text
LAPSE: characteristic 2, 3 sources, 5 probeable wires
beta(p=0.1, U=a[0]#0) = 0.0919
Effective marginal rate p' = 0.18271
```

From Python:

```python
from pathlib import Path
from lapse import build_ir, estimate_beta, parse_program

ir = build_ir(parse_program(Path("two_share.circ").read_text()), {"n": 2})
result = estimate_beta(ir, 0.1, "a[0]", method="exact")
print(result.beta, result.p_effective)  # about 0.0919 and 0.18271
```

## Which circuits can be analyzed

β is defined for refresh gadgets. LAPSE checks that the circuit

- has nonempty input and output sharings of the same size, with every output
  share assigned;
- preserves the sum: the input shares and the output shares add up to the same
  native;
- has correctness as its only boundary relation when nothing leaks, as a fresh
  refresh should (formally, the $`2n`$ boundary vectors have rank $`2n-1`$).

A circuit that fails a check is rejected with an explanatory error. Of the
shipped examples, `two_share.circ`, `multi_refresh.circ`, `SNI_Refresh.circ`,
and `SNI_Refresh_field.circ` qualify. `Compress.circ` does not, because it maps
$`n^2`$ inputs to $`n`$ outputs.

The boundary shares are SSA wires. The input share is the first version of
each input element, so `a[0]` means `a[0]#0`. The output share is the last
version assigned to each output element. You can pass either the short name
(`a[0]`, `b[1]`) or the full SSA name.

## When is a share affected?

### The rule in words

Given the leaked wires, LAPSE finds every linear relation among the boundary
shares whose value the attacker can compute. It sorts these relations using
the paper's terms:

- **Input-only and output-only relations** affect exactly the shares they
  involve. For example, if $`X_2`$ leaks, the relation $`X_2=c`$ affects only
  $`X_2`$.
- **Bridges** involve both input and output shares. If a bridge can be written
  as a combination of one-sided relations and the correctness relation, it
  adds no new information and affects no further shares. Otherwise LAPSE treats
  it as an *essential* bridge and, conservatively, marks **every** boundary
  share as affected.
- **The correctness relation** $`\sum_i X_i-\sum_i Y_i=0`$ is public and never
  counts as leakage.

If the leakage reveals the native, then $`\sum_i X_i`$ and $`\sum_i Y_i`$ are
both known one-sided relations, so every boundary share is affected.

These three outcomes follow the structure of the paper's analysis of
$`\mathsf R`$-multi (the corollary bounding $`\beta(n,p)`$ in Section 6.3):
the native is revealed, an essential bridge occurs, or $`U`$ lies in a
one-sided relation.

### Why correctness is set aside

Combining a known relation with correctness does not change what the attacker
knows, but it does change which shares appear. For $`n=2`$, if $`X_1`$ leaks,
the attacker knows $`X_1=c`$. Subtracting correctness gives
$`-X_2+Y_1+Y_2=c`$, which involves the three other shares. Counting such
combinations would make every share look affected whenever anything leaks.
LAPSE therefore compares relations modulo correctness, and the result does
not depend on which equivalent form of a relation one writes down.

### Formal definition

Let $`B`$ be the matrix whose $`2n`$ rows are the vectors of
$`X_1,\ldots,X_n,Y_1,\ldots,Y_n`$ over the source variables (input shares and
fresh randomness). For a leakage instance whose leaked wire vectors form the
matrix $`L`$, the known boundary relations are the coefficient vectors

```math
H_L=\{a\in\mathbb F_q^{2n} : aB\in\operatorname{rowspan}(L)\}.
```

Their right-hand sides are known values. Public constants in the circuit do not
change $`H_L`$. Inside $`H_L`$, LAPSE distinguishes:

| Space | Meaning |
| --- | --- |
| $`H_0=\operatorname{span}\{(1,\ldots,1,-1,\ldots,-1)\}`$ | correctness, always contained in $`H_L`$ |
| $`I_L`$ | relations in $`H_L`$ with a zero output part (input-only) |
| $`O_L`$ | relations in $`H_L`$ with a zero input part (output-only) |
| $`C_L=H_L/(H_0+I_L+O_L)`$ | essential bridges: what correctness and one-sided relations do not explain |

The affected set is

```math
\operatorname{Aff}(L)=
\begin{cases}
\operatorname{supp}(I_L)\cup\operatorname{supp}(O_L), & \dim C_L=0,\\[0pt]
\mathbf X\cup\mathbf Y, & \dim C_L>0.
\end{cases}
```

Here the support of a subspace is the set of coordinates on which at least one
of its vectors is nonzero. All of these spaces depend only on what the leakage
determines, not on the order of the leaked wires or on elimination pivots.
LAPSE computes them exactly by Gaussian elimination over the circuit's field.

Why the rule is sufficient: given the native, the boundary shares are uniformly
distributed over the solutions of the known relations. When $`\dim C_L=0`$,
those relations split into separate input and output constraints, so the
conditional distribution is a product of an input part and an output part, and
exposing the shares in the supports of $`I_L`$ and $`O_L`$ reveals everything
those constraints say. When
$`\dim C_L>0`$, the constraints do not split, and exposing all boundary shares
covers the instance.

## Worked example: the two-share refresh

Take $`Y_1=X_1+R`$ and $`Y_2=X_2-R`$. (`two_share.circ` is the
$`\mathbb F_2`$ version, where subtraction is addition.) The circuit has five
leakage sites: $`X_1, X_2, R, Y_1, Y_2`$. Choose $`U=X_1`$. Conditioning on
$`X_1`$ not leaking leaves four sites, each leaking independently with
probability $`p`$.

The table shows several leakage instances and the affected sets LAPSE assigns:

| Leaked wires | Known boundary relations | Outcome | Affected shares |
| --- | --- | --- | --- |
| none | correctness only | nothing learned | none |
| $`X_2`$ | $`X_2=c`$ | input-only | $`X_2`$ |
| $`R`$ | $`Y_1-X_1=c`$, $`Y_2-X_2=-c`$ | essential bridge | all four |
| $`R, Y_1`$ | $`X_1=c`$, $`Y_1=c'`$ | one-sided | $`X_1, Y_1`$ |
| $`X_2, Y_1`$ | $`X_2=c`$, $`Y_1=c'`$ | one-sided; their sum is a bridge explained by them | $`X_2, Y_1`$ |
| $`X_2, R`$ | $`X_2=c`$, $`Y_2=c'`$, $`Y_1-X_1=c''`$ | the bridge equals $`X_2-Y_2`$ by correctness, so it is explained | $`X_2, Y_2`$ |
| $`Y_1, Y_2`$ | $`Y_1`$, $`Y_2`$, hence the native | native revealed | all four |

Collecting every instance in which $`X_1`$ is affected:

| Leaked wires among $`X_2,R,Y_1,Y_2`$ | Reason | Probability |
| --- | --- | --- |
| $`R`$ only | essential bridge | $`p(1-p)^3`$ |
| $`R, Y_1`$ | $`X_1`$ becomes known | $`p^2(1-p)^2`$ |
| $`Y_1, Y_2`$ | native revealed | $`p^2(1-p)^2`$ |
| any three of them, except $`X_2, R, Y_2`$ | native revealed | $`3p^3(1-p)`$ |
| all four | native revealed | $`p^4`$ |

Adding these probabilities gives

```math
\beta_{X_1}(p)=p-p^2+2p^3-p^4 .
```

At $`p=0.1`$ this is $`0.0919`$, so $`p'=0.1+0.9\cdot 0.0919=0.18271`$, as
reported in the quick start. The same value results over any field.

Compare the rows for $`R`$ and for $`R, Y_1`$: the additional leaked wire turns
the bridge into two one-sided relations, so the affected set becomes smaller.
This is expected under the rule above, and the estimators below handle it.

## How LAPSE estimates β

- **Conditioning.** LAPSE removes exactly the site of $`U`$ from sampling, and
  every other allowed site leaks independently at rate $`p`$. This is the
  conditional distribution given that $`U`$ does not leak. At $`p=1`$ the
  condition has probability zero, so β is undefined and LAPSE reports an
  error.
- **Copies are separate sites.** Another wire with the same value as $`U`$
  remains a separate leakage site. In `multi_refresh.circ`, each output `b[i]`
  is a copy of the last version of `a[i]`; when that wire leaks, `b[i]` is
  known and therefore affected. For that circuit, β of an output share
  includes this event, which alone has probability $`p`$. If the output share
  and the last-round wire should be a single wire, write the last round's
  result directly into the output array.
- **Leakage roles.** `--allow-roles` restricts which roles may leak. The role of
  $`U`$ must be among them.
- **Methods.** β uses the same estimators as $`\mathsf E_X`$: `exact`, `mc`,
  `is`, and `auto`, with the same sample, seed, and pilot options. On the CLI,
  `--beta-q` sets the importance-sampling proposal rate for β. In `lapse-sweep`,
  `--mc-q` sets it for both quantities. Monte Carlo and importance sampling
  remain unbiased even though the affected set can shrink when more wires leak.
  The automatic choice of the proposal rate is a heuristic, so check the
  reported warnings and effective sample size.

## Using β in the composition bound

β is computed for the share you choose, and different shares of the same
gadget can have different values. For `multi_refresh.circ` with $`n=3`$, one
round, and $`p=0.05`$, exact enumeration gives:

| Share | `a[0]` | `a[1]`, `a[2]` | `b[0]` | `b[1]`, `b[2]` |
| --- | --- | --- | --- | --- |
| $`\beta_U`$ | 0.173 | 0.167 | 0.214 | 0.208 |

To obtain a single rate $`p'`$ for a refresh, compute β for every boundary
share and use the largest value, as the paper's bound for $`\mathsf R`$-multi
takes the maximum over shares:

```python
from lapse import boundary_shares

p = 0.05
ir = build_ir(parse_program(Path("multi_refresh.circ").read_text()), {"n": 3, "rounds": 1})
worst = max(estimate_beta(ir, p, u, method="exact").beta for u in boundary_shares(ir))
p_prime = p + (1 - p) * worst  # about 0.214 and 0.253
```

As in the paper's composition theorem, $`p'`$ is then applied to each boundary
share, with the exposures of different shares treated as independent. LAPSE
reports the per-share rate. When you publish values, include the circuit source
and parameters, the field, the selected share, the leakage roles, and the
method and seed. The JSON output records the parameters and options together
with a hash of the circuit source; sweep metadata also stores the source.

## Reading the output

The CLI JSON `beta` object contains:

| Field | Meaning |
| --- | --- |
| `estimate` | $`\beta_U(p)`$ |
| `p_effective` | $`p'=p+(1-p)\beta_U(p)`$ |
| `boundary_share` | the resolved SSA name of $`U`$ |
| `wire_count` | number of sampled sites, excluding $`U`$ |
| `method`, `samples`, `hits`, `standard_error`, `confidence_interval`, `proposal_rate`, `warnings` | estimator diagnostics, as for $`\mathsf E_X`$ |

`lapse-sweep --metric beta` (or `recovery-beta` for both quantities) writes the
estimate to the CSV column `beta` and the diagnostics to `beta_*` columns.

To inspect a single leakage instance, call
`inspect_boundary_relations(ir, leaked)` with full SSA names, for example
`["r#0", "b[0]#0"]`. It returns the dimensions
of $`H_L`$, $`I_L`$, $`O_L`$, and $`C_L`$ (`relation_dimension`,
`input_only_dimension`, `output_only_dimension`, `coupled_dimension`), whether
the native is recovered, bases for the relation spaces, and `affected_shares`.
