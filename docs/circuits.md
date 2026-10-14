# Circuit catalog and leakage conventions

Each assignment creates an SSA wire. A copy is a separate probeable wire, even
when its vector equals an existing vector. The default leakage set includes
input shares, random draws, intermediate assignments, and output assignments;
the target `const S` is excluded. Expressions are linearized as a single
assignment: intermediate operations inside one expression are not separate
probeable wires. Circuit source and assignment order are therefore part of an
experiment's specification.

The supported examples are:

| File | Parameters | Meaning |
| --- | --- | --- |
| `two_share.circ` | `n=2` | Minimal two-share GF(2) refresh. The array size is fixed at two. |
| `multi_refresh.circ` | `n >= 1`, `rounds >= 1` | GF(2) implementation of the paper's Algorithm R-multi; `rounds` is its parameter `k`. |
| `SNI_Refresh.circ` | `n >= 1` | GF(2) pairwise additive refresh. LAPSE's target-recovery check is not an SNI verification procedure. |
| `SNI_Refresh_field.circ` | `n >= 1` | The same pairwise refresh over GF(257), with opposite signs for the paired updates. |
| `Compress.circ` | `n >= 1` | The paper's linear Compress stage on a row-major array of `n*n` independent input symbols. |

For example:

```bash
python3 -m lapse --file multi_refresh.circ --n 2 --params rounds=1 --delta-p 0.01
python3 -m lapse --file SNI_Refresh_field.circ --n 3 --t 1
```

For R-multi, the explicit assignment `r[n-1] = 0` produces a deterministic wire.
For `n=2, rounds=1`, the paper lists eight nonconstant wires; this source also
records the zero initialization. Probing that zero wire supplies no information
and does not change the recovery probability. Equal-valued random, intermediate,
and final-output wires remain distinct probe locations, matching the paper's
example. Its exact recovery probability at `p=1/2` is `190/256`.

`Compress.circ` records the target sum of all matrix entries before masking,
masks each pair of off-diagonal entries with the same fresh random value, and
sums each row. These operations preserve the target. Treating its matrix entries
as independent source variables only models this linear subcircuit; it does not
prove security of nonlinear multiplications that might supply those entries.
