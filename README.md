# LAPSE

**Linear Algebraic Probing Security Estimator** analyzes linear circuits in the random probing model. It compiles a circuit into wire vectors over a finite field and computes the probability that leaked wires span a specified target `S`. The same exact enumeration, Monte Carlo (MC), and importance sampling (IS) methods apply to all supported circuits.

This target-recovery probability is denoted $E_X(n,p)$ in the paper, where $n$ is the number of shares in the paper's circuit families and $p$ is the independent per-wire probing probability. For refresh gadgets, LAPSE also computes the boundary exposure rate $\beta$ used in the paper's composition theorem (see [Boundary exposure β](#boundary-exposure-β)).

Repository: <https://github.com/vahid-jahandideh/LAPSE>

## Install

Python 3.10 or newer is required. From a source checkout:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

For CSV parameter sweeps and plots, install the optional dependencies with `python -m pip install -e '.[experiments]'`.

## Usage

```bash
# Exact target-recovery probability E_X(2, 0.1): 0.02152
lapse --file two_share.circ --n 2 --recovery-p 0.1 --method exact

# Check whether any set of at most one probe recovers S
lapse --file two_share.circ --n 2 --t 1

# Importance sampling, with a fixed seed and statistical diagnostics
lapse --file multi_refresh.circ --n 4 --params rounds=3 \
  --recovery-p 0.01 --method is --mc-samples 100000 --seed 2026 --json

# Parameter sweep; requires the experiments extra
lapse-sweep --file multi_refresh.circ --params rounds=1 \
  --n-list 2 3 --p-list 0.01 0.1 --method exact \
  --outfile results/example.csv --plots-out results/plots
```

The `--recovery-p` option sets the probing rate $p$ at which LAPSE computes $E_X(n,p)$. The recovery value appears as `recovery.estimate` in JSON output. In sweeps, `--metric delta` selects it and the CSV column `delta` holds it.

The module entry points are `python -m lapse` and `python -m lapse.experiments`. Use `--help` to see all options. The Python API accepts circuit source directly:

```python
from pathlib import Path
from lapse import build_ir, estimate_recovery, parse_program

source = Path("two_share.circ").read_text()
ir = build_ir(parse_program(source), {"n": 2})
result = estimate_recovery(ir, 0.1, method="exact")
print(result.estimate)  # E_X(2, 0.1)
```

## Circuit language and leakage model

```text
circuit two_share(n);
input gf2m(8, 0x11B) a[2];
rand gf2m(8, 0x11B) r;
output gf2m(8, 0x11B) b[2];
r = rand();
b[0] = a[0] + r;
b[1] = a[1] - r;
const gf2m(8, 0x11B) S;
S = sum_reduce(a, 0, 2);
```

Fields are `bit`/`gf2`, prime fields such as `gf(257)`, and binary extension fields such as `gf2m(8, 0x11B)`. One circuit uses one field; prime moduli and extension polynomials are validated. Expressions support addition, subtraction, characteristic-two `xor`, and multiplication by a numeric field literal. Nonlinear products are rejected.

Roles are `input`, `output`, `rand`, `temp`, and `const`. Arrays and loops accept integer parameters. Loop and reduction ranges are half-open: `for i in [0 to n) { b[i] = a[i]; }` and `sum_reduce(a, 0, n)` both visit indices `0` through `n-1`.

Each modeled wire leaks independently with probability `p`. Input declarations and every assignment create leakage sites, including copies, overwritten values, and constant-valued assignments. Internal operations within one expression are not separate sites; write separate assignments when they should leak independently. Each `rand()` call creates a fresh independent source. A bare `rand` variable is sampled on its first read, so explicit sampling is recommended. All versions of the target `S` are excluded from probes; other `const` values must be public. `--allow-roles` restricts the probeable roles.

Recovery means that the target vector belongs to the span of the leaked vectors. The `--t` check concerns this target only; it does not establish general SNI or composition security. See the [example circuits and their conventions](docs/circuits.md).

## Estimation

`--method exact` exhaustively accounts for leakage subsets, using floating-point probabilities. It groups identical recovery vectors and integrates zero vectors; these reductions are valid for general linear circuits. The default `--max-enum 20` guards the original probeable wire count.

`--method mc` samples independent probes at rate `p`. `--method is` samples at a proposal rate and applies likelihood-ratio weights; set `--mc-q` to choose that rate, or let a discarded pilot tune it. The default `--method auto` uses exact enumeration for small circuits and a discarded pilot to select MC or IS otherwise. `--mc-samples` sets the main sample budget, separate from pilot samples.

JSON and sweep CSV output include sample and hit counts, standard errors, confidence intervals, proposal rates, effective sample sizes, and warnings. MC uses a Wilson interval; IS uses an approximate normal interval and provides no interval for zero hits. At extreme probabilities, a few dominant weights or missed events can make IS estimates and intervals unreliable. Zero hits do not establish zero probability, and a fixed seed establishes repeatability, not accuracy. Raw IS estimates and reliability warnings are retained. Sweeps also write a metadata sidecar with the circuit, parameters, versions, and seeds.

## Boundary exposure β

Leakage inside a refresh gadget can give the attacker linear relations between its input shares and output shares. A boundary share that appears in such a relation is *affected*. For a boundary share $U$, LAPSE computes

$$\beta_U(p) = \Pr[\,U \text{ is affected} \mid U \text{ does not leak}\,], \qquad p' = p + (1-p)\,\beta_U(p).$$

Here $p'$ is the effective rate at which $U$ is exposed, directly or through a relation. These are the quantities $\beta(n,p)$ and $p'$ of the paper's composition theorem.

LAPSE decides which shares are affected as follows:

- Input-only and output-only relations affect the shares they involve.
- A relation that links inputs and outputs, and cannot be explained by one-sided relations and correctness, makes every boundary share affected.
- The public correctness relation $\sum_i X_i = \sum_i Y_i$ never counts as leakage.

```bash
# beta for input share a[0] at p = 0.1: prints beta = 0.0919 and p' = 0.18271
lapse --file two_share.circ --n 2 --beta-p 0.1 --boundary-share 'a[0]' --method exact
```

- **Share names.** `a[0]` selects the original input share `a[0]#0`, and `b[0]` selects the final output share. Full SSA names are also accepted.
- **Supported circuits.** The circuit must be a refresh: equal-size input and output sharings with equal sums, and no relation between them other than correctness when nothing leaks. `two_share.circ`, `multi_refresh.circ`, `SNI_Refresh.circ`, and `SNI_Refresh_field.circ` qualify; `Compress.circ` does not.
- **One value per share.** Shares of the same gadget can have different β. To obtain a single $p'$ for a refresh, compute β for every boundary share and use the largest.
- **Estimation.** `--method`, sample counts, and `--seed` work as for $E_X$, and `--beta-q` sets the importance-sampling proposal rate for β. JSON output reports `beta.estimate` and `beta.p_effective`.
- **Sweeps.** `lapse-sweep ... --metric beta --boundary-share 'a[0]'` writes a `beta` column, and `--metric recovery-beta` computes both quantities. In sweeps, `--mc-q` sets the proposal rate for both. The exponent printed by `--fit` is a curve-fit parameter, not the boundary exposure β.
- **Python.** `estimate_beta`, `boundary_shares`, and `inspect_boundary_relations` are exported by `lapse`.

See [Boundary exposure β](docs/beta-definition.md) for the precise definition, a worked example, and guidance on using β in composition bounds.
