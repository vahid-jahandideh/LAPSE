# LAPSE

**Linear Algebraic Probing Security Estimator** analyzes linear circuits in the random probing model. It compiles a circuit into wire vectors over a finite field and computes the probability that leaked wires span a specified target `S`. The same exact enumeration, Monte Carlo (MC), and importance sampling (IS) methods apply to all supported circuits.

This target-recovery probability is denoted $E_X(n,p)$ in the paper, where $n$ is the number of shares in the paper's circuit families and $p$ is the independent per-wire probing probability.

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
lapse --file two_share.circ --n 2 --delta-p 0.1 --method exact

# Check whether any set of at most one probe recovers S
lapse --file two_share.circ --n 2 --t 1

# Importance sampling, with a fixed seed and statistical diagnostics
lapse --file multi_refresh.circ --n 4 --params rounds=3 \
  --delta-p 0.01 --method is --mc-samples 100000 --seed 2026 --json

# Parameter sweep; requires the experiments extra
lapse-sweep --file multi_refresh.circ --params rounds=1 \
  --n-list 2 3 --p-list 0.01 0.1 --method exact \
  --outfile results/example.csv --plots-out results/plots
```

The `--delta-p` option sets $p$ for computing $E_X(n,p)$. The recovery value appears as `recovery.estimate` in JSON output and as `delta` in sweep CSV files.

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

## Estimation and boundary exposure

`--method exact` exhaustively accounts for leakage subsets, using floating-point probabilities. It groups identical recovery vectors and integrates zero vectors; these reductions are valid for general linear circuits. The default `--max-enum 20` guards the original probeable wire count.

`--method mc` samples independent probes at rate `p`. `--method is` samples at a proposal rate and applies likelihood-ratio weights; set `--mc-q` to choose that rate, or let a discarded pilot tune it. The default `--method auto` uses exact enumeration for small circuits and a discarded pilot to select MC or IS otherwise. `--mc-samples` sets the main sample budget, separate from pilot samples.

JSON and sweep CSV output include sample and hit counts, standard errors, confidence intervals, proposal rates, effective sample sizes, and warnings. MC uses a Wilson interval; IS uses an approximate normal interval and provides no interval for zero hits. At extreme probabilities, a few dominant weights or missed events can make IS estimates and intervals unreliable. Zero hits do not establish zero probability, and a fixed seed establishes repeatability, not accuracy. Raw IS estimates and reliability warnings are retained. Sweeps also write a metadata sidecar with the circuit, parameters, versions, and seeds.

For sum-preserving refresh circuits, LAPSE also estimates conditional boundary exposure β under an explicit reduction policy. The default `separated` policy excludes public correctness consequences and conservatively exposes all boundary shares when cross-boundary dependence remains. This per-share quantity does not establish independence between shares. See the [β definition and limitations](docs/beta-definition.md).

```bash
lapse --file two_share.circ --n 2 --beta-p 0.1 \
  --boundary-share 'a[0]' --method exact --json
```

## Tests

```bash
python -m pip install -e '.[experiments]'
MPLBACKEND=Agg python -m unittest discover -s tests -v
```

Tests cover finite-field arithmetic, circuit semantics, algebraic probability oracles, estimator statistics, and command-line behavior. CI runs on Python 3.10–3.13.

