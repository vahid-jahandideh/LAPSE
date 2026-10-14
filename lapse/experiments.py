#!/usr/bin/env python3
"""Reproducible LAPSE sweeps, uncertainty output, and exploratory fits."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import sys
import warnings
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lark.exceptions import LarkError

from lapse import core as rpt
from lapse.cli_utils import clean_json, finite_probability, parse_params, parse_roles, provenance, result_dict


@dataclass
class EstimatorCfg:
    samples: int = 200_000
    max_enum: int = 20
    mc_q: Optional[float] = None
    allow_roles: Optional[list[str]] = None
    seed: int = 0
    method: str = 'auto'
    pilot_samples: int = 1000
    boundary_share: Optional[str] = None
    beta_policy: str = 'separated'


def build_ir_for_n(n, circuit_text, extra_params=None):
    if n < 1:
        raise ValueError('n must be positive')
    params = dict(extra_params or {})
    if 'n' in params:
        raise ValueError('use the n grid, not --params n=...')
    params['n'] = n
    return rpt.build_ir(rpt.parse_program(circuit_text), params)


def _point_seed(seed, n, p, metric):
    # Independent of grid ordering, and stable across processes/Python hash seeds.
    key = f'{seed}:{n}:{float(p).hex()}:{metric}'.encode('ascii')
    return int.from_bytes(hashlib.sha256(key).digest()[:8], 'big')


def _recovery(ir, p, n, cfg):
    return rpt.estimate_recovery(
        ir, p, set(cfg.allow_roles) if cfg.allow_roles is not None else None,
        samples=cfg.samples, max_enum=cfg.max_enum, q=cfg.mc_q,
        method='is' if cfg.mc_q is not None else cfg.method,
        seed=_point_seed(cfg.seed, n, p, 'delta'), pilot_samples=cfg.pilot_samples)


def grid_values(n_values, p_values, cfg, circuit_text, do_delta=True, extra_params=None, do_beta=False):
    if not n_values or not p_values:
        raise ValueError('parameter grids must be nonempty')
    if cfg.samples < 1 or cfg.max_enum < 0 or cfg.pilot_samples < 1:
        raise ValueError('sample counts must be positive and max-enum nonnegative')
    rows = []
    for n in n_values:
        ir = build_ir_for_n(n, circuit_text, extra_params)
        for p in p_values:
            finite_probability(p)
            rec = {'n': n, 'p': p}
            if do_delta:
                result = result_dict(_recovery(ir, p, n, cfg))
                rec['delta'] = result.pop('estimate')
                ci = result.pop('confidence_interval')
                rec['ci95_low'], rec['ci95_high'] = ci if ci is not None else (None, None)
                rec['warnings'] = '; '.join(result.pop('warnings'))
                rec.update(result)
            if do_beta:
                from lapse.boundary import estimate_beta
                if not cfg.boundary_share:
                    raise ValueError('beta requires --boundary-share')
                result = result_dict(estimate_beta(
                    ir, p, cfg.boundary_share, policy=cfg.beta_policy,
                    allow_roles=set(cfg.allow_roles) if cfg.allow_roles is not None else None,
                    method='is' if cfg.mc_q is not None else cfg.method,
                    samples=cfg.samples, max_enum=cfg.max_enum, q=cfg.mc_q,
                    seed=_point_seed(cfg.seed, n, p, 'beta'), pilot_samples=cfg.pilot_samples))
                rec['beta'] = result.pop('estimate')
                ci = result.pop('confidence_interval')
                rec['beta_ci95_low'], rec['beta_ci95_high'] = ci if ci is not None else (None, None)
                result['warnings'] = '; '.join(result['warnings'])
                rec.update({'beta_'+key: value for key, value in result.items()})
            rows.append(rec)
    return pd.DataFrame(rows)


@dataclass
class FitResult:
    alpha0: float
    beta: float
    r2_stage1: float
    r2_stage2: float


def _linear_fit(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 2 or len(np.unique(x)) < 2:
        raise ValueError('a fit requires at least two distinct x values')
    slope, intercept = np.linalg.lstsq(np.column_stack([x, np.ones_like(x)]), y, rcond=None)[0]
    residual = np.sum((y - slope*x - intercept)**2)
    total = np.sum((y - y.mean())**2)
    r2 = 1 - residual/total if total > 0 else float('nan')
    return float(slope), float(intercept), float(r2)


def fit_alpha_beta(df, value_col):
    """Exploratory two-stage fit. Zero/no-hit estimates are never imputed.

    The fitted exponent beta is unrelated to the paper's boundary exposure beta.
    Excluding unresolved points can still introduce selection bias.
    """
    valid = np.isfinite(df[value_col]) & (df[value_col] > 0) & np.isfinite(df.p) & (df.p > 0)
    if (~valid).any():
        warnings.warn('Excluded zero/nonfinite estimates from fit; no-hit estimates are unresolved, not tiny positive values.', stacklevel=2)
    slopes, logps, scores = [], [], []
    for p, group in df[valid].groupby('p'):
        if group.n.nunique() < 2:
            continue
        slope, _, score = _linear_fit(group.n, np.log(group[value_col]))
        slopes.append(slope)
        logps.append(math.log(p))
        scores.append(score)
    if len(slopes) < 2:
        raise ValueError('fit requires at least two positive p values, each with two distinct n and positive estimates')
    beta, intercept, score = _linear_fit(logps, slopes)
    exponent = intercept / beta if abs(beta) > 1e-14 else float('nan')
    alpha = math.exp(exponent) if math.isfinite(exponent) and exponent < 709 else float('nan')
    finite_scores = [x for x in scores if math.isfinite(x)]
    return FitResult(alpha, beta, float(np.mean(finite_scores)) if finite_scores else float('nan'), score)


def plot_curves(df, value_col, outdir=None):
    """Draw gaps for unresolved zeros rather than fabricated 1e-300 values."""
    path = Path(outdir) if outdir else None
    if path:
        path.mkdir(parents=True, exist_ok=True)
    data = df.copy()
    data.loc[data[value_col] <= 0, value_col] = np.nan
    for x, groupcol, suffix in [('n', 'p', 'vs_n'), ('p', 'n', 'vs_p')]:
        fig, ax = plt.subplots()
        for group, values in data.groupby(groupcol):
            values = values.sort_values(x)
            if x == 'p':
                values = values[values.p > 0]
            ax.plot(values[x], values[value_col], marker='.', label=f'{groupcol}={group:g}')
        if data[value_col].notna().any():
            ax.set_yscale('log')
        if x == 'p' and (data.p > 0).any():
            ax.set_xscale('log')
        ax.set(xlabel=x, ylabel=value_col, title=f'{value_col} vs {x} (unresolved zeros omitted)')
        ax.legend()
        fig.tight_layout()
        if path:
            fig.savefig(path / f'{value_col}_{suffix}.png', dpi=160)
        else:
            plt.show()
        plt.close(fig)
    pivot = data.pivot(index='n', columns='p', values=value_col)
    fig, ax = plt.subplots()
    values = np.ma.masked_invalid(np.log10(pivot.to_numpy()))
    im = ax.imshow(values, aspect='auto', origin='lower')
    ax.set_xticks(range(len(pivot.columns)), [f'{p:g}' for p in pivot.columns], rotation=45)
    ax.set_yticks(range(len(pivot.index)), [str(n) for n in pivot.index])
    ax.set(xlabel='p', ylabel='n', title=f'log10 {value_col} (unresolved zeros omitted)')
    fig.colorbar(im, ax=ax, label=f'log10 {value_col}')
    fig.tight_layout()
    if path:
        fig.savefig(path / f'{value_col}_heatmap.png', dpi=160)
    else:
        plt.show()
    plt.close(fig)


def _inclusive_range(start, stop, step):
    if not all(math.isfinite(v) for v in (start, stop, step)) or step <= 0 or stop < start:
        raise ValueError('ranges require finite START <= STOP and STEP > 0')
    count = int(math.floor((stop-start)/step + 1e-10)) + 1
    if count > 100000:
        raise ValueError('range has more than 100000 points; use a coarser grid')
    return [start + i*step for i in range(count)]


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description='LAPSE recovery sweeps with seed, uncertainty, and provenance')
    ap.add_argument('--file', required=True, help='Circuit source (.circ)')
    ap.add_argument('--params', default='')
    ap.add_argument('--metric', choices=['delta', 'beta', 'recovery-beta'], default='delta')
    ap.add_argument('--boundary-share', help='Boundary share for beta, e.g. a[0]')
    ap.add_argument('--beta-policy', choices=['separated', 'literal'], default='separated')
    ns = ap.add_mutually_exclusive_group()
    ns.add_argument('--n-list', type=int, nargs='+')
    ns.add_argument('--n-range', type=int, nargs=3, metavar=('START', 'STOP', 'STEP'))
    ps = ap.add_mutually_exclusive_group()
    ps.add_argument('--p-list', type=float, nargs='+')
    ps.add_argument('--p-range', type=float, nargs=3, metavar=('START', 'STOP', 'STEP'))
    ps.add_argument('--p-logspace', type=float, nargs=3, metavar=('START_EXP', 'STOP_EXP', 'COUNT'))
    ap.add_argument('--samples', type=int, default=200000)
    ap.add_argument('--max-enum', type=int, default=20)
    ap.add_argument('--pilot-samples', type=int, default=1000)
    ap.add_argument('--method', choices=['auto', 'exact', 'mc', 'is'], default='auto')
    ap.add_argument('--mc-q', type=float)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--allow-roles')
    ap.add_argument('--outfile')
    ap.add_argument('--plots-out')
    ap.add_argument('--fit', action='store_true', help='Exploratory exponent fit (not boundary beta)')
    return ap, ap.parse_args(argv)


def main(argv=None):
    ap, args = parse_args(argv)
    try:
        text = Path(args.file).read_text(encoding='utf-8')
        params = parse_params(args.params)
        roles = parse_roles(args.allow_roles)
        ns = args.n_list or ([int(n) for n in _inclusive_range(*args.n_range)] if args.n_range else list(range(2, 9)))
        if args.p_list:
            ps = args.p_list
        elif args.p_range:
            ps = _inclusive_range(*args.p_range)
        elif args.p_logspace:
            start, stop, count = args.p_logspace
            if not all(math.isfinite(x) for x in (start, stop, count)) or count != int(count) or not 1 <= count <= 100000:
                raise ValueError('p-logspace COUNT must be a positive integer <= 100000')
            ps = np.logspace(start, stop, int(count)).tolist()
        else:
            ps = np.logspace(-4, -1, 7).tolist()
        ns, ps = sorted(set(ns)), sorted(set(ps))
        if args.mc_q is not None and args.method not in ('auto', 'is'):
            raise ValueError('--mc-q requires --method auto or is')
        cfg = EstimatorCfg(args.samples, args.max_enum, args.mc_q,
                           sorted(roles) if roles is not None else None,
                           args.seed, args.method, args.pilot_samples, args.boundary_share, args.beta_policy)
        df = grid_values(ns, ps, cfg, text,
                         do_delta=args.metric in ('delta', 'recovery-beta'),
                         extra_params=params, do_beta=args.metric in ('beta', 'recovery-beta'))
        if args.outfile:
            path = Path(args.outfile)
            path.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(path, index=False)
            metadata = provenance(text, params, vars(args))
            metadata.update({'n_grid': ns, 'p_grid': ps, 'circuit_source': text,
                             'source_sha256': {Path(m.__file__).name: hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest()
                                               for m in [rpt, sys.modules[__name__]]}})
            path.with_suffix(path.suffix + '.metadata.json').write_text(
                json.dumps(clean_json(metadata), indent=2, allow_nan=False) + '\n', encoding='utf-8')
            print(f'Wrote {path} and {path.name}.metadata.json ({len(df)} rows)')
        else:
            print(df.to_string(index=False))
        for col in ('delta', 'beta'):
            if col not in df:
                continue
            if args.plots_out:
                plot_curves(df, col, args.plots_out)
            if args.fit:
                fit = fit_alpha_beta(df, col)
                print(f'Exploratory {col} fit: alpha0={fit.alpha0:.6g}, exponent_beta={fit.beta:.6g}, '
                      f'R2(stage1)={fit.r2_stage1:.4f}, R2(stage2)={fit.r2_stage2:.4f}')
    except (ValueError, OSError, rpt.EvalError, LarkError) as exc:
        ap.error(str(exc))


if __name__ == '__main__':
    main()
