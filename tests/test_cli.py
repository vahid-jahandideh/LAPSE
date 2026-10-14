"""Public command-line and experiment regression tests."""
import contextlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import warnings

import numpy as np
import pandas as pd

from lapse import experiments as sweep
from lapse import core as rpt

ROOT = Path(__file__).resolve().parents[1]


class CLITests(unittest.TestCase):
    def run_cli(self, *args):
        return subprocess.run([sys.executable, '-m', 'lapse', *args], cwd=ROOT,
                              text=True, capture_output=True, timeout=30)

    def test_exact_json_reports_recovery(self):
        args = ['--demo_gf2', '--delta-p', '.1', '--method', 'exact', '--json']
        result = self.run_cli(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        data = json.loads(result.stdout)
        self.assertEqual(data['recovery']['method'], 'exact')
        self.assertAlmostEqual(data['recovery']['estimate'], .02152)

    def test_explicit_importance_sampling_reproducible(self):
        args = ['--demo_gf2', '--delta-p', '.01', '--mc-q', '.3',
                '--mc-samples', '1000', '--seed', '123', '--json']
        first, second = self.run_cli(*args), self.run_cli(*args)
        self.assertEqual(first.returncode, 0, first.stderr)
        self.assertEqual(first.stdout, second.stdout)
        self.assertEqual(json.loads(first.stdout)['recovery']['method'], 'is')

    def test_user_errors_are_actionable(self):
        for args in [('--demo_gf2', '--delta-p', 'nan'),
                     ('--demo_gf2', '--n', '0'),
                     ('--demo_gf2', '--params', 'n=2,n=3'),
                     ('--demo_gf2', '--allow-roles', 'typo'),
                     ('--demo_gf2', '--demo_gfp')]:
            with self.subTest(args=args):
                result = self.run_cli(*args)
                self.assertEqual(result.returncode, 2)
                self.assertNotIn('Traceback', result.stderr)

    def test_sweep_metadata_and_logspace_option(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'nested' / 'grid.csv'
            with contextlib.redirect_stdout(io.StringIO()):
                sweep.main(['--file', str(ROOT/'two_share.circ'), '--n-list', '2',
                            '--p-logspace', '-2', '-1', '2', '--method', 'exact',
                            '--outfile', str(path)])
            data = pd.read_csv(path)
            self.assertEqual(len(data), 2)
            self.assertAlmostEqual(data.p.iloc[0], .01)
            metadata = json.loads(path.with_suffix('.csv.metadata.json').read_text())
            self.assertIn('circuit_sha256', metadata)
            self.assertEqual(metadata['options']['seed'], 0)
            self.assertIn('standard_error', data)

    def test_grid_seed_does_not_depend_on_order(self):
        cfg = sweep.EstimatorCfg(samples=100, method='mc', seed=42)
        first = sweep.grid_values([2], [.1, .2], cfg, rpt.DEMO_GF2)
        second = sweep.grid_values([2], [.2, .1], cfg, rpt.DEMO_GF2)
        pd.testing.assert_frame_equal(first.sort_values('p').reset_index(drop=True),
                                      second.sort_values('p').reset_index(drop=True))

    def test_fitting_recovers_known_model_without_zero_imputation(self):
        df = pd.DataFrame([{'n': n, 'p': p, 'delta': (2*p)**(.75*n)}
                           for p in [.01, .03, .1] for n in [2, 3, 4]])
        df = pd.concat([df, pd.DataFrame([{'n': 5, 'p': .01, 'delta': 0}])])
        with warnings.catch_warnings(record=True) as caught:
            fit = sweep.fit_alpha_beta(df, 'delta')
        self.assertTrue(caught)
        self.assertAlmostEqual(fit.alpha0, 2)
        self.assertAlmostEqual(fit.beta, .75)
        with self.assertRaises(ValueError):
            sweep.fit_alpha_beta(df[df.n == 2], 'delta')

    def test_plots_support_singleton_and_unresolved_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            sweep.plot_curves(pd.DataFrame([{'n': 2, 'p': 0, 'delta': 0}]), 'delta', directory)
            self.assertEqual(len(list(Path(directory).glob('*.png'))), 3)


if __name__ == '__main__':
    unittest.main()

class BetaCLITests(unittest.TestCase):
    def test_beta_cli_exact_and_policy_provenance(self):
        result = subprocess.run([sys.executable, '-m', 'lapse', '--file', 'two_share.circ',
                                 '--n', '2', '--beta-p', '.1', '--boundary-share', 'a[0]',
                                 '--beta-policy', 'separated', '--method', 'exact', '--json'],
                                cwd=ROOT, capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        beta = json.loads(result.stdout)['beta']
        self.assertEqual(beta['policy'], 'separated')
        self.assertEqual(beta['boundary_share'], 'a[0]#0')
        self.assertAlmostEqual(beta['estimate'], .0919)
        self.assertAlmostEqual(beta['p_effective'], .18271)

    def test_beta_sweep(self):
        cfg = sweep.EstimatorCfg(method='exact', boundary_share='a[0]')
        frame = sweep.grid_values([2], [.1], cfg, rpt.DEMO_GF2, do_delta=True, do_beta=True)
        self.assertAlmostEqual(frame.beta.iloc[0], .0919)
        self.assertAlmostEqual(frame.delta.iloc[0], .02152)
        self.assertEqual(frame.beta_policy.iloc[0], 'separated')
