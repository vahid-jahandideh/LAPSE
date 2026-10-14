"""Check general estimation against algebraic probability oracles across fields."""
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import lapse


ROOT = Path(__file__).resolve().parents[1]
FIELDS = [('bit', 1, 1), ('gf(3)', 2, 2),
          ('gf2m(2,0b111)', 2, 3), ('gf(257)', 17, 256)]


def scaled_copies_source(field, first, second, name='scaled_copies'):
    # Recovery occurs if b[2] leaks, or at least one copy of each input leaks.
    # Nonzero field scaling does not change the information in a share copy.
    return f'''circuit {name}(n);
input {field} a[n];
output {field} b[3];
const {field} S;
S = a[0] + a[1];
b[0] = {first} * a[0];
b[1] = {second} * a[1];
b[2] = a[0] + a[1];
'''


def recovery_probability(p):
    return p + (1 - p) * (p * (2 - p)) ** 2


class GeneralCircuitTests(unittest.TestCase):
    def test_exact_matches_independent_analytic_oracle_across_fields(self):
        for field, first, second in FIELDS:
            ir = lapse.build_ir(lapse.parse_program(
                scaled_copies_source(field, first, second)), {'n': 2})
            for p in (0, .01, .2, .5, 1):
                with self.subTest(field=field, p=p):
                    result = lapse.estimate_recovery(ir, p, method='exact')
                    self.assertEqual(result.wire_count, 5)
                    self.assertAlmostEqual(result.estimate, recovery_probability(p), places=14)

    def test_mc_and_is_match_analytic_probability_across_fields(self):
        p, samples = .2, 10000
        expected = recovery_probability(p)
        for field, first, second in FIELDS:
            ir = lapse.build_ir(lapse.parse_program(
                scaled_copies_source(field, first, second)), {'n': 2})
            for method, proposal in [('mc', p), ('is', .45)]:
                with self.subTest(field=field, method=method):
                    result = lapse.estimate_recovery(
                        ir, p, method=method, q=proposal if method == 'is' else None,
                        samples=samples, seed=6173)
                    # Analytic second moment of 1[event] * P(pattern)/Q(pattern).
                    # a/b are the leaked/unleaked contributions of one site;
                    # factorization follows the independent recovery alternatives.
                    a, b = p**2 / proposal, (1-p)**2 / (1-proposal)
                    second_moment = a * (a+b)**4 + b * ((a+b)**2 - b**2)**2
                    true_standard_error = math.sqrt((second_moment - expected**2) / samples)
                    self.assertEqual(result.method, method)
                    self.assertLess(abs(result.estimate - expected), 6 * true_standard_error)

    def test_circuit_name_does_not_select_an_estimator(self):
        for method in ('exact', 'auto', 'mc', 'is'):
            estimates = []
            for name in ('scaled_copies', 'multi_refresh'):
                ir = lapse.build_ir(lapse.parse_program(
                    scaled_copies_source('gf(3)', 2, 2, name)), {'n': 2})
                estimates.append(lapse.estimate_recovery(
                    ir, .1, method=method, max_enum=5 if method == 'exact' else 0,
                    samples=200, pilot_samples=40, seed=91))
            with self.subTest(method=method):
                self.assertEqual(estimates[0], estimates[1])
                self.assertIn(estimates[0].method, {'exact', 'mc', 'is'})

    def test_custom_circuit_sweep_accepts_more_than_three_rounds(self):
        source = '''circuit random_mixing(n, rounds);
input gf(3) a[n];
rand gf(3) r;
output gf(3) b[n];
const gf(3) S;
S = sum_reduce(a, 0, n);
for step in [0 to rounds) {
  r = rand();
  for i in [0 to n-1) {
    a[i] = a[i] + r;
    a[n-1] = a[n-1] - r;
  }
}
for i in [0 to n) { b[i] = a[i]; }
'''
        with tempfile.TemporaryDirectory() as directory:
            circuit = Path(directory)/'custom.circ'
            output = Path(directory)/'result.csv'
            circuit.write_text(source)
            result = subprocess.run(
                [sys.executable, '-m', 'lapse.experiments', '--file', str(circuit),
                 '--n-list', '3', '--p-list', '.1', '--params', 'rounds=4',
                 '--method', 'mc', '--samples', '2000', '--seed', '91',
                 '--outfile', str(output)], cwd=ROOT,
                text=True, capture_output=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)
            with output.open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]['method'], 'mc')
            self.assertEqual(int(rows[0]['samples']), 2000)
            self.assertGreater(int(rows[0]['wire_count']), 20)
            metadata = json.loads(output.with_suffix('.csv.metadata.json').read_text())
            self.assertEqual(metadata['parameters']['rounds'], 4)


if __name__ == '__main__':
    unittest.main()
