"""Analytic small-circuit checks for exact, MC and importance-sampling recovery."""
import math
import random
import unittest

from lapse import core as rpt


def sharing_ir(n=2, zero_target=False):
    field = rpt.PrimeField(2)
    names = [f'a[{i}]#0' for i in range(n)]
    vectors = {name: [int(i == j) for j in range(n)] for i, name in enumerate(names)}
    vectors['S#1'] = [0 if zero_target else 1] * n
    return rpt.IR({}, {}, field, names, vectors,
                  [(name, 'input') for name in names] + [('S#1', 'const')], 'S#1')


class DrawSequence:
    def __init__(self, values):
        self.values = iter(values)
        self.draws = 0

    def random(self):
        self.draws += 1
        return next(self.values)


class EstimatorTests(unittest.TestCase):
    def test_exact_matches_all_shares_probability(self):
        for n in (1, 2, 3, 5):
            for p in (0.0, 0.01, 0.2, 0.5, 1.0):
                with self.subTest(n=n, p=p):
                    result = rpt.estimate_recovery(sharing_ir(n), p, method='exact')
                    self.assertAlmostEqual(result.estimate, p ** n)
                    self.assertEqual(result.wire_count, n)

    def test_exact_groups_duplicate_vectors_and_integrates_zero_sites(self):
        ir = sharing_ir()
        ir.created.extend([('b#0', 'output'), ('c#0', 'output'), ('z#0', 'temp')])
        ir.vectors.update({'b#0': [1, 0, 0], 'c#0': [0, 1], 'z#0': [0, 0]})
        for p in (1e-150, 0.01, 0.2, 0.9):
            result = rpt.estimate_recovery(ir, p, method='exact')
            expected = (p * (2 - p)) ** 2
            self.assertAlmostEqual(result.estimate / expected, 1.0)
            self.assertEqual(result.wire_count, 5)
        # The enumeration budget counts physical leakage sites.
        with self.assertRaises(ValueError):
            rpt.estimate_recovery(ir, 0.2, method='exact', max_enum=4)

    def test_public_target_includes_empty_pattern(self):
        ir = sharing_ir(2, zero_target=True)
        for method in ('exact', 'mc', 'is', 'auto'):
            result = rpt.estimate_recovery(ir, 0.0, method=method)
            self.assertEqual(result.estimate, 1.0)
        result = rpt.t_probing(ir, 0, None)
        self.assertFalse(result.ok)
        self.assertEqual(result.counterexample, [])

    def test_no_leakable_wires_and_missing_target(self):
        result = rpt.estimate_recovery(sharing_ir(), 0.5, set(), method='exact')
        self.assertEqual(result.estimate, 0.0)
        self.assertEqual(result.wire_count, 0)
        ir = sharing_ir()
        ir.secret_ssa = None
        with self.assertRaisesRegex(ValueError, 'No target'):
            rpt.estimate_recovery(ir, 0.2, method='mc', samples=5)

    def test_all_versions_of_target_are_excluded(self):
        ir = sharing_ir()
        ir.created.insert(0, ('S#0', 'const'))
        ir.vectors['S#0'] = [1, 1]
        self.assertEqual(rpt.filter_leakables(ir.created, ir.secret_ssa, {'const'}), [])

    def test_probability_and_budget_validation(self):
        ir = sharing_ir()
        for p in (-0.1, 1.1, math.inf, math.nan):
            with self.subTest(p=p), self.assertRaises(ValueError):
                rpt.estimate_recovery(ir, p, method='exact')
        for q in (0.0, 1.0, -0.1, math.inf, math.nan):
            with self.subTest(q=q), self.assertRaises(ValueError):
                rpt.estimate_recovery(ir, 0.2, method='is', samples=100, q=q)
        for samples in (0, -1, 1.5, True):
            with self.subTest(samples=samples), self.assertRaises(ValueError):
                rpt.estimate_recovery(ir, 0.2, method='mc', samples=samples)
        with self.assertRaisesRegex(ValueError, 'max_enum'):
            rpt.estimate_recovery(ir, 0.2, method='exact', max_enum=1)
        with self.assertRaises(ValueError):
            rpt.t_probing(ir, -1, None)
        with self.assertRaises(ValueError):
            rpt.filter_leakables(ir.created, ir.secret_ssa, {'typo'})
        with self.assertRaises(ValueError):
            rpt.estimate_recovery(ir, 0.2, seed=1, rng=random.Random(1))

    def test_mc_is_seed_reproducible(self):
        ir = sharing_ir()
        for method in ('mc', 'is', 'auto'):
            kwargs = dict(method=method, samples=200, seed=413, max_enum=0, pilot_samples=40)
            first = rpt.estimate_recovery(ir, 0.1, **kwargs)
            second = rpt.estimate_recovery(ir, 0.1, **kwargs)
            self.assertEqual(first, second)

    def test_mc_zero_hit_interval_is_not_zero_width(self):
        result = rpt.estimate_recovery(sharing_ir(), 0.01, method='mc', samples=10,
                                       rng=DrawSequence([0.9] * 20))
        self.assertEqual(result.estimate, 0.0)
        self.assertGreater(result.confidence_interval[1], 0.0)
        self.assertTrue(result.warnings)

    def test_is_weight_and_unbiased_sample_variance(self):
        # Exactly one of each of the four equiprobable patterns under q=1/2.
        draws = [0.9, 0.9, 0.1, 0.9, 0.9, 0.1, 0.1, 0.1]
        result = rpt.estimate_recovery(sharing_ir(), 0.1, method='is', q=0.5,
                                       samples=4, rng=DrawSequence(draws))
        self.assertAlmostEqual(result.estimate, 0.01)
        self.assertAlmostEqual(result.variance, 0.0001)
        self.assertAlmostEqual(result.standard_error, 0.01)
        self.assertAlmostEqual(result.effective_sample_size, 1.0)
        self.assertEqual(result.hits, 1)

    def test_is_zero_hit_interval_unavailable(self):
        result = rpt.estimate_recovery(sharing_ir(), 0.01, method='is', q=0.5,
                                       samples=10, rng=DrawSequence([0.9] * 20))
        self.assertEqual(result.estimate, 0.0)
        self.assertIsNone(result.confidence_interval)
        self.assertTrue(any('unresolved' in warning for warning in result.warnings))

    def test_is_does_not_clamp_target_rate_or_raw_estimate(self):
        result = rpt.estimate_recovery(sharing_ir(1), 1e-200, method='is', q=0.5,
                                       samples=2, rng=DrawSequence([0.1, 0.9]))
        self.assertAlmostEqual(result.estimate / 1e-200, 1.0)
        self.assertAlmostEqual(result.standard_error / 1e-200, 1.0)
        result = rpt.estimate_recovery(sharing_ir(), 0.9, method='is', q=0.1,
                                       samples=2, rng=DrawSequence([0.01, 0.01, 0.9, 0.9]))
        self.assertAlmostEqual(result.estimate, 40.5)

    def test_auto_discards_pilot_observations(self):
        # Pilot sees 10 hits and selects MC; fresh main observations have no hits.
        rng = DrawSequence([0.1] * 20 + [0.9] * 20)
        result = rpt.estimate_recovery(sharing_ir(), 0.5, method='auto', max_enum=0,
                                       samples=10, pilot_samples=10, rng=rng)
        self.assertEqual(result.method, 'mc')
        self.assertEqual(result.estimate, 0.0)
        self.assertEqual(result.pilot_samples, 10)
        self.assertEqual(rng.draws, 40)

    def test_sampling_below_input_subset_bound_warns(self):
        result = rpt.estimate_recovery(sharing_ir(), 0.5, method='mc', samples=10,
                                       rng=DrawSequence([0.9] * 20))
        self.assertTrue(any('input-subset' in warning for warning in result.warnings))

    def test_lower_bound_uses_original_inputs_only(self):
        ir = sharing_ir()
        ir.created.extend([('a[0]#1', 'input'), ('a[1]#1', 'input')])
        ir.vectors.update({'a[0]#1': [1, 0], 'a[1]#1': [0, 1]})
        # Estimate 0.1 is below p**2=0.25 but above the incorrect p**4=0.0625.
        draws = [0.1, 0.1, 0.9, 0.9] + [0.9] * 36
        result = rpt.estimate_recovery(ir, 0.5, method='mc', samples=10,
                                       rng=DrawSequence(draws))
        self.assertEqual(result.estimate, 0.1)
        self.assertTrue(any('bound 0.25' in warning for warning in result.warnings))

    def test_prepacked_recovery_matches_span_for_all_subsets(self):
        ir = sharing_ir(3)
        ir.vectors.update({'b#0': [1, 1], 'c#0': [0, 1, 1], 'd#0': [0]})
        ir.created.extend([('b#0', 'output'), ('c#0', 'output'), ('d#0', 'temp')])
        sites, event, _ = rpt._recovery_problem(ir, None)
        for mask in range(1 << len(sites)):
            indices = [i for i in range(len(sites)) if (mask >> i) & 1]
            expected = rpt.in_span([ir.vectors[sites[i]] for i in indices],
                                   ir.vectors[ir.secret_ssa], ir.field)
            self.assertEqual(event(indices), expected)

    def test_extension_field_recovery_uses_field_arithmetic(self):
        ir = sharing_ir(2)
        ir.field = rpt.GF2m(2, 0b111)
        ir.vectors = {'a[0]#0': [2, 1], 'a[1]#0': [0, 0], 'S#1': [1, 3]}
        # In GF(4), 3*(2,1)=(1,3), so leaking a[0] alone reveals S.
        result = rpt.estimate_recovery(ir, 0.3, method='exact')
        self.assertAlmostEqual(result.estimate, 0.3)

    def test_fixed_seed_mc_agrees_with_exact(self):
        result = rpt.estimate_recovery(sharing_ir(), 0.3, method='mc', samples=20000, seed=42)
        self.assertLess(abs(result.estimate - 0.09), 5 * result.standard_error)

    def test_generic_event_exact_matches_analytic_threshold(self):
        result = rpt.estimate_event(lambda subset: len(subset) >= 2, 3, 0.2,
                                    method='exact')
        self.assertAlmostEqual(result.estimate, 3 * 0.2 ** 2 * 0.8 + 0.2 ** 3)
        self.assertEqual(result.wire_count, 3)

    def test_generic_event_sampling_matches_recovery(self):
        for method in ('mc', 'is', 'auto'):
            arguments = dict(method=method, samples=1000, max_enum=0,
                             seed=42, q=0.4, pilot_samples=40)
            recovery = rpt.estimate_recovery(sharing_ir(), 0.2, **arguments)
            generic = rpt.estimate_event(lambda subset: len(subset) == 2, 2, 0.2,
                                         **arguments)
            self.assertEqual(generic.estimate, recovery.estimate)
            self.assertEqual(generic.variance, recovery.variance)
            self.assertEqual(generic.confidence_interval, recovery.confidence_interval)

    def test_generic_event_deterministic_and_invalid_count(self):
        self.assertEqual(rpt.estimate_event(lambda _: True, 0, 0.2).estimate, 1.0)
        self.assertEqual(rpt.estimate_event(lambda _: False, 0, 0.2).estimate, 0.0)
        self.assertEqual(rpt.estimate_event(lambda x: len(x) >= 1, 3, 1.0).estimate, 1.0)
        with self.assertRaises(ValueError):
            rpt.estimate_event(lambda _: False, -1, 0.2)

    def test_generic_nonmonotone_event_with_true_empty_pattern(self):
        event = lambda subset: len(subset) % 2 == 0
        result = rpt.estimate_event(event, 2, 0.25, method='exact')
        self.assertAlmostEqual(result.estimate, 0.625)
        self.assertEqual(rpt.estimate_event(event, 2, 0.0).estimate, 1.0)
        self.assertEqual(rpt.estimate_event(event, 1, 1.0).estimate, 0.0)
        # One of every equiprobable proposal pattern; the IS estimator sees both
        # the empty event and the complete pair, weighted at the target rate.
        draws = [0.9, 0.9, 0.1, 0.9, 0.9, 0.1, 0.1, 0.1]
        result = rpt.estimate_event(event, 2, 0.25, method='is', q=0.5,
                                   samples=4, rng=DrawSequence(draws))
        self.assertAlmostEqual(result.estimate, 0.625)


if __name__ == '__main__':
    unittest.main()
