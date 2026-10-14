"""Regression tests for the DSL compiler and exact finite-field algebra."""

import itertools
import random
import unittest
from unittest.mock import patch

from lapse import core as lapse


def compile_circuit(body, params=None):
    return lapse.build_ir(lapse.parse_program("circuit test(n);\n" + body), params or {})


class FieldTests(unittest.TestCase):
    def test_reject_composite_moduli(self):
        for modulus in (0, 1, 4, 9, 341, 561, 3215031751, 341550071728321):
            with self.subTest(modulus=modulus), self.assertRaises(ValueError):
                lapse.PrimeField(modulus)

    def test_prime_field_inverses(self):
        for modulus in (2, 3, 17, 257):
            field = lapse.PrimeField(modulus)
            for value in range(1, modulus):
                self.assertEqual(field.mul(value, field.inv(value)), 1)

    def test_reject_reducible_polynomial(self):
        for degree, polynomial in ((2, 0b101), (3, 0b1000), (4, 0b10101)):
            with self.subTest(polynomial=polynomial), self.assertRaises(ValueError):
                lapse.GF2m(degree, polynomial)

    def test_extension_field_inverses(self):
        for degree, polynomial in ((1, 3), (2, 7), (3, 11), (8, 0x11B)):
            field = lapse.GF2m(degree, polynomial)
            for value in range(1, 1 << degree):
                self.assertEqual(field.mul(value, field.inv(value)), 1)
            with self.assertRaises(ZeroDivisionError):
                field.inv(0)

    def test_aes_product(self):
        self.assertEqual(lapse.GF2m(8, 0x11B).mul(0x57, 0x83), 0xC1)

    def test_negative_extension_literal_is_additive_inverse(self):
        field = lapse.GF2m(8, 0x11B)
        self.assertEqual(field.from_int(-1), field.neg(field.from_int(1)))
        self.assertEqual(field.from_int(-2), field.neg(field.from_int(2)))


class LinearAlgebraTests(unittest.TestCase):
    def test_extension_field_is_not_binary_elimination(self):
        field = lapse.GF2m(2, 7)
        self.assertFalse(lapse.in_span([[2, 1]], [0, 1], field))
        self.assertTrue(lapse.in_span([[2, 1]], [1, 3], field))

    def test_span_against_exhaustive_linear_combinations(self):
        rng = random.Random(73)
        for field, order in ((lapse.PrimeField(2), 2), (lapse.PrimeField(3), 3),
                             (lapse.GF2m(2, 7), 4)):
            for _ in range(20):
                rows = [[rng.randrange(order) for _ in range(3)] for _ in range(2)]
                span = set()
                for coefficients in itertools.product(range(order), repeat=len(rows)):
                    value = [0] * 3
                    for coefficient, row in zip(coefficients, rows):
                        value = [field.add(a, field.mul(coefficient, b))
                                 for a, b in zip(value, row)]
                    span.add(tuple(value))
                expected_rank = 0
                while order ** expected_rank < len(span):
                    expected_rank += 1
                self.assertEqual(lapse.rank_over_field(rows, field), expected_rank)
                for target in itertools.product(range(order), repeat=3):
                    self.assertEqual(lapse.in_span(rows, list(target), field), target in span)

    def test_numpy_and_python_packed_paths_agree(self):
        rng = random.Random(91)
        for width in (0, 1, 63, 64, 65, 129):
            rows = [[rng.randrange(2) for _ in range(width)] for _ in range(4)]
            target = [rng.randrange(2) for _ in range(width)]
            expected = lapse._in_rowspan_gf2_packed(rows, target)
            with patch.object(lapse, "_np", None):
                self.assertEqual(lapse._in_rowspan_gf2_packed(rows, target), expected)

    def test_prepacked_subset_event_matches_span(self):
        rng = random.Random(274)
        for width in (0, 1, 5, 64, 129):
            rows = [[rng.randrange(2) for _ in range(width)] for _ in range(5)]
            targets = [[0] * width, rows[0], [rng.randrange(2) for _ in range(width)]]
            for target in targets:
                event = lapse._make_gf2_span_event(rows, target)
                for mask in range(1 << len(rows)):
                    indices = [i for i in range(len(rows)) if mask >> i & 1]
                    self.assertEqual(event(indices),
                                     lapse._in_rowspan_gf2_packed([rows[i] for i in indices], target))

    def test_zero_padding_and_empty_span(self):
        for field in (lapse.PrimeField(2), lapse.PrimeField(257), lapse.GF2m(2, 7)):
            self.assertTrue(lapse.in_span([], [0, 0], field))
            self.assertFalse(lapse.in_span([], [0, 1], field))
            self.assertTrue(lapse.in_span([[1]], [1, 0, 0], field))
            self.assertFalse(lapse.in_span([[1]], [0, 1], field))


class CompilerTests(unittest.TestCase):
    def test_all_builtin_demos(self):
        for name in ("DEMO_GF2", "DEMO_GFP", "DEMO_GF2M"):
            with self.subTest(demo=name):
                ir = lapse.build_ir(lapse.parse_program(getattr(lapse, name)), {"n": 2})
                self.assertIsNotNone(ir.secret_ssa)
                self.assertTrue(all(len(row) == len(ir.sources) for row in ir.vectors.values()))

    def test_initializer_runs_once_and_marks_secret(self):
        ir = compile_circuit("input bit a; rand bit r = rand(); const bit S = a xor r;")
        self.assertEqual([ssa for ssa, role in ir.created if role == "rand"], ["r#0"])
        self.assertEqual(ir.sources, ["a#0", "__randcall_0"])
        self.assertEqual(ir.secret_ssa, "S#0")
        self.assertEqual(ir.vectors[ir.secret_ssa], [1, 1])

    def test_initializer_follows_preceding_assignment(self):
        ir = compile_circuit("input bit a; temp bit x; x = a; const bit S = x;")
        self.assertEqual(ir.vectors[ir.secret_ssa], [1])

    def test_public_offsets_do_not_mask_input(self):
        ir = compile_circuit("input gf(257) a; output gf(257) b = a + 17; const gf(257) S = a;")
        self.assertEqual(ir.sources, ["a#0"])
        self.assertTrue(lapse.in_span([ir.vectors["b#0"]], ir.vectors[ir.secret_ssa], ir.field))

    def test_subtraction_preserves_new_randomness_column(self):
        ir = compile_circuit("input gf(257) a; rand gf(257) r; output gf(257) b = a - r;")
        self.assertEqual(ir.vectors["b#0"], [1, 256])

    def test_invalid_accesses_and_declarations_fail(self):
        bodies = (
            "input bit a[2]; output bit b = a[2];",
            "input bit a[2]; output bit b = a[-1];",
            "input bit a[2]; output bit b = a;",
            "input bit a; output bit b = a[0];",
            "output bit b; b[0] = 1;",
            "output bit b[2]; b[2] = 1;",
            "temp bit a; output bit b = a;",
            "const bit c; output bit b = c;",
            "input bit a; output gf(257) b = a;",
            "input bit a[2] = 1;",
            "input bit a = 0;",
            "input bit S;",
            "const bit S[2];",
            "rand bit __r;",
            "const bit c = 0; c = 1;",
            "input bit a; const bit c = a;",
            "for i in [0 .. 2) { input bit a; }",
            "input bit a[pow2(-1)];",
            "rand bit r = rand(1);",
            "input gf(3) a[2]; const gf(3) S = xor_reduce(a, 0, 2);",
            "input bit a[2]; const bit S = sum_reduce(a, 0, 3);",
        )
        for body in bodies:
            with self.subTest(body=body), self.assertRaises(lapse.EvalError):
                compile_circuit(body)

    def test_loop_local_randomness_is_fresh(self):
        ir = compile_circuit("""
            input bit a[2]; output bit b[2];
            for i in [0 .. 2) { rand bit r; b[i] = a[i] xor r; }
            const bit S = sum_reduce(a, 0, 2);
        """)
        self.assertEqual(ir.sources, ["a[0]#0", "a[1]#0", "r#0", "r#1"])
        self.assertEqual(ir.vectors["b[0]#0"], [1, 0, 1, 0])
        self.assertEqual(ir.vectors["b[1]#0"], [0, 1, 0, 1])

    def test_nested_loop_index_restored(self):
        ir = compile_circuit("""
            input bit a[2]; output bit b[2]; temp bit t;
            for i in [0 .. 2) {
                for i in [0 .. 1) { t = a[i]; }
                b[i] = a[i];
            }
        """, {"i": 9})
        self.assertEqual(ir.vectors["b[1]#0"], [0, 1])

    def test_reverse_loop_and_integer_literals(self):
        ir = compile_circuit("""
            input gf2m(02, 0b111) a[2]; output gf2m(2, 0x7) b[2];
            for i in [1 to -1 step -1) { b[i] = 0x02 * a[i]; }
        """)
        self.assertEqual(ir.vectors["b[0]#0"], [2, 0])
        self.assertEqual(ir.vectors["b[1]#0"], [0, 2])

    def test_refresh_range_preserves_secret_and_ssa_history(self):
        for field in ("bit", "gf(257)", "gf2m(2, 7)"):
            with self.subTest(field=field):
                ir = compile_circuit(f"""
                    input {field} a[3]; refresh_range(a, 0, 3);
                    const {field} S = sum_reduce(a, 0, 3);
                """)
                self.assertEqual(ir.vectors[ir.secret_ssa], [1, 1, 1] + [0] * (len(ir.sources) - 3))
                for i in range(3):
                    expected = [0] * len(ir.sources)
                    expected[i] = 1
                    self.assertEqual(ir.vectors[f"a[{i}]#0"], expected)
                self.assertTrue(any(role == "rand" for _, role in ir.created))

    def test_invalid_refresh_arguments(self):
        for call in ("refresh_range(a, -1, 2)", "refresh_range(a, 0, 4)"):
            with self.subTest(call=call), self.assertRaises(lapse.EvalError):
                compile_circuit(f"input bit a[3]; {call};")


if __name__ == "__main__":
    unittest.main()
