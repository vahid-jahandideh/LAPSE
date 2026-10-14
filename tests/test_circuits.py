"""Algebraic checks for circuits included in the public release."""

from pathlib import Path
import unittest

from lapse import core as lapse


ROOT = Path(__file__).resolve().parents[1]


def build(name, n, **params):
    text = (ROOT / name).read_text(encoding="utf-8")
    return lapse.build_ir(lapse.parse_program(text), {"n": n, **params})


def vector_sum(ir, names):
    result = [ir.field.zero()] * len(ir.sources)
    for name in names:
        row = ir.vectors[name]
        for index, coefficient in enumerate(row):
            result[index] = ir.field.add(result[index], coefficient)
    return result


class CircuitTests(unittest.TestCase):
    def assert_conservation(self, ir, output_count):
        inputs, outputs = lapse.candidates_inputs_outputs(ir)
        self.assertEqual(len(outputs), output_count)
        native = vector_sum(ir, inputs)
        self.assertEqual(vector_sum(ir, outputs), native)
        self.assertEqual(vector_sum(ir, [ir.secret_ssa]), native)

    def test_supported_refreshes_preserve_native(self):
        for n in (1, 2, 3, 4):
            for name, params in (
                ("SNI_Refresh.circ", {}),
                ("SNI_Refresh_field.circ", {}),
                ("multi_refresh.circ", {"rounds": 1}),
                ("multi_refresh.circ", {"rounds": 3}),
            ):
                with self.subTest(circuit=name, n=n, **params):
                    self.assert_conservation(build(name, n, **params), n)

    def test_two_share_example(self):
        self.assert_conservation(build("two_share.circ", 2), 2)

    def test_compression_preserves_all_matrix_entries(self):
        for n in (1, 2, 3, 4):
            with self.subTest(n=n):
                ir = build("Compress.circ", n)
                self.assert_conservation(ir, n)
                inputs, _ = lapse.candidates_inputs_outputs(ir)
                self.assertEqual(len(inputs), n * n)

    def test_paper_two_share_one_round_wire_vectors(self):
        ir = build("multi_refresh.circ", 2, rounds=1)
        leakables = lapse.filter_leakables(ir.created, ir.secret_ssa, None)
        nonzero = [ir.vectors[name] for name in leakables if any(ir.vectors[name])]
        # The paper omits the harmless deterministic R_n = 0 initialization.
        self.assertEqual(len(nonzero), 8)
        source_columns = [i for i, source in enumerate(ir.sources) if source != "CONST_ONE"]
        self.assertEqual(
            [[row[i] for i in source_columns] for row in nonzero],
            [[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 1],
             [1, 0, 1], [0, 1, 1], [1, 0, 1], [0, 1, 1]],
        )
        # Of the 2^8 equally likely subsets of these eight rows, 190 span X.
        result = lapse.estimate_recovery(ir, 0.5, method='exact')
        self.assertEqual(result.estimate, 190 / 256)


if __name__ == "__main__":
    unittest.main()
