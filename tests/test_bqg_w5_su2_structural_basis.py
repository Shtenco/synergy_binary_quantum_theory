import unittest

from bqg_w5_su2_structural_basis import (
    allowed_pair2, all_total_J2, charged_labels, charged_node_basis, covariant_basis_size,
    gauss_labels, iter_covariant_keys, local_spins, validate_covariant_key,
    w5_gauss_seed, w5_seed_as_covariant, spin_support_for_sequence, t_sequence_spin_support,
)


class StructuralBasisTests(unittest.TestCase):
    def test_triangle_and_parity(self):
        self.assertEqual(allowed_pair2(4, 4), (0, 2, 4, 6, 8))
        self.assertEqual(allowed_pair2(1, 2), (1, 3))
        self.assertEqual(charged_labels((4, 4, 4, 4), 1), ())
        self.assertEqual(len(charged_labels((4, 4, 4, 4), 0)), 5)

    def test_complete_total_J_content(self):
        js = all_total_J2((4, 4, 4, 4))
        self.assertEqual(js, tuple(range(0, 17, 2)))
        self.assertEqual(len(gauss_labels((4, 4, 4, 4))), 5)
        self.assertEqual(len(charged_node_basis((4, 4, 4, 4))), 625)
        self.assertEqual(covariant_basis_size((4,) * 10, 0), 390625)

    def test_w5_seed_embedding(self):
        for K in (0, 2, 4, 6, 8):
            spins, gauss = w5_gauss_seed(K)
            self.assertEqual(spins, (4,) * 10)
            self.assertEqual(gauss, (K, 0, 0, 0, 0))
            key = w5_seed_as_covariant(K, 0)
            self.assertTrue(validate_covariant_key(key, 0))
            self.assertEqual(key[1][0], ("C", 0, 0, K, K))

    def test_lazy_basis_count_matches_iteration(self):
        spins = (1,) * 10
        factors_size = covariant_basis_size(spins, 0)
        # Small-spin sector is small enough for exhaustive enumeration.
        keys = list(iter_covariant_keys(spins, 0))
        self.assertEqual(len(keys), factors_size)
        self.assertEqual(len(set(keys)), factors_size)
        self.assertTrue(all(validate_covariant_key(k, 0) for k in keys))

    def test_word_spin_support_respects_edge_hits(self):
        spins = (4,) * 10
        # Two consecutive fundamental hits on edge (0,1).
        seq = (("P", (0, 1), 0, 0), ("P", (1, 0), 0, 0))
        got = spin_support_for_sequence(spins, seq, 10)
        edge_index = 0  # K5 edge ordering starts with (0,1).
        self.assertEqual({state[edge_index] for state in got}, {2, 4, 6})
        self.assertTrue(all(sum(s != 4 for s in state) <= 1 for state in got))

    def test_repository_t_sequence_support(self):
        support = t_sequence_spin_support((4,) * 10, 0, 1, 2, 3, 10)
        # The frozen T word hits edge (0,3) twice and three other edges once.
        self.assertEqual(len(support), 3 * (2 ** 3))
        self.assertEqual({sum(s != 4 for s in state) for state in support}, {3, 4})

    def test_source_node_local_spin_order(self):
        spins = tuple(range(10))
        self.assertEqual(local_spins(spins, 0), (0, 1, 2, 3))


if __name__ == "__main__":
    unittest.main()
