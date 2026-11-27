import unittest

from bqg_w5_su2_structural_basis import (
    allowed_pair2, all_total_J2, charged_labels, covariant_basis_size,
    gauss_labels, iter_covariant_keys, local_spins, validate_covariant_key,
    w5_gauss_seed, w5_seed_as_covariant,
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

    def test_source_node_local_spin_order(self):
        spins = tuple(range(10))
        self.assertEqual(local_spins(spins, 0), (0, 1, 2, 3))


if __name__ == "__main__":
    unittest.main()
