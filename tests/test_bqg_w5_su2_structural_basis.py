import unittest

from bqg_w5_su2_structural_basis import (
    allowed_pair2, all_total_J2, charged_labels, charged_node_basis, covariant_basis_size,
    gauss_labels, iter_covariant_keys, local_spins, validate_covariant_key,
    w5_gauss_seed, w5_seed_as_covariant, spin_support_for_sequence, t_sequence_spin_support, spin_sector_basis_profile,
    charged_nodes_basis_size, charged_nodes_basis_profile, volume_label_targets, lift_key_to_all_charged,
    project_all_charged_key, advance_hit_label_support, structural_sequence_support,
    C_volume_structural_support,
)


class StructuralBasisTests(unittest.TestCase):
    # CI exercises the real repository word templates and structural label rules.
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
        self.assertEqual(charged_nodes_basis_size((4,) * 10, (0, 1)), 48828125)

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
        profile = spin_sector_basis_profile(support, 0)
        sizes = [n for _, n in profile]
        self.assertEqual(len(profile), 24)
        self.assertEqual((min(sizes), max(sizes), sum(sizes)), (57600, 630000, 5759800))
        two_charge = charged_nodes_basis_profile(support, (0, 1))
        two_sizes = [n for _, n in two_charge]
        self.assertEqual((min(two_sizes), max(two_sizes), sum(two_sizes)), (5760000, 113400000, 878540000))

    def test_volume_selection_preserves_total_J_and_M(self):
        local = (4, 4, 4, 4)
        label = ("C", 4, 0, 0, 4)
        targets = volume_label_targets(local, label)
        self.assertTrue(targets)
        self.assertTrue(all(t[1:3] == (4, 0) for t in targets))
        self.assertEqual({(t[3], t[4]) for t in targets}, set(charged_labels(local, 4)))

    def test_gauss_lift_and_final_projection_roundtrip(self):
        key = ((0,) * 10, tuple(("G", 0) for _ in range(5)))
        lifted = lift_key_to_all_charged(key)
        self.assertTrue(all(lab == ("C", 0, 0, 0, 0) for lab in lifted[1]))
        projected = project_all_charged_key(lifted, (0, 1))
        self.assertEqual(projected, (key[0], (("C", 0, 0, 0, 0), ("C", 0, 0, 0, 0), ("G", 0), ("G", 0), ("G", 0))))

    def test_hit_support_updates_only_endpoint_labels_and_link(self):
        key = ((0,) * 10, tuple(("G", 0) for _ in range(5)))
        states = list(advance_hit_label_support(key, 0, 1, 2))
        self.assertEqual(len(states), 4)
        self.assertTrue(all(s[0][0] == 1 for s in states))
        self.assertTrue(all(s[0][1:] == (0,) * 9 for s in states))
        self.assertTrue(all(s[1][2:] == tuple(("C", 0, 0, 0, 0) for _ in range(3)) for s in states))

    def test_complete_sequence_support_projects_only_at_end(self):
        key = ((0,) * 10, tuple(("G", 0) for _ in range(5)))
        seq = (("P", (0, 1), 0, 0),)
        support = structural_sequence_support(key, seq, 0, 2, charged_nodes=(0, 1))
        self.assertEqual(len(support), 4)
        self.assertTrue(all(k[1][0][0] == "C" and k[1][1][0] == "C" for k in support))
        self.assertTrue(all(k[1][2:] == tuple(("G", 0) for _ in range(3)) for k in support))

    def test_structural_support_guard_raises_without_truncating(self):
        key = ((0,) * 10, tuple(("G", 0) for _ in range(5)))
        seq = (("P", (0, 1), 0, 0),)
        with self.assertRaises(RuntimeError):
            structural_sequence_support(key, seq, 0, 2, charged_nodes=(0, 1), max_states=1)

    def test_CV_structural_support_contains_direct_scalar_volume_branch(self):
        key = ((0,) * 10, (("C", 0, 0, 0, 0),) + tuple(("G", 0) for _ in range(4)))
        support = C_volume_structural_support(key, 0, 1, 2, 0, 0, max_states=100)
        self.assertIn(key, support)
        self.assertTrue(all(k[1][0][0] == "C" for k in support))

    def test_fundamental_hit_changes_total_J_by_half_spin(self):
        key = ((2,) * 10, (("C", 2, 0, 0, 0),) + tuple(("G", 0) for _ in range(4)))
        states = list(advance_hit_label_support(key, 0, 1, 3))
        self.assertTrue(states)
        for _, labels in states:
            self.assertIn(labels[0][1], (1, 3))
            self.assertEqual(labels[1][1], 1)

    def test_source_node_local_spin_order(self):
        spins = tuple(range(10))
        self.assertEqual(local_spins(spins, 0), (0, 1, 2, 3))


if __name__ == "__main__":
    unittest.main()
