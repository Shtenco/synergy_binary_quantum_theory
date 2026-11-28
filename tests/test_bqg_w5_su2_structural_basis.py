import unittest

import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_covariant_K_leg_gate as CK

from bqg_w5_su2_structural_basis import (
    allowed_pair2, all_total_J2, charged_labels, charged_node_basis, covariant_basis_size,
    gauss_labels, iter_covariant_keys, local_spins, validate_covariant_key,
    w5_gauss_seed, w5_seed_as_covariant, spin_support_for_sequence, t_sequence_spin_support, spin_sector_basis_profile,
    charged_nodes_basis_size, charged_nodes_basis_profile, volume_label_targets, lift_key_to_all_charged,
    project_all_charged_key, advance_hit_label_support, structural_sequence_support,
    C_volume_structural_support, factorized_t_sequence_support,
    factorized_sequence_support, factorized_contains_key, factorized_support_profile,
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
        self.assertEqual({s[1][0][2] for s in states}, {-1, 1})
        self.assertEqual({s[1][1][2] for s in states}, {-1, 1})
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
            structural_sequence_support(key, seq, 0, 2, charged_nodes=(0, 1), max_states=0)

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

    def test_factorized_support_contains_explicit_one_hit_support(self):
        key = ((0,) * 10, tuple(("G", 0) for _ in range(5)))
        seq = (("P", (0, 1), 0, 0),)
        explicit = structural_sequence_support(key, seq, 0, 2, charged_nodes=(0, 1))
        factored = factorized_sequence_support(key, seq, 0, 2, charged_nodes=(0, 1))
        self.assertTrue(all(factorized_contains_key(factored, item) for item in explicit))
        profile = factorized_support_profile(factored)
        self.assertEqual(profile["spin_sector_count"], 1)
        self.assertGreaterEqual(profile["product_basis_dimension_upper_bound"], len(explicit))

    def test_factorized_support_contains_explicit_path_support(self):
        key = ((0,) * 10, tuple(("G", 0) for _ in range(5)))
        seq = (("P", (0, 1, 2), 0, 1),)
        explicit = structural_sequence_support(key, seq, 0, 2, charged_nodes=(0, 2))
        factored = factorized_sequence_support(key, seq, 0, 2, charged_nodes=(0, 2))
        self.assertTrue(all(factorized_contains_key(factored, item) for item in explicit))

    def test_source_node_local_spin_order(self):
        spins = tuple(range(10))
        self.assertEqual(local_spins(spins, 0), (0, 1, 2, 3))


    def test_K_complete_structural_support_tiny_sector(self):
        # Replace the real H_E word family with one empty primitive word.
        # This isolates the commutator's two support branches without claiming
        # anything about the full W5 KKV expansion.
        import bqg_w5_su2_structural_basis as SB
        from unittest.mock import patch

        key = ((0,) * 10,
               (("C", 0, 0, 0, 0),) +
               tuple(("G", 0) for _ in range(4)))
        with patch.object(PW, "oriented_specs",
                          return_value=((1, (0, 1, 2, 3)),)), \
             patch.object(PW, "T_sequences",
                          return_value=((1.0, ()),)):
            support = SB.k_complete_structural_support(
                key, 0, 1, charged_nodes=(0,), max_states=100
            )
        self.assertIn(key, support)
        self.assertTrue(all(item[1][0][0] == "C" for item in support))
        self.assertTrue(all(
            all(lab[0] == "G" for lab in item[1][1:])
            for item in support
        ))

    def test_CK_structural_support_tiny_sector_and_offdiagonal(self):
        # A single empty H_E word lets us test outer inverse/closure ordering
        # on the smallest nontrivial fundamental edge-spin transition.
        import bqg_w5_su2_structural_basis as SB
        from unittest.mock import patch

        key = ((0,) * 10,
               (("C", 0, 0, 0, 0),) +
               tuple(("G", 0) for _ in range(4)))
        with patch.object(PW, "oriented_specs",
                          return_value=((1, (0, 1, 2, 3)),)), \
             patch.object(PW, "T_sequences",
                          return_value=((1.0, ()),)):
            diagonal = SB.C_K_structural_support(
                key, 0, 1, 1, 0, 0, max_states=100
            )
            offdiagonal = SB.C_K_structural_support(
                key, 0, 1, 1, 0, 1, max_states=100
            )
        self.assertTrue(diagonal)
        self.assertTrue(offdiagonal)
        self.assertTrue(all(item[1][0][0] == "C" for item in diagonal))
        self.assertTrue(all(item[1][1][0] == "G" for item in diagonal))
        self.assertTrue(all(item[1][0][0] == "C" for item in offdiagonal))
        self.assertTrue(all(item[1][1][0] == "G" for item in offdiagonal))

    def test_CK_structural_support_guard_is_fail_closed(self):
        import bqg_w5_su2_structural_basis as SB
        from unittest.mock import patch

        key = ((0,) * 10,
               (("C", 0, 0, 0, 0),) +
               tuple(("G", 0) for _ in range(4)))
        with patch.object(PW, "oriented_specs",
                          return_value=((1, (0, 1, 2, 3)),)), \
             patch.object(PW, "T_sequences",
                          return_value=((1.0, ()),)):
            with self.assertRaises(RuntimeError):
                SB.C_K_structural_support(
                    key, 0, 1, 1, 0, 0, max_states=0
                )


    def test_factorized_operator_K_contains_exact_tiny_support(self):
        import bqg_w5_su2_structural_basis as SB
        import bqg_w5_factorized_operator_support as FO
        from unittest.mock import patch

        key = ((0,) * 10,
               (("C", 0, 0, 0, 0),) +
               tuple(("G", 0) for _ in range(4)))
        lifted = SB.lift_key_to_all_charged(key)
        initial = {lifted[0]: tuple(frozenset((x,)) for x in lifted[1])}
        with (
            patch.object(PW, "oriented_specs",
                         return_value=((1, (0, 1, 2, 3)),)),
            patch.object(PW, "T_sequences",
                         return_value=((1.0, ()),)),
        ):
            exact = SB.k_complete_structural_support(
                key, 0, 1, charged_nodes=(0,), max_states=100
            )
            factored = FO.factorized_k(
                initial, 0, 1, charged_nodes=(0,), max_factor_labels=1000
            )
        self.assertTrue(all(FO.SB.factorized_contains_key(factored, x)
                            for x in exact))

    def test_factorized_operator_CK_contains_exact_tiny_support(self):
        import bqg_w5_su2_structural_basis as SB
        import bqg_w5_factorized_operator_support as FO
        from unittest.mock import patch

        key = ((0,) * 10,
               (("C", 0, 0, 0, 0),) +
               tuple(("G", 0) for _ in range(4)))
        lifted = SB.lift_key_to_all_charged(key)
        initial = {lifted[0]: tuple(frozenset((x,)) for x in lifted[1])}
        with (
            patch.object(PW, "oriented_specs",
                         return_value=((1, (0, 1, 2, 3)),)),
            patch.object(PW, "T_sequences",
                         return_value=((1.0, ()),)),
        ):
            exact = SB.C_K_structural_support(
                key, 0, 1, 1, 0, 0, max_states=100
            )
            factored = FO.factorized_c_k(
                initial, 0, 1, 1, 0, 0, max_factor_labels=1000
            )
        self.assertTrue(all(FO.SB.factorized_contains_key(factored, x)
                            for x in exact))

    def test_factorized_operator_guard_fails_closed(self):
        import bqg_w5_factorized_operator_support as FO
        support = {tuple([0] * 10): tuple(
            frozenset((("C", 0, 0, 0, 0),)) for _ in range(5)
        )}
        with self.assertRaises(RuntimeError):
            FO._volume(support, 0, 0)


    def test_factorized_forward_KKV_contains_exact_tiny_support(self):
        import bqg_w5_su2_structural_basis as SB
        import bqg_w5_factorized_operator_support as FO
        from unittest.mock import patch

        key = ((0,) * 10,
               (("C", 0, 0, 0, 0),) +
               tuple(("G", 0) for _ in range(4)))
        lifted = SB.lift_key_to_all_charged(key)
        initial = {lifted[0]: tuple(frozenset((x,)) for x in lifted[1])}
        with (
            patch.object(PW, "oriented_specs",
                         return_value=((1, (0, 1, 2, 3)),)),
            patch.object(PW, "T_sequences",
                         return_value=((1.0, ()),)),
        ):
            exact = SB.triple_forward_structural_support(
                key, 0, 1, 2, 3, 1, max_states=1000
            )
            factored = FO.factorized_triple(
                initial, 0, 1, 2, 3, 1, adjoint=False,
                max_factor_labels=5000
            )
        self.assertTrue(all(SB.factorized_contains_key(factored, x)
                            for x in exact))

    def test_factorized_reversed_adjoint_KKV_contains_exact_tiny_support(self):
        import bqg_w5_su2_structural_basis as SB
        import bqg_w5_factorized_operator_support as FO
        from unittest.mock import patch

        key = ((0,) * 10,
               (("C", 0, 0, 0, 0),) +
               tuple(("G", 0) for _ in range(4)))
        lifted = SB.lift_key_to_all_charged(key)
        initial = {lifted[0]: tuple(frozenset((x,)) for x in lifted[1])}
        with (
            patch.object(PW, "oriented_specs",
                         return_value=((1, (0, 1, 2, 3)),)),
            patch.object(PW, "T_sequences",
                         return_value=((1.0, ()),)),
        ):
            exact = SB.triple_adjoint_structural_support(
                key, 0, 1, 2, 3, 1, max_states=1000
            )
            factored = FO.factorized_triple(
                initial, 0, 1, 2, 3, 1, adjoint=True,
                max_factor_labels=5000
            )
        self.assertTrue(all(SB.factorized_contains_key(factored, x)
                            for x in exact))


    def test_factorized_node_KKV_contains_exact_24_triple_union(self):
        import bqg_w5_su2_structural_basis as SB
        import bqg_w5_factorized_operator_support as FO
        from unittest.mock import patch

        key = SB.w5_seed_as_covariant(0, 0)
        lifted = SB.lift_key_to_all_charged(key)
        initial = {lifted[0]: tuple(frozenset((x,)) for x in lifted[1])}
        with (
            patch.object(PW, "oriented_specs",
                         return_value=((1, (0, 1, 2, 3)),)),
            patch.object(PW, "T_sequences",
                         return_value=((1.0, ()),)),
        ):
            exact = SB.node_kkv_structural_support(
                key, 0, 2, max_states=5000
            )
            factored = FO.factorized_node_kkv(
                initial, 0, 2, max_factor_labels=10000
            )
        self.assertTrue(all(SB.factorized_contains_key(factored, x)
                            for x in exact))


    def test_pairwise_factor_graph_contains_exact_one_hit_support(self):
        import bqg_w5_pairwise_factor_graph as FG
        import bqg_w5_su2_structural_basis as SB
        key = ((0,) * 10, tuple(("G", 0) for _ in range(5)))
        exact = SB.structural_sequence_support(
            key, (("P", (0, 1)),), 0, 2, charged_nodes=(0, 1)
        )
        pairwise = FG._sequence(
            FG._initial(key), (("P", (0, 1)),), 0, 2,
            (0, 1), 10000
        )
        self.assertTrue(all(FG.contains_joint(pairwise, x) for x in exact))
        self.assertGreater(
            sum(len(x) for s in pairwise.values() for x in s["pairs"].values()), 0
        )

    def test_pairwise_factor_graph_contains_exact_volume_support(self):
        import bqg_w5_pairwise_factor_graph as FG
        import bqg_w5_su2_structural_basis as SB
        key = ((4,) * 10,
               (("C", 0, 0, 0, 0),) + tuple(("G", 0) for _ in range(4)))
        lifted = SB.lift_key_to_all_charged(key)
        exact = set(SB.advance_volume_label_support(lifted, 0))
        pairwise = FG._volume(FG._initial(key), 0, 10000)
        self.assertTrue(all(FG.contains_joint(pairwise, x) for x in exact))

    def test_pairwise_factor_graph_K_and_CK_include_exact_mock_support(self):
        import bqg_w5_pairwise_factor_graph as FG
        import bqg_w5_su2_structural_basis as SB
        from unittest.mock import patch
        key = ((0,) * 10,
               (("C", 0, 0, 0, 0),) + tuple(("G", 0) for _ in range(4)))
        initial = FG._initial(key)
        with (
            patch.object(PW, "oriented_specs",
                         return_value=((1, (0, 1, 2, 3)),)),
            patch.object(PW, "T_sequences",
                         return_value=((1.0, ()),)),
        ):
            exact_k = SB.k_complete_structural_support(
                key, 0, 1, charged_nodes=(0,), max_states=100
            )
            pair_k = FG.factor_graph_k(initial, 0, 1, (0,), 10000)
            self.assertTrue(all(FG.contains_joint(pair_k, x) for x in exact_k))
            exact_ck = SB.C_K_structural_support(
                key, 0, 1, 1, 0, 0, max_states=100
            )
            pair_ck = FG.factor_graph_c_k(initial, 0, 1, 1, 0, 0, 10000)
            self.assertTrue(all(FG.contains_joint(pair_ck, x) for x in exact_ck))

    def test_pairwise_factor_graph_forward_and_adjoint_KKV_include_exact_mock(self):
        import bqg_w5_pairwise_factor_graph as FG
        import bqg_w5_su2_structural_basis as SB
        from unittest.mock import patch
        key = ((0,) * 10,
               (("C", 0, 0, 0, 0),) + tuple(("G", 0) for _ in range(4)))
        initial = FG._initial(key)
        with (
            patch.object(PW, "oriented_specs",
                         return_value=((1, (0, 1, 2, 3)),)),
            patch.object(PW, "T_sequences",
                         return_value=((1.0, ()),)),
        ):
            exact_f = SB.triple_forward_structural_support(
                key, 0, 1, 2, 3, 1, max_states=1000
            )
            pair_f = FG.factor_graph_triple(
                initial, 0, 1, 2, 3, 1, adjoint=False, limit=50000
            )
            self.assertTrue(all(FG.contains_joint(pair_f, x) for x in exact_f))
            exact_a = SB.triple_adjoint_structural_support(
                key, 0, 1, 2, 3, 1, max_states=1000
            )
            pair_a = FG.factor_graph_triple(
                initial, 0, 1, 2, 3, 1, adjoint=True, limit=50000
            )
            self.assertTrue(all(FG.contains_joint(pair_a, x) for x in exact_a))

    def test_pairwise_factor_graph_guard_is_fail_closed(self):
        import bqg_w5_pairwise_factor_graph as FG
        key = ((0,) * 10,
               (("C", 0, 0, 0, 0),) + tuple(("G", 0) for _ in range(4)))
        with self.assertRaises(RuntimeError):
            FG._volume(FG._initial(key), 0, 0)


    def test_triplet_cluster_contains_exact_one_hit_support(self):
        import bqg_w5_triplet_cluster_support as TG
        import bqg_w5_su2_structural_basis as SB
        key = ((0,) * 10, tuple(("G", 0) for _ in range(5)))
        exact = SB.structural_sequence_support(
            key, (("P", (0, 1)),), 0, 2, charged_nodes=(0, 1)
        )
        got = TG._sequence(TG._initial(key), (("P", (0, 1)),),
                           0, 2, (0, 1), 20000)
        self.assertTrue(all(TG.contains_joint(got, x) for x in exact))

    def test_triplet_cluster_forward_and_adjoint_KKV_include_exact_mock(self):
        import bqg_w5_triplet_cluster_support as TG
        import bqg_w5_su2_structural_basis as SB
        from unittest.mock import patch
        key = ((0,) * 10,
               (("C", 0, 0, 0, 0),) + tuple(("G", 0) for _ in range(4)))
        initial = TG._initial(key)
        with (
            patch.object(PW, "oriented_specs",
                         return_value=((1, (0, 1, 2, 3)),)),
            patch.object(PW, "T_sequences",
                         return_value=((1.0, ()),)),
        ):
            exact_f = SB.triple_forward_structural_support(
                key, 0, 1, 2, 3, 1, max_states=1000
            )
            got_f = TG.factor_graph_triple(
                initial, 0, 1, 2, 3, 1, adjoint=False, limit=100000
            )
            self.assertTrue(all(TG.contains_joint(got_f, x) for x in exact_f))
            exact_a = SB.triple_adjoint_structural_support(
                key, 0, 1, 2, 3, 1, max_states=1000
            )
            got_a = TG.factor_graph_triple(
                initial, 0, 1, 2, 3, 1, adjoint=True, limit=100000
            )
            self.assertTrue(all(TG.contains_joint(got_a, x) for x in exact_a))

    def test_triplet_cluster_guard_is_fail_closed(self):
        import bqg_w5_triplet_cluster_support as TG
        key = ((0,) * 10,
               (("C", 0, 0, 0, 0),) + tuple(("G", 0) for _ in range(4)))
        with self.assertRaises(RuntimeError):
            TG._volume(TG._initial(key), 0, 0)


if __name__ == "__main__":
    unittest.main()
