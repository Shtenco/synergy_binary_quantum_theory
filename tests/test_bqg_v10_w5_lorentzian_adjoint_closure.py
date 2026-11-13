import sys
import unittest
from unittest import mock
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import bqg_v10_w5_lorentzian_adjoint_closure as gate


def random_blocks(rng, dim, antihermitian):
    raw = rng.normal(size=(2 * dim, 2 * dim)) + 1j * rng.normal(
        size=(2 * dim, 2 * dim)
    )
    full = raw - raw.conj().T if antihermitian else raw + raw.conj().T
    return full.reshape(2, dim, 2, dim).transpose(0, 2, 1, 3)


def callback(blocks, antihermitian):
    def apply(state, source_v, edge, i, j, cutoff):
        vec = np.zeros(blocks[edge].shape[2], dtype=complex)
        for key, amp in state.items():
            vec[key[0]] += amp
        out_vec = blocks[edge][i, j] @ vec
        out = {
            (n, 0, 0): z
            for n, z in enumerate(out_vec)
            if abs(z) > 1e-14
        }
        if antihermitian:
            return out, {
                "outer_complete_basis_leakage": 0.0,
                "internal_volume_sector_leakage": 0.0,
                "complete_charge_basis_leakage": 0.0,
            }
        return out, 0.0
    return apply


def dense(state, dim):
    out = np.zeros(dim, dtype=complex)
    for key, amp in state.items():
        out[key[0]] += amp
    return out


class TestW5AdjointClosureAlgebra(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(20261110)
        self.dim = 4
        self.K = {e: random_blocks(rng, self.dim, True) for e in (1, 2, 3, 4)}
        self.V = {e: random_blocks(rng, self.dim, False) for e in (1, 2, 3, 4)}
        self.ck = callback(self.K, True)
        self.cv = callback(self.V, False)
        self.psi = {
            (0, 0, 0): 0.2 + 0.1j,
            (1, 0, 0): -0.3 + 0.2j,
            (2, 0, 0): 0.5 - 0.4j,
            (3, 0, 0): 0.1 + 0.3j,
        }

    def test_run_closure_memoizes_exact_repeated_CV_calls_and_replays_events(self):
        import bqg_pruning_audit as audit
        import peter_weyl_covariant_volume_leg_gate as CV
        import peter_weyl_lorentzian_ordered_triple_gate as raw
        import peter_weyl_lorentzian_sine_ordered_triple_gate as lsine

        original_component = raw.COMP.C_volume_component
        original_forward = gate.triple_forward
        original_adjoint = gate.triple_adjoint
        original_gauss_to_covariant = CV.gauss_to_covariant
        original_lsine_run = lsine.run
        component_calls = []
        component_results = []

        def fake_component(state, source_v, target_v, i, j, cutoff):
            component_calls.append((
                tuple(sorted(state.items(), key=lambda kv: repr(kv[0]))),
                source_v, target_v, i, j, cutoff,
            ))
            audit.record_threshold(
                "test.memoized_CV", (source_v, target_v, i, j),
                1e-15 + 0j, 1e-14,
            )
            return ({("memo", source_v, target_v, i, j): 0.5 + 0.25j}, 0.0)

        def fake_gauss_to_covariant(state, source_v):
            seed_key = next(iter(state))
            return {("base", seed_key[1][0]): 1.0 + 0j}

        def fake_forward(psi_cov, source_v, a, b, c, Jmax2, cv_component, ck_component):
            for i in range(2):
                for j in range(2):
                    first = cv_component(psi_cov, source_v, c, i, j, Jmax2)
                    second = cv_component(psi_cov, source_v, c, i, j, Jmax2)
                    self.assertEqual(first, second)
                    component_results.append(first)
            return {}, 0.0, {}

        def fake_adjoint(psi_cov, source_v, a, b, c, Jmax2, cv_component, ck_component):
            intermediate = {("adjoint-intermediate", a, b, c): 1.0 + 0j}
            cv_component(intermediate, source_v, c, 0, 0, Jmax2)
            return {}, 0.0, {}

        def fake_lsine_run(source_v):
            raw.run()
            return {"passed": False}

        try:
            raw.COMP.C_volume_component = fake_component
            with mock.patch.object(CV, "gauss_to_covariant", fake_gauss_to_covariant), \
                 mock.patch.object(gate, "triple_forward", fake_forward), \
                 mock.patch.object(gate, "triple_adjoint", fake_adjoint), \
                 mock.patch.object(lsine, "run", fake_lsine_run):
                cached = gate.run_closure(
                    source_nodes=(0,), Jmax2=10, max_states=500,
                    max_depth=0, amp_tol=1e-12, cache_mode="on",
                )
        finally:
            raw.COMP.C_volume_component = original_component
            gate.triple_forward = original_forward
            gate.triple_adjoint = original_adjoint
            CV.gauss_to_covariant = original_gauss_to_covariant
            lsine.run = original_lsine_run

        cache = cached["lower_level_pruning_ledger"]["component_cache"]["C_volume_component"]
        self.assertEqual(cache["maxsize"], 16)
        self.assertEqual(cache["misses"], 80)
        self.assertEqual(cache["hits"], 880)
        self.assertEqual(len(component_calls), 200)
        cached_results = list(component_results)
        cached_events = [
            e for e in cached["lower_level_pruning_ledger"]["events"]
            if e["site_id"] == "test.memoized_CV"
        ]
        self.assertEqual(len(cached_results), 480)
        self.assertEqual(len(cached_events), 1080)

        component_calls.clear()
        component_results.clear()
        try:
            raw.COMP.C_volume_component = fake_component
            with mock.patch.object(CV, "gauss_to_covariant", fake_gauss_to_covariant), \
                 mock.patch.object(gate, "triple_forward", fake_forward), \
                 mock.patch.object(gate, "triple_adjoint", fake_adjoint), \
                 mock.patch.object(lsine, "run", fake_lsine_run):
                uncached = gate.run_closure(
                    source_nodes=(0,), Jmax2=10, max_states=500,
                    max_depth=0, amp_tol=1e-12, cache_mode="off",
                )
        finally:
            raw.COMP.C_volume_component = original_component
            gate.triple_forward = original_forward
            gate.triple_adjoint = original_adjoint
            CV.gauss_to_covariant = original_gauss_to_covariant
            lsine.run = original_lsine_run

        off_cache = uncached["lower_level_pruning_ledger"]["component_cache"]["C_volume_component"]
        self.assertEqual(off_cache["hits"], 0)
        self.assertEqual(off_cache["misses"], 0)
        self.assertEqual(len(component_calls), 1080)
        self.assertEqual(component_results, cached_results)
        off_events = [
            e for e in uncached["lower_level_pruning_ledger"]["events"]
            if e["site_id"] == "test.memoized_CV"
        ]
        self.assertEqual(off_events, cached_events)

    def test_single_triple_direct_adjoint_matches_conjugate_transpose(self):
        forward, _, _ = gate.triple_forward(
            self.psi, 0, 1, 2, 3, 10, self.cv, self.ck
        )
        direct, _, _ = gate.triple_adjoint(
            self.psi, 0, 1, 2, 3, 10, self.cv, self.ck
        )
        x = dense(self.psi, self.dim)
        y = dense(forward, self.dim)
        z = dense(direct, self.dim)
        # Compare the two maps as matrices, not only on one vector.
        F = np.column_stack([
            dense(gate.triple_forward(
                {(n, 0, 0): 1.0}, 0, 1, 2, 3, 10, self.cv, self.ck
            )[0], self.dim)
            for n in range(self.dim)
        ])
        D = np.column_stack([
            dense(gate.triple_adjoint(
                {(n, 0, 0): 1.0}, 0, 1, 2, 3, 10, self.cv, self.ck
            )[0], self.dim)
            for n in range(self.dim)
        ])
        np.testing.assert_allclose(D, F.conj().T, rtol=2e-13, atol=2e-13)
        np.testing.assert_allclose(z, F.conj().T @ x, rtol=2e-13, atol=2e-13)
        np.testing.assert_allclose(y, F @ x, rtol=2e-13, atol=2e-13)

    def test_complete_epsilon_node_adjoint_matches_matrix_adjoint(self):
        triples = gate.epsilon_triples((1, 2, 3, 4))
        F = np.zeros((self.dim, self.dim), dtype=complex)
        D = np.zeros_like(F)
        for sign, (a, b, c) in triples:
            for n in range(self.dim):
                basis = {(n, 0, 0): 1.0}
                f, _, _ = gate.triple_forward(
                    basis, 0, a, b, c, 10, self.cv, self.ck
                )
                d, _, _ = gate.triple_adjoint(
                    basis, 0, a, b, c, 10, self.cv, self.ck
                )
                F[:, n] += sign * dense(f, self.dim)
                D[:, n] += sign * dense(d, self.dim)
        np.testing.assert_allclose(D, F.conj().T, rtol=5e-13, atol=5e-13)

    def test_closure_rejects_unknown_cache_mode_before_operator_execution(self):
        with self.assertRaisesRegex(ValueError, "cache_mode"):
            gate.run_closure(cache_mode="sometimes")

    def test_zeroaware_spectral_factor_reports_eigenvalue_sensitivity(self):
        import peter_weyl_zeroaware_volume_migration_experiment as ZVM

        original_factor = ZVM.ZERO_FACTOR
        original_audit = ZVM.SPECTRAL_AUDIT
        try:
            matrix = np.diag([0.0, 1.0e-12, 1.0]).astype(complex)
            ZVM.ZERO_FACTOR = 1000.0
            ZVM.SPECTRAL_AUDIT = {}
            _, tau_low, _ = ZVM.zeroaware_sqrt_abs(matrix)
            low = ZVM.spectral_audit_snapshot()
            ZVM.ZERO_FACTOR = 2000.0
            ZVM.SPECTRAL_AUDIT = {}
            _, tau_high, _ = ZVM.zeroaware_sqrt_abs(matrix)
            high = ZVM.spectral_audit_snapshot()

            self.assertGreater(tau_high, tau_low)
            self.assertEqual(low["zeroed_eigenvalue_count"], 1)
            self.assertEqual(low["retained_eigenvalue_count"], 2)
            self.assertEqual(high["zeroed_eigenvalue_count"], 2)
            self.assertEqual(high["retained_eigenvalue_count"], 1)
            self.assertEqual(high["by_dimension"]["3"]["matrix_count"], 1)
        finally:
            ZVM.ZERO_FACTOR = original_factor
            ZVM.SPECTRAL_AUDIT = original_audit

    def test_state_and_depth_caps_are_reported_independently(self):
        self.assertEqual(
            gate._closure_cap_reasons(500, 8, 500, 8),
            {"max_states": 500, "max_depth": 8},
        )
        self.assertEqual(
            gate._closure_cap_reasons(500, 7, 500, 8),
            {"max_states": 500},
        )
        self.assertEqual(
            gate._closure_cap_reasons(499, 8, 500, 8),
            {"max_depth": 8},
        )
        self.assertEqual(gate._closure_cap_reasons(1000, 100, None, None), {})

    def test_final_closure_gate_requires_lsine_gate(self):
        common = dict(
            forward_norm=1.0,
            adjoint_norm=1.0,
            adjoint_relative_residual=0.0,
            scalar_projection_rejected_terms=0,
            max_component_leakage=0.0,
            max_spin_reached=3.5,
            spin_wall_omitted_branches=0,
            max_intermediate_spin_reached=3.5,
        )
        self.assertTrue(gate.closure_acceptance(
            **common, sine_engine_result={"passed": True}
        ))
        self.assertFalse(gate.closure_acceptance(
            **common, sine_engine_result={"passed": False}
        ))
        self.assertFalse(gate.closure_acceptance(
            **common, sine_engine_result=None
        ))

    def test_final_closure_gate_rejects_spin_wall_violation(self):
        common = dict(
            forward_norm=1.0,
            adjoint_norm=1.0,
            adjoint_relative_residual=0.0,
            scalar_projection_rejected_terms=0,
            max_component_leakage=0.0,
            sine_engine_result={"passed": True},
        )
        self.assertFalse(gate.closure_acceptance(
            **common, max_spin_reached=3.5000001
        ))


    def test_final_closure_gate_rejects_unobserved_spin_wall_branches(self):
        common = dict(
            forward_norm=1.0,
            adjoint_norm=1.0,
            adjoint_relative_residual=0.0,
            scalar_projection_rejected_terms=0,
            max_component_leakage=0.0,
            max_spin_reached=3.5,
            sine_engine_result={"passed": True},
        )
        self.assertFalse(gate.closure_acceptance(
            **common, spin_wall_omitted_branches=1
        ))

    def test_spin_wall_probe_records_complex_amplitude_before_rejection(self):
        stats = {
            "count": 0, "l1": 0.0, "l2_sq": 0.0, "max_abs": 0.0,
            "by_edge": {}, "samples": [],
            "max_intermediate_spin_reached": 0.0,
            "intermediate_above_3_5_count": 0,
            "intermediate_above_3_5_l1": 0.0,
            "intermediate_above_3_5_l2_sq": 0.0,
            "intermediate_samples": [],
        }

        def primitive(branch, x, y, i, j, wall):
            spins, tensors, amp = branch
            edge_spin = spins[0]
            candidates = []
            for next_spin in (edge_spin - 1, edge_spin + 1):
                if 0 <= next_spin <= wall:
                    next_spins = (next_spin,)
                    candidates.append((next_spins, tensors, amp * (1 + 1j)))
            return candidates

        probe = gate._make_spin_wall_probe(
            primitive, {(0, 1): 0}, 4, stats
        )
        actual = probe(((4,), (), 2 + 0j), 0, 1, 0, 1, 4)
        self.assertEqual([item[0] for item in actual], [(3,)])
        self.assertEqual(stats["count"], 1)
        self.assertAlmostEqual(stats["l1"], abs(2 + 2j))
        self.assertAlmostEqual(stats["l2_sq"], abs(2 + 2j) ** 2)
        self.assertEqual(stats["samples"][0]["omitted_output_spins"], [5])


    def test_final_closure_gate_rejects_intermediate_spin_excursion(self):
        common = dict(
            forward_norm=1.0,
            adjoint_norm=1.0,
            adjoint_relative_residual=0.0,
            scalar_projection_rejected_terms=0,
            max_component_leakage=0.0,
            max_spin_reached=3.5,
            sine_engine_result={"passed": True},
            spin_wall_omitted_branches=0,
        )
        self.assertFalse(gate.closure_acceptance(
            **common, max_intermediate_spin_reached=4.0
        ))

    def test_spin_probe_tracks_retained_intermediate_spin_above_3_5(self):
        stats = {
            "count": 0, "l1": 0.0, "l2_sq": 0.0, "max_abs": 0.0,
            "by_edge": {}, "samples": [],
            "max_intermediate_spin_reached": 0.0,
            "intermediate_above_3_5_count": 0,
            "intermediate_above_3_5_l1": 0.0,
            "intermediate_above_3_5_l2_sq": 0.0,
            "intermediate_samples": [],
        }

        def primitive(branch, x, y, i, j, wall):
            spins, tensors, amp = branch
            return [
                ((next_spin,), tensors, amp * (1 + 1j))
                for next_spin in (spins[0] - 1, spins[0] + 1)
                if 0 <= next_spin <= wall
            ]

        probe = gate._make_spin_wall_probe(
            primitive, {(0, 1): 0}, 10, stats
        )
        probe(((7,), (), 1 + 0j), 0, 1, 0, 1, 10)
        self.assertEqual(stats["max_intermediate_spin_reached"], 4.0)
        self.assertEqual(stats["intermediate_above_3_5_count"], 1)


    def test_factorized_branch_metrics_include_tensor_coefficients(self):
        import numpy as np

        metrics = gate._factorized_branch_metrics(
            (np.array([1.0 + 0j, 0.0 + 2j]),),
            2.0 + 0j,
            include_payload=True,
        )
        self.assertAlmostEqual(metrics["branch_l1"], 6.0)
        self.assertAlmostEqual(metrics["branch_l2"], 2.0 * np.sqrt(5.0))
        self.assertAlmostEqual(metrics["branch_max_abs"], 4.0)
        self.assertEqual(metrics["tensor_shapes"], [[2]])
        self.assertEqual(len(metrics["tensor_sha256"]), 64)
        self.assertEqual(metrics["tensor_payload"][0]["real"], [1.0, 0.0])
        self.assertEqual(metrics["tensor_payload"][0]["imag"], [0.0, 2.0])

    def test_epsilon_orbit_contains_24_unique_terms(self):
        triples = gate.epsilon_triples((1, 2, 3, 4))
        self.assertEqual(len(triples), 24)
        self.assertEqual(len({triple for _, triple in triples}), 24)

    def test_sparse_closed_matrix_adjoint_indexing(self):
        keys = [("basis", 0), ("basis", 1), ("basis", 2)]
        # Columns of a complex 3x3 forward matrix, represented sparsely.
        F = [
            {keys[0]: 1 + 2j, keys[2]: -0.5j},
            {keys[1]: 3 - 1j},
            {keys[0]: 2 + 0.25j, keys[1]: -4j},
        ]
        index = {key: i for i, key in enumerate(keys)}
        expected = gate._matrix_adjoint_expected_columns(F, keys, index)
        dense_F = np.column_stack([
            np.array([col.get(k, 0j) for k in keys], dtype=complex)
            for col in F
        ])
        dense_expected = np.column_stack([
            np.array([col.get(k, 0j) for k in keys], dtype=complex)
            for col in expected
        ])
        np.testing.assert_allclose(dense_expected, dense_F.conj().T)
        err, scale = gate.sparse_matrix_residual(expected, expected)
        self.assertEqual(err, 0.0)
        self.assertGreater(scale, 0.0)

    def test_antihermitian_commutator_to_hermitian_K_sign_and_hbar(self):
        rng = np.random.default_rng(20261110)
        x = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        y = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        V = x + x.conj().T
        HE = y + y.conj().T
        hbar = 1.37
        Kcomm = V @ HE - HE @ V
        self.assertTrue(np.allclose(Kcomm.conj().T, -Kcomm, atol=2e-13))
        KH = Kcomm / (1j * hbar)
        self.assertTrue(np.allclose(KH.conj().T, KH, atol=2e-13))
        np.testing.assert_allclose(
            KH @ KH,
            -(Kcomm @ Kcomm) / (hbar ** 2),
            rtol=2e-13,
            atol=2e-13,
        )

    def test_symmetric_completion_is_hermitian(self):
        rng = np.random.default_rng(20261111)
        B = rng.normal(size=(5, 5)) + 1j * rng.normal(size=(5, 5))
        Bsym = 0.5 * (B + B.conj().T)
        np.testing.assert_allclose(Bsym, Bsym.conj().T, atol=2e-13)

    def test_sparse_ordering_parts_have_correct_adjoint_signs(self):
        rng = np.random.default_rng(20261112)
        B = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        keys = [("state", i) for i in range(4)]
        forward = [
            {keys[row]: B[row, col] for row in range(4) if B[row, col] != 0}
            for col in range(4)
        ]
        Bdag = B.conj().T
        adjoint = [
            {keys[row]: Bdag[row, col] for row in range(4) if Bdag[row, col] != 0}
            for col in range(4)
        ]
        sym, anti, kh = gate.sparse_ordering_parts(forward, adjoint)

        def dense_columns(columns):
            return np.column_stack([
                np.array([col.get(key, 0j) for key in keys], dtype=complex)
                for col in columns
            ])

        S = dense_columns(sym)
        A = dense_columns(anti)
        KH = dense_columns(kh)
        np.testing.assert_allclose(S, S.conj().T, rtol=2e-13, atol=2e-13)
        np.testing.assert_allclose(A, -A.conj().T, rtol=2e-13, atol=2e-13)
        np.testing.assert_allclose(KH, -S, rtol=2e-13, atol=2e-13)
        np.testing.assert_allclose(S + A, B, rtol=2e-13, atol=2e-13)

    def test_covariant_scalar_mapping_rejects_open_channels(self):
        spins = (4,) * 10
        valid = (spins, (-1, 0, 0, 0, 0), 0, 0, 0, 0)
        mapped = gate.covariant_scalar_to_gauss(valid, 0)
        self.assertEqual(mapped, (spins, (0, 0, 0, 0, 0)))
        invalid = (spins, (-1, 0, 0, 0, 0), 2, 0, 0, 0)
        with self.assertRaises(ValueError):
            gate.covariant_scalar_to_gauss(invalid, 0)


if __name__ == "__main__":
    unittest.main()
