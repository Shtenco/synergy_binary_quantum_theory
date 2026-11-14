import sys
import unittest
import itertools
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import bqg_pruning_audit as audit
import charged_intertwiner_recoupling_gate as CHARGED
import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_covariant_K_leg_gate as CK
import peter_weyl_covariant_volume_leg_gate as CV
import peter_weyl_covariant_composition_gate as COMP
import peter_weyl_lorentzian_ordered_triple_gate as RAW
import peter_weyl_covariant_K_projection_audit_gate as KPA


class TestPruningAuditOperatorWiring(unittest.TestCase):
    def test_prune_state_and_complete_to_covariant_projection_are_observed(self):
        spins = (4,) * 10
        with audit.capture() as collector:
            pruned = PW.prune_state({"small": 0.5e-8}, 1e-8)
            source_label = ("C", 0, 0, 0, 0)
            labels = (source_label,) + tuple(("G", 0) for _ in range(4))
            projected = COMP.complete_to_covariant(
                {(spins, labels): 0.5e-11}, 0, tol=1e-11
            )

        self.assertEqual(pruned, {})
        self.assertEqual(projected, {})
        sites = collector.snapshot()["by_site"]
        self.assertIn("k5.prune_state", sites)
        self.assertIn("covariant_composition.complete_to_covariant", sites)

    def test_shadow_mode_preserves_subthreshold_coefficients_at_wired_sites(self):
        with audit.capture() as collector:
            with audit.bypass():
                a = {}
                PW.add_dict(a, {"k": 0.5e-13})
                b = {}
                CK.add_state(b, {"k": 0.5e-11})
                c = {}
                COMP.add(c, {"k": 0.5e-11})
                d = {}
                RAW.add(d, {"k": 0.5e-11})
                e = PW.prune_state({"k": 0.5e-8}, 1e-8)
                source_label = ("C", 0, 0, 0, 0)
                labels = (source_label,) + tuple(("G", 0) for _ in range(4))
                f = COMP.complete_to_covariant(
                    {((4,) * 10, labels): 0.5e-11}, 0, tol=1e-11
                )

        self.assertEqual(a, {"k": 0.5e-13})
        self.assertEqual(b, {"k": 0.5e-11})
        self.assertEqual(c, {"k": 0.5e-11})
        self.assertEqual(d, {"k": 0.5e-11})
        self.assertEqual(e, {"k": 0.5e-8})
        self.assertEqual(f, {((4,) * 10, (-1, 0, 0, 0, 0), 0, 0, 0, 0): 0.5e-11})
        events = collector.snapshot()["events"]
        self.assertEqual(len(events), 6)
        self.assertTrue(all(x["kind"] == "would_prune_bypassed" for x in events))

    def test_independent_K_projection_ladder_prune_is_instrumented(self):
        source = {"small": 0.5e-12, "large": 2e-12}
        with audit.capture() as collector:
            production = KPA.prune(source, 1e-12)
            with audit.bypass():
                shadow = KPA.prune(source, 1e-12)
        self.assertEqual(production, {"large": 2e-12})
        self.assertEqual(shadow, source)
        events = collector.snapshot()["events"]
        self.assertEqual(len(events), 2)
        self.assertTrue(all(x["site_id"] == "covariant_K_projection_audit.ladder_prune" for x in events))
        self.assertEqual(
            [x["kind"] for x in events],
            ["coefficient_threshold", "would_prune_bypassed"],
        )

    def test_core_operator_merge_sites_are_observed_without_semantic_change(self):
        with audit.capture() as collector:
            a = {}
            PW.add_dict(a, {"k": 0.5e-13})
            b = {}
            CK.add_state(b, {"k": 0.5e-11})
            c = {}
            COMP.add(c, {"k": 0.5e-11})
            d = {}
            RAW.add(d, {"k": 0.5e-11})

        self.assertEqual(a, {})
        self.assertEqual(b, {})
        self.assertEqual(c, {})
        self.assertEqual(d, {})
        by_site = collector.snapshot()["by_site"]
        self.assertEqual(
            set(by_site),
            {
                "k5.add_dict.1e-13",
                "covariant_K_leg.add_state",
                "covariant_composition.add",
                "lorentzian_ordered_triple.add",
            },
        )
        self.assertTrue(all(row["coefficient_count"] == 1 for row in by_site.values()))



class TestW5TransitivePruningCoverage(unittest.TestCase):
    """Static guard over every repository-local import reachable from V10.

    Threshold-looking comparisons must either use the shared PRUNE policy (or
    delegate to PW.prune_state), or be explicitly classified as an exact
    selection rule, exact-nonzero bookkeeping, or spectral operator convention.
    A newly imported local module automatically enters this audit.
    """

    def test_v10_local_import_closure_and_threshold_classification(self):
        import ast

        scripts = Path(__file__).resolve().parents[1] / "scripts"
        pending = ["bqg_v10_w5_lorentzian_adjoint_closure"]
        modules = set()
        while pending:
            module_name = pending.pop()
            if module_name in modules:
                continue
            source_path = scripts / (module_name + ".py")
            if not source_path.is_file():
                continue
            modules.add(module_name)
            tree = ast.parse(source_path.read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    pending.extend(
                        alias.name for alias in node.names
                        if (scripts / (alias.name.split(".")[0] + ".py")).is_file()
                    )
                elif isinstance(node, ast.ImportFrom) and node.module:
                    local_name = node.module.split(".")[0]
                    if (scripts / (local_name + ".py")).is_file():
                        pending.append(local_name)

        expected_modules = {
            "bqg_pruning_audit",
            "bqg_pruning_coverage_inventory",
            "bqg_sparse_linear_accumulator",
            "bqg_v10_w5_lorentzian_adjoint_closure",
            "charged_intertwiner_recoupling_gate",
            "k5_peter_weyl_safe_hda_column",
            "peter_weyl_covariant_K_composition_cached_gate",
            "peter_weyl_covariant_K_composition_gate",
            "peter_weyl_covariant_K_leg_gate",
            "peter_weyl_covariant_K_projection_audit_gate",
            "peter_weyl_covariant_K_sine_composition_gate",
            "peter_weyl_covariant_composition_gate",
            "peter_weyl_covariant_volume_leg_gate",
            "peter_weyl_euclidean_sine_ordering_gate",
            "peter_weyl_lorentzian_K_block_gate",
            "peter_weyl_lorentzian_ordered_triple_gate",
            "peter_weyl_lorentzian_sine_ordered_triple_gate",
            "peter_weyl_q_block_gauge_audit",
            "peter_weyl_volume_cross_rep_audit",
            "peter_weyl_zeroaware_volume_migration_experiment",
        }
        self.assertEqual(modules, expected_modules)

        plain_lru_caches = []
        for module_name in sorted(modules):
            source_path = scripts / (module_name + ".py")
            tree = ast.parse(source_path.read_text(encoding="utf-8"))
            for function in (
                node for node in ast.walk(tree)
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
            ):
                # This is the private primitive inside audit_lru_cache itself;
                # the public wrapper explicitly bypasses it when cache_mode=off.
                if module_name == "bqg_pruning_audit" and function.name == "cached":
                    continue
                for decorator in function.decorator_list:
                    name = ast.unparse(decorator)
                    if "lru_cache" in name and "audit_lru_cache" not in name:
                        plain_lru_caches.append((module_name, function.name, name))
        self.assertEqual(plain_lru_caches, [])

        # Every sparse-amplitude decision must emit an event directly, use
        # the instrumented sparse-state filter, or write V10's own ledger row.
        unrecorded_decisions = []
        for module_name in sorted(modules):
            tree = ast.parse((scripts / (module_name + ".py")).read_text(encoding="utf-8"))
            for function in (
                node for node in ast.walk(tree)
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
            ):
                call_names = {
                    ast.unparse(node.func)
                    for node in ast.walk(function) if isinstance(node, ast.Call)
                }
                has_decision = any(
                    name.endswith(("PRUNE.keep", "PRUNE.filter_state", "prune_state"))
                    for name in call_names
                )
                has_audit = any(
                    name.endswith(("PRUNE.record_threshold", "PRUNE.record_merge",
                                   "PRUNE.filter_state", "prune_state", "record_dropped"))
                    for name in call_names
                ) or (module_name == "bqg_v10_w5_lorentzian_adjoint_closure"
                      and function.name == "record_dropped")
                if has_decision and not has_audit:
                    unrecorded_decisions.append((module_name, function.name, function.lineno))
        self.assertEqual(unrecorded_decisions, [])

        # The shared sparse-state filter must itself record every nonzero
        # coefficient at/below threshold before applying PRUNE.keep.
        pw_tree = ast.parse((scripts / "k5_peter_weyl_safe_hda_column.py").read_text(encoding="utf-8"))
        prune_defs = [
            node for node in ast.walk(pw_tree)
            if isinstance(node, ast.FunctionDef) and node.name == "prune_state"
        ]
        self.assertEqual(len(prune_defs), 1)
        prune_calls = {
            ast.unparse(node.func)
            for node in ast.walk(prune_defs[0]) if isinstance(node, ast.Call)
        }
        self.assertIn("PRUNE.record_threshold", prune_calls)
        self.assertIn("PRUNE.keep", prune_calls)

        audit_tree = ast.parse((scripts / "bqg_pruning_audit.py").read_text(encoding="utf-8"))
        filter_defs = [
            node for node in ast.walk(audit_tree)
            if isinstance(node, ast.FunctionDef) and node.name == "filter_state"
        ]
        self.assertEqual(len(filter_defs), 1)
        filter_calls = {
            ast.unparse(node.func)
            for node in ast.walk(filter_defs[0]) if isinstance(node, ast.Call)
        }
        self.assertIn("record_threshold", filter_calls)
        self.assertIn("keep", filter_calls)

        # V10's sparse accumulator is allowed to be uninstrumented only
        # because its default tolerance is exactly zero (exact-zero removal).
        v10_path = scripts / "bqg_v10_w5_lorentzian_adjoint_closure.py"
        v10_text = v10_path.read_text(encoding="utf-8")
        v10_tree = ast.parse(v10_text)
        tol_assignments = [
            node.value for node in ast.walk(v10_tree)
            if isinstance(node, ast.Assign)
            and any(isinstance(target, ast.Name) and target.id == "TOL"
                    for target in node.targets)
        ]
        self.assertTrue(any(
            isinstance(value, ast.Constant) and value.value == 0.0
            for value in tol_assignments
        ))
        add_sparse_defs = [
            node for node in ast.walk(v10_tree)
            if isinstance(node, ast.FunctionDef) and node.name == "add_sparse"
        ]
        self.assertEqual(len(add_sparse_defs), 1)
        self.assertEqual(ast.unparse(add_sparse_defs[0].args.defaults[-1]), "TOL")
        add_sparse_calls = [
            node for node in ast.walk(v10_tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name) and node.func.id == "add_sparse"
        ]
        self.assertTrue(add_sparse_calls)
        self.assertTrue(all(
            len(node.args) <= 3 and all(keyword.arg != "tol" for keyword in node.keywords)
            for node in add_sparse_calls
        ))

        # These are deliberately not amplitude pruning:
        # - add_sparse has tol=0 at every V10 call site: exact-zero cancellation;
        # - spin probes test exact nonzero branch existence before wall accounting;
        # - source-J and charged-label bounds are exact SU(2) selection rules;
        # - zero-aware eigenvalue masks define the volume operator and are a
        #   separate spectral sensitivity axis, not sparse-amplitude deletion.
        exact_or_separate = {
            ("bqg_sparse_linear_accumulator", "sparsify_matrix_columns"),
            ("bqg_sparse_linear_accumulator", "sparsify_matrix_columns_with_error_bound"),
            ("bqg_v10_w5_lorentzian_adjoint_closure", "add_sparse"),
            ("bqg_v10_w5_lorentzian_adjoint_closure", "_make_spin_wall_probe"),
            ("bqg_v10_w5_lorentzian_adjoint_closure", "probed"),
            ("peter_weyl_covariant_volume_leg_gate", "run"),
            ("peter_weyl_lorentzian_ordered_triple_gate", "run"),
            ("charged_intertwiner_recoupling_gate", "allowed_charged_labels"),
            ("peter_weyl_zeroaware_volume_migration_experiment", "zeroaware_sqrt_abs"),
            ("peter_weyl_volume_cross_rep_audit", "stable_zeroaware_volume"),
            ("peter_weyl_q_block_gauge_audit", "sqrt_abs"),
            ("bqg_pruning_audit", "keep"),
            ("bqg_pruning_audit", "record_threshold"),
            ("bqg_pruning_audit", "record_merge"),
        }

        unclassified = []
        for module_name in sorted(modules):
            source_path = scripts / (module_name + ".py")
            source_text = source_path.read_text(encoding="utf-8")
            tree = ast.parse(source_text)
            for function in (
                node for node in ast.walk(tree)
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
            ):
                function_text = ast.get_source_segment(source_text, function) or ""
                threshold_comparisons = []
                for node in ast.walk(function):
                    if not isinstance(node, ast.Compare):
                        continue
                    left = ast.unparse(node.left)
                    if ("abs(" in left or "np.abs(" in left):
                        threshold_comparisons.append(ast.unparse(node))
                if (
                    threshold_comparisons
                    and "PRUNE.keep" not in function_text
                    and "prune_state" not in function_text
                    and (module_name, function.name) not in exact_or_separate
                ):
                    unclassified.append(
                        (module_name, function.name, function.lineno, threshold_comparisons)
                    )

        self.assertEqual(unclassified, [])


class TestChargedIntertwinerTensorContraction(unittest.TestCase):
    @staticmethod
    def reference_tensor(spins, K12, K34, J2, M2):
        spins = tuple(spins)
        dims = tuple(s + 1 for s in spins)
        tensor = np.zeros(dims, complex)
        magnetic = [PW.m2vals_t(s) for s in spins]
        for indices in itertools.product(*[range(d) for d in dims]):
            m1, m2, m3, m4 = [magnetic[a][indices[a]] for a in range(4)]
            value = 0.0
            for MK in PW.m2vals_t(K12):
                ML = M2 - MK
                if ML not in PW.m2vals_t(K34):
                    continue
                value += (
                    PW.cg2(spins[0], spins[1], K12, m1, m2, MK)
                    * PW.cg2(spins[2], spins[3], K34, m3, m4, ML)
                    * PW.cg2(K12, K34, J2, MK, ML, M2)
                )
            tensor[indices] = value
        return tensor

    def test_vectorized_basis_projection_matches_original_vdot_sum(self):
        rng = np.random.default_rng(20261113)
        basis = [rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2)) for _ in range(4)]
        X = rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2))
        expected_coeff = np.asarray([np.vdot(B, X) for B in basis], complex)
        expected_recon = sum(
            (coefficient * B for coefficient, B in zip(expected_coeff, basis)),
            np.zeros_like(X),
        )
        actual_coeff, actual_recon = CHARGED.project_to_basis(X, basis)
        np.testing.assert_allclose(actual_coeff, expected_coeff, rtol=0.0, atol=1e-14)
        np.testing.assert_allclose(actual_recon, expected_recon, rtol=0.0, atol=1e-14)

    def test_vectorized_complete_multiplet_projection_matches_explicit_vdot_loop(self):
        rng = np.random.default_rng(20261114)
        for spins in ((1, 1, 1, 1), (2, 1, 1, 2)):
            dims = tuple(s + 1 for s in spins)
            X = rng.normal(size=dims) + 1j * rng.normal(size=dims)
            for J2 in CV.all_total_J2(spins):
                matrix, magnetic_count, channel_count = CHARGED.charged_basis_multiplet_matrix(
                    spins, J2
                )
                if channel_count == 0:
                    continue
                expected_coefficients = []
                expected_reconstruction = np.zeros_like(X)
                for M2 in PW.m2vals_t(J2):
                    for K12, K34 in CHARGED.allowed_charged_labels(spins, J2):
                        basis = CHARGED.charged_tensor(spins, K12, K34, J2, M2)
                        coefficient = np.vdot(basis, X)
                        expected_coefficients.append(coefficient)
                        expected_reconstruction += coefficient * basis
                actual_coefficients = matrix.conj() @ X.reshape(-1)
                actual_reconstruction = (actual_coefficients @ matrix).reshape(X.shape)
                np.testing.assert_allclose(
                    actual_coefficients, expected_coefficients, rtol=2e-14, atol=2e-13
                )
                np.testing.assert_allclose(
                    actual_reconstruction, expected_reconstruction, rtol=2e-13, atol=2e-12
                )
                self.assertEqual(magnetic_count, len(PW.m2vals_t(J2)))

    def test_cached_multiplet_basis_stacks_every_M_sector(self):
        for spins in ((1, 1, 1, 1), (4, 4, 4, 4), (2, 4, 4, 2)):
            for J2 in CV.all_total_J2(spins):
                expected = np.concatenate([
                    CHARGED.charged_basis_matrix(spins, J2, M2)
                    for M2 in PW.m2vals_t(J2)
                ], axis=0)
                actual, magnetic_count, channel_count = CHARGED.charged_basis_multiplet_matrix(spins, J2)
                np.testing.assert_allclose(actual, expected, rtol=0.0, atol=2e-14)
                self.assertEqual(magnetic_count, J2 + 1)
                self.assertEqual(channel_count, len(CHARGED.allowed_charged_labels(spins, J2)))
                self.assertFalse(actual.flags.writeable)
        CHARGED.charged_basis_multiplet_matrix.cache_clear()

    def test_batched_all_M_volume_action_matches_original_sector_loop(self):
        spins = (1, 1, 1, 1)
        rng = np.random.default_rng(20261113)
        T = rng.normal(size=(2, 2, 2, 2)) + 1j * rng.normal(size=(2, 2, 2, 2))
        X = CV.unorient_local(T, spins, 0)
        expected_Y = np.zeros_like(X)
        expected_recon = np.zeros_like(X)
        for J2 in CV.all_total_J2(spins):
            Vb = CK.canonical_volume_block_general(spins, J2)
            labels = CHARGED.allowed_charged_labels(spins, J2)
            for M2 in PW.m2vals_t(J2):
                basis = [
                    CHARGED.charged_tensor(spins, a, b, J2, M2)
                    for a, b in labels
                ]
                coeff = np.asarray([np.vdot(B, X) for B in basis], complex)
                for coefficient, B in zip(coeff, basis):
                    expected_recon += coefficient * B
                output = Vb @ coeff
                for coefficient, B in zip(output, basis):
                    expected_Y += coefficient * B
        expected_leak = float(
            np.linalg.norm(X - expected_recon) / max(np.linalg.norm(X), 1e-30)
        )
        expected = CV.orient_local(expected_Y, spins, 0)
        actual, actual_leak = CK.apply_volume_allJ_oriented(T, spins, 0)
        np.testing.assert_allclose(actual, expected, rtol=0.0, atol=1e-12)
        self.assertAlmostEqual(actual_leak, expected_leak, delta=1e-14)

    def test_cached_flattened_basis_matches_charged_tensor_stack(self):
        spins = (1, 1, 1, 1)
        for J2 in CV.all_total_J2(spins):
            for M2 in PW.m2vals_t(J2):
                labels = CHARGED.allowed_charged_labels(spins, J2)
                expected = np.stack([
                    CHARGED.charged_tensor(spins, a, b, J2, M2).reshape(-1)
                    for a, b in labels
                ], axis=0)
                actual = CHARGED.charged_basis_matrix(spins, J2, M2)
                np.testing.assert_allclose(actual, expected, rtol=0.0, atol=1e-14)
                self.assertFalse(actual.flags.writeable)
        CHARGED.charged_basis_matrix.cache_clear()

    def test_batched_q_block_matches_original_per_basis_action(self):
        cases = (((4, 4, 4, 4), 0, 0), ((4, 4, 4, 4), 4, 4), ((2, 4, 4, 2), 2, 2))
        for spins, J2, M2 in cases:
            with self.subTest(spins=spins, J2=J2, M2=M2):
                basis = [
                    CHARGED.charged_tensor(spins, a, b, J2, M2)
                    for a, b in CHARGED.allowed_charged_labels(spins, J2)
                ]
                expected = np.zeros((len(basis), len(basis)), complex)
                abs_leak = rel_leak = 0.0
                for j, vector in enumerate(basis):
                    q_vector = CHARGED.apply_q_four(vector, spins)
                    coeff = np.asarray([np.vdot(item, q_vector) for item in basis], complex)
                    expected[:, j] = coeff
                    recon = sum(
                        (value * item for value, item in zip(coeff, basis)),
                        np.zeros_like(q_vector),
                    )
                    local_abs = float(np.linalg.norm(q_vector - recon))
                    local_rel = local_abs / max(float(np.linalg.norm(q_vector)), 1e-30)
                    abs_leak = max(abs_leak, local_abs)
                    rel_leak = max(rel_leak, local_rel)
                expected = 0.5 * (expected + expected.conj().T)
                _, actual, actual_abs, actual_rel = CHARGED.q_block(spins, J2, M2)
                np.testing.assert_allclose(actual, expected, rtol=0.0, atol=1e-12)
                self.assertAlmostEqual(actual_abs, abs_leak, delta=1e-12)
                self.assertAlmostEqual(actual_rel, rel_leak, delta=1e-12)

    def test_highest_weight_volume_block_matches_all_magnetic_sectors(self):
        spins = (1, 1, 1, 1)
        for J2 in CV.all_total_J2(spins):
            blocks = [
                CHARGED.q_block(spins, J2, M2)[1]
                for M2 in PW.m2vals_t(J2)
            ]
            for block in blocks[1:]:
                np.testing.assert_allclose(block, blocks[0], rtol=0.0, atol=1e-12)

            averaged = sum(blocks) / len(blocks)
            expected = CHARGED.canonical_volume_block(averaged)
            actual = CK.canonical_volume_block_general(spins, J2)
            np.testing.assert_allclose(actual, expected, rtol=0.0, atol=1e-12)
            CK.canonical_volume_block_general.cache_clear()

    def test_separable_contraction_matches_original_magnetic_sum(self):
        cases = (
            ((4, 4, 4, 4), 4, 4, 0, 0),
            ((4, 4, 4, 4), 2, 2, 2, 0),
            ((4, 4, 4, 4), 4, 2, 2, 0),
            ((2, 4, 4, 2), 2, 4, 2, 0),
            ((4, 2, 4, 2), 2, 2, 2, 2),
        )
        for case in cases:
            with self.subTest(case=case):
                CHARGED.charged_tensor.cache_clear()
                actual = CHARGED.charged_tensor(*case)
                expected = self.reference_tensor(*case)
                np.testing.assert_allclose(actual, expected, rtol=0.0, atol=1e-14)
                self.assertEqual(actual.shape, expected.shape)
                with audit.cache_bypass():
                    uncached = CHARGED.charged_tensor.__wrapped__(*case)
                np.testing.assert_allclose(uncached, expected, rtol=0.0, atol=1e-14)
        CHARGED.charged_tensor.cache_clear()


class TestW5CacheTransparency(unittest.TestCase):
    def test_context_local_cache_bypass_reexecutes_operator_and_keeps_audit(self):
        calls = []

        @audit.audit_lru_cache(None)
        def operator(value):
            calls.append(value)
            audit.record_threshold("cache_bypass.fixture", value, 0.5e-12, 1e-12)
            return value * 3

        with audit.capture() as collector:
            self.assertEqual(operator(7), 21)
            self.assertEqual(operator(7), 21)
            with audit.cache_bypass():
                self.assertEqual(operator(7), 21)
                self.assertEqual(operator(7), 21)

        self.assertEqual(calls, [7, 7, 7])
        rows = collector.snapshot()["events"]
        self.assertEqual(len(rows), 4)
        self.assertEqual(
            [row["kind"] for row in rows],
            ["coefficient_threshold", "coefficient_threshold",
             "coefficient_threshold", "coefficient_threshold"],
        )
        self.assertEqual(operator.cache_info().hits, 1)
        self.assertEqual(operator.cache_info().misses, 1)
        operator.cache_clear()

    def test_gauss_volume_column_cached_and_uncached_full_outputs_match(self):
        import peter_weyl_lorentzian_K_block_gate as KG

        key = PW.basis_full_jhalf()[0]
        KG.local_volume_column.cache_clear()
        with audit.capture() as cold:
            uncached = KG.local_volume_column.__wrapped__(key, 0)
        with audit.capture() as warm_miss:
            cached = KG.local_volume_column(key, 0)
        with audit.capture() as warm_hit:
            cached_again = KG.local_volume_column(key, 0)

        self.assertEqual(cached, uncached)
        self.assertEqual(cached_again, uncached)
        self.assertEqual(cold.snapshot()["events"], warm_miss.snapshot()["events"])
        self.assertEqual(warm_miss.snapshot()["events"], warm_hit.snapshot()["events"])
        self.assertEqual(cold.snapshot()["event_count"], 2)
        self.assertEqual(warm_hit.snapshot()["event_count"], 2)
        self.assertEqual(
            [row["site_id"] for row in warm_hit.snapshot()["events"]],
            [
                "lorentzian_K_block.local_volume.1e-12",
                "lorentzian_K_block.local_volume.propagated_loss",
            ],
        )
        info = KG.local_volume_column.cache_info()
        self.assertEqual(info.misses, 1)
        self.assertEqual(info.hits, 1)
        KG.local_volume_column.cache_clear()

    def test_local_volume_loss_is_scaled_by_actual_incoming_branch_amplitude(self):
        import peter_weyl_lorentzian_K_block_gate as KG
        from unittest.mock import patch

        key = PW.basis_full_jhalf()[0]
        output_key = ("omitted-volume-output",)
        coefficient = 0.5e-12 - 0.25e-12j
        incoming = -2.0 + 3.0j
        with patch.object(
            KG, "_local_volume_column_full",
            return_value=((), ((output_key, coefficient),)),
        ), patch.object(PW, "prune_state", side_effect=lambda state, tol: state):
            with audit.capture() as collector:
                result = KG.apply_V_local({key: incoming}, 0)

        self.assertEqual(result, {})
        event = collector.snapshot()["events"][0]
        propagated = incoming * coefficient
        self.assertEqual(event["site_id"], "lorentzian_K_block.local_volume.propagated_loss")
        self.assertAlmostEqual(event["incoming_real"], incoming.real)
        self.assertAlmostEqual(event["incoming_imag"], incoming.imag)
        self.assertAlmostEqual(event["propagated_real"], propagated.real)
        self.assertAlmostEqual(event["propagated_imag"], propagated.imag)
        self.assertAlmostEqual(event["propagated_abs"], abs(propagated))

    def test_complete_HE_cache_preserves_full_complex_state_and_diagnostics(self):
        import peter_weyl_covariant_K_leg_gate as KLEG
        from unittest.mock import patch

        key = ("cache-fixture", (1 + 2j, -3 + 0.5j))
        expected_state = {
            ("out-a",): 0.125 - 0.75j,
            ("out-b",): -2.0 + 1e-14j,
        }
        expected = (tuple(expected_state.items()), 0.125, 0.375)

        def deterministic_operator(input_key, source_v, Jmax2):
            self.assertEqual(input_key, key)
            self.assertEqual(source_v, 2)
            self.assertEqual(Jmax2, 10)
            return expected_state, 0.125, 0.375

        KLEG.HE_complete_cached.cache_clear()
        with patch.object(KLEG, "apply_HE_complete_key", side_effect=deterministic_operator):
            direct = KLEG.HE_complete_cached.__wrapped__(key, 2, 10)
            cached = KLEG.HE_complete_cached(key, 2, 10)
            cached_again = KLEG.HE_complete_cached(key, 2, 10)

        self.assertEqual(direct, expected)
        self.assertEqual(cached, direct)
        self.assertEqual(cached_again, direct)
        info = KLEG.HE_complete_cached.cache_info()
        self.assertEqual(info.misses, 1)
        self.assertEqual(info.hits, 1)
        KLEG.HE_complete_cached.cache_clear()


    def test_complete_HE_cache_replays_nested_audit_trace_on_warm_hit(self):
        from unittest.mock import patch
        import peter_weyl_covariant_K_leg_gate as KLEG

        key = ("audit-cache-fixture", (1 + 2j,))
        expected = ({("out",): 0.25 - 0.5j}, 0.0, 0.0)

        def deterministic_operator(input_key, source_v, Jmax2):
            audit.record_threshold("nested.fixture.threshold", input_key, 0.5e-12, 1e-12)
            return expected

        KLEG.HE_complete_cached.cache_clear()
        with patch.object(KLEG, "apply_HE_complete_key", side_effect=deterministic_operator):
            with audit.capture() as cold:
                first = KLEG.HE_complete_cached(key, 1, 10)
            with audit.capture() as warm:
                second = KLEG.HE_complete_cached(key, 1, 10)

        self.assertEqual(first, second)
        self.assertEqual(cold.snapshot()["events"], warm.snapshot()["events"])
        self.assertEqual(cold.snapshot()["event_count"], 1)
        self.assertEqual(warm.snapshot()["event_count"], 1)
        self.assertEqual(KLEG.HE_complete_cached.cache_info().hits, 1)
        KLEG.HE_complete_cached.cache_clear()

    def test_audit_lru_cache_replays_nested_and_bypass_event_kinds(self):
        @audit.audit_lru_cache(None)
        def cached_operator(key):
            audit.record_threshold("nested", key, 0.25e-9, 1e-9)
            return {key: 1 + 2j}

        with audit.capture() as cold:
            normal = cached_operator("x")
        with audit.capture() as warm:
            normal_again = cached_operator("x")
        self.assertEqual(normal, normal_again)
        self.assertEqual(cold.snapshot()["events"], warm.snapshot()["events"])

        with audit.capture() as shadow_mode:
            with audit.bypass():
                shadow = cached_operator("x")
        with audit.capture() as shadow_warm:
            with audit.bypass():
                shadow_again = cached_operator("x")
        self.assertEqual(shadow, shadow_again)
        self.assertEqual(shadow_mode.snapshot()["events"], shadow_warm.snapshot()["events"])
        self.assertEqual(shadow_mode.snapshot()["events"][0]["kind"], "would_prune_bypassed")
        self.assertEqual(cached_operator.cache_info().misses, 2)
        self.assertEqual(cached_operator.cache_info().hits, 2)


    def test_nested_cached_operator_wrappers_replay_pruning_audit(self):
        scripts = Path(__file__).resolve().parents[1] / "scripts"
        wrappers = {
            "peter_weyl_covariant_K_composition_cached_gate.py": 5,
            "peter_weyl_lorentzian_sine_ordered_triple_gate.py": 7,
        }
        for filename, expected in wrappers.items():
            source = (scripts / filename).read_text(encoding="utf-8")
            self.assertIn("import bqg_pruning_audit as PRUNE", source)
            self.assertEqual(
                source.count("@PRUNE.audit_lru_cache(None)"), expected,
                msg=f"every nested numerical cache in {filename} must replay audit rows",
            )
            self.assertNotIn("@functools.lru_cache(maxsize=None)", source)



if __name__ == "__main__":
    unittest.main()
