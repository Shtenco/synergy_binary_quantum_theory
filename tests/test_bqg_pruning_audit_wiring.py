import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import bqg_pruning_audit as audit
import k5_peter_weyl_safe_hda_column as PW
import peter_weyl_covariant_K_leg_gate as CK
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
            ("bqg_v10_w5_lorentzian_adjoint_closure", "add_sparse"),
            ("bqg_v10_w5_lorentzian_adjoint_closure", "_make_spin_wall_probe"),
            ("bqg_v10_w5_lorentzian_adjoint_closure", "probed"),
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


class TestW5CacheTransparency(unittest.TestCase):
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
        info = KG.local_volume_column.cache_info()
        self.assertEqual(info.misses, 1)
        self.assertEqual(info.hits, 1)
        KG.local_volume_column.cache_clear()

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



if __name__ == "__main__":
    unittest.main()
