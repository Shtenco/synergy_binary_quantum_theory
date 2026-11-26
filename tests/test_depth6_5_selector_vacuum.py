import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "depth6_5_selector_vacuum",
    ROOT / "scripts" / "audit_depth6_5_selector_vacuum.py",
)
assert SPEC and SPEC.loader
PROBE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = PROBE
SPEC.loader.exec_module(PROBE)


def test_exact_trivial_projector_path_reaches_generic_master_row_on_vacuum():
    result = PROBE.run_probe()
    assert result["status"] == "EXACT_TRIVIAL_SELECTOR_VACUUM_ROW_ZERO"
    assert result["selector"]["dimension"] == 1
    assert result["selector"]["max_deviation_from_identity"] == 0
    assert result["generic_master_row"]["output_support"] == 0
    assert result["generic_master_row"]["output_norm"] == 0
