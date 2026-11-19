import json
import numpy as np
import pytest

from scripts.q2_g3_matrix_free_feshbach import run


def test_g3_feshbach_refuses_non_g3_export(tmp_path):
    np.save(tmp_path / "g3_dual_endpoints.npy",
            np.array([[0, 1], [1, 2]], dtype=np.uint32))
    np.save(tmp_path / "g3_fibre_labels.npy",
            np.array([0, 0, 1], dtype=np.uint16))
    (tmp_path / "G3_CHUNKED_MANIFEST.json").write_text(
        json.dumps({"passed": True}), encoding="utf-8"
    )
    with pytest.raises(ValueError, match="unexpected G3 graph dimensions"):
        run(tmp_path)


def test_g3_feshbach_module_has_explicit_probe_scope():
    from scripts import q2_g3_matrix_free_feshbach as module

    assert "not the full 512x512 self-energy matrix" in module.__doc__


def test_feshbach_report_is_valid_json_with_real_line_ending(tmp_path):
    from scripts.q2_g3_matrix_free_feshbach import write_report

    target = tmp_path / "report.json"
    write_report(target, {"passed": True, "generation": 3})
    raw = target.read_text(encoding="utf-8")
    assert raw.endswith("\n")
    assert json.loads(raw) == {"passed": True, "generation": 3}
