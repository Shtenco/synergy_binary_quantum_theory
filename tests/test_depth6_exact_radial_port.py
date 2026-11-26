import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "depth6_collapsed_exact_port",
    ROOT / "scripts" / "verify_depth6_vacuum_collapsed_exact_equivalence.py",
)
assert SPEC and SPEC.loader
PORT = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = PORT
SPEC.loader.exec_module(PORT)


def test_exact_far_radial_weights_match_pinned_source():
    for s, so in ((0, 1), (1, 0), (1, 2), (2, 1)):
        for v, c in ((0, 4), (4, 0)):
            exact = PORT.exact_far_radial_weight(s, so, v, c)
            source = PORT.A.far_radial_weight(s, so, v, c)
            assert abs(complex(exact.evalf()) - source) < 1e-14
