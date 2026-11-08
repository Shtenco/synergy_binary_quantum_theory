import unittest


def validate_path_metadata(meta):
    """Isolated copy of the strict path-metadata guard in the V21 collector."""
    if not isinstance(meta, dict):
        raise ValueError("INVALID_PATH_METADATA")
    mode, order, aux = meta.get("mode"), meta.get("order"), meta.get("aux")
    if (mode not in ("forward", "adjoint") or order not in ("p", "q")
            or not isinstance(aux, list) or len(aux) != 3
            or any(type(value) is not int or value not in (0, 1) for value in aux)):
        raise ValueError("INVALID_PATH_INDEX")
    for key in ("source_node", "omit", "cycle", "Jmax2"):
        if type(meta.get(key)) is not int:
            raise ValueError("BAD_STAGE_INTEGER_METADATA_" + key.upper())
    return mode, order, tuple(aux)


def validate_manifest_indices(manifest):
    """Isolated copy of the strict top-level manifest selector guard."""
    if not isinstance(manifest, dict):
        raise ValueError("BAD_V21_MANIFEST")
    values = []
    for key in ("source_node", "omit", "cycle", "Jmax2"):
        value = manifest.get(key)
        if type(value) is not int:
            raise ValueError("BAD_MANIFEST_INTEGER_" + key.upper())
        values.append(value)
    return tuple(values[:3])


def validate_support(value, predecessor=False):
    if type(value) is not int or value <= (0 if predecessor else -1):
        raise ValueError("BAD_PREDECESSOR_STATE_SUPPORT" if predecessor
                         else "BAD_STATE_SUPPORT_METADATA")
    return value


class MetadataGuardAudit(unittest.TestCase):
    def test_boolean_aux_is_rejected(self):
        meta = {"mode": "forward", "order": "p", "aux": [False, 0, 1],
                "source_node": 0, "omit": 0, "cycle": 0, "Jmax2": 7}
        with self.assertRaisesRegex(ValueError, "INVALID_PATH_INDEX"):
            validate_path_metadata(meta)

    def test_float_aux_is_rejected(self):
        meta = {"mode": "forward", "order": "p", "aux": [0, 1.0, 1],
                "source_node": 0, "omit": 0, "cycle": 0, "Jmax2": 7}
        with self.assertRaisesRegex(ValueError, "INVALID_PATH_INDEX"):
            validate_path_metadata(meta)

    def test_boolean_manifest_selector_is_rejected(self):
        manifest = {"source_node": True, "omit": 0, "cycle": 0, "Jmax2": 7}
        with self.assertRaisesRegex(ValueError, "BAD_MANIFEST_INTEGER_SOURCE_NODE"):
            validate_manifest_indices(manifest)

    def test_boolean_checkpoint_support_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "BAD_STATE_SUPPORT_METADATA"):
            validate_support(True)
        with self.assertRaisesRegex(ValueError, "BAD_PREDECESSOR_STATE_SUPPORT"):
            validate_support(True, predecessor=True)

    def test_zero_predecessor_support_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "BAD_PREDECESSOR_STATE_SUPPORT"):
            validate_support(0, predecessor=True)

    def test_canonical_metadata_is_accepted(self):
        meta = {"mode": "adjoint", "order": "q", "aux": [1, 0, 1],
                "source_node": 0, "omit": 2, "cycle": 1, "Jmax2": 7}
        self.assertEqual(validate_path_metadata(meta), ("adjoint", "q", (1, 0, 1)))
        self.assertEqual(validate_manifest_indices(meta), (0, 2, 1))
        self.assertEqual(validate_support(1), 1)
        self.assertEqual(validate_support(1, predecessor=True), 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
