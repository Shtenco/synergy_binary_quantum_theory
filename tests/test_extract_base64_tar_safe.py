import base64
import io
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from extract_base64_tar_safe import extract_base64_tar


def make_bundle(name: str, payload: bytes) -> bytes:
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        info = tarfile.TarInfo(name)
        info.size = len(payload)
        archive.addfile(info, io.BytesIO(payload))
    return base64.b64encode(buffer.getvalue())


class Base64TarExtractionTests(unittest.TestCase):
    def write_bundle(self, root: Path, name: str, encoded: bytes) -> Path:
        path = root / name
        path.write_bytes(encoded)
        return path

    def test_extracts_regular_file(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive = self.write_bundle(root, "safe.b64", make_bundle("pkg/mod.py", b"x = 1\n"))
            destination = root / "out"
            extract_base64_tar(archive, destination)
            self.assertEqual((destination / "pkg/mod.py").read_bytes(), b"x = 1\n")

    def test_rejects_parent_traversal(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive = self.write_bundle(root, "bad.b64", make_bundle("../escape.py", b"x = 1\n"))
            with self.assertRaisesRegex(ValueError, "escapes extraction root"):
                extract_base64_tar(archive, root / "out")
            self.assertFalse((root / "escape.py").exists())

    def test_rejects_uncompressed_size_bomb(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive = self.write_bundle(root, "large.b64", make_bundle("large.bin", b"x" * 128))
            with self.assertRaisesRegex(ValueError, "uncompressed archive exceeds limit"):
                extract_base64_tar(archive, root / "out", max_uncompressed_bytes=64)

    def test_rejects_invalid_base64(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive = self.write_bundle(root, "invalid.b64", b"not base64!")
            with self.assertRaises(Exception):
                extract_base64_tar(archive, root / "out")


if __name__ == "__main__":
    unittest.main()
