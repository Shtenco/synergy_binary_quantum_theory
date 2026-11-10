import io
import stat
import tarfile
import tempfile
import unittest
import zipfile
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "proof_tools"))
from safe_archive import safe_artifact_directory, safe_extract_tar, safe_extract_zip


class SafeArchiveTests(unittest.TestCase):
    def test_extracts_regular_tar_members(self):
        with tempfile.TemporaryDirectory() as tmp:
            archive_bytes = io.BytesIO()
            with tarfile.open(fileobj=archive_bytes, mode="w:gz") as tf:
                payload = b"trusted"
                info = tarfile.TarInfo("subdir/data.bin")
                info.size = len(payload)
                tf.addfile(info, io.BytesIO(payload))
            archive_bytes.seek(0)
            dest = Path(tmp) / "out"
            with tarfile.open(fileobj=archive_bytes, mode="r:gz") as tf:
                safe_extract_tar(tf, dest)
            self.assertEqual((dest / "subdir" / "data.bin").read_bytes(), b"trusted")

    def test_rejects_tar_parent_traversal(self):
        with tempfile.TemporaryDirectory() as tmp:
            archive_bytes = io.BytesIO()
            with tarfile.open(fileobj=archive_bytes, mode="w:gz") as tf:
                payload = b"escape"
                info = tarfile.TarInfo("../escape.txt")
                info.size = len(payload)
                tf.addfile(info, io.BytesIO(payload))
            archive_bytes.seek(0)
            with tarfile.open(fileobj=archive_bytes, mode="r:gz") as tf:
                with self.assertRaisesRegex(ValueError, "escapes extraction root"):
                    safe_extract_tar(tf, Path(tmp) / "out")
            self.assertFalse((Path(tmp) / "escape.txt").exists())

    def test_rejects_tar_symlink(self):
        with tempfile.TemporaryDirectory() as tmp:
            archive_bytes = io.BytesIO()
            with tarfile.open(fileobj=archive_bytes, mode="w:gz") as tf:
                link = tarfile.TarInfo("link")
                link.type = tarfile.SYMTYPE
                link.linkname = "../outside"
                tf.addfile(link)
            archive_bytes.seek(0)
            with tarfile.open(fileobj=archive_bytes, mode="r:gz") as tf:
                with self.assertRaisesRegex(ValueError, "non-regular tar member"):
                    safe_extract_tar(tf, Path(tmp) / "out")

    def test_extracts_regular_zip_members(self):
        with tempfile.TemporaryDirectory() as tmp:
            archive_bytes = io.BytesIO()
            with zipfile.ZipFile(archive_bytes, "w") as zf:
                zf.writestr("subdir/data.bin", b"trusted")
            archive_bytes.seek(0)
            dest = Path(tmp) / "out"
            with zipfile.ZipFile(archive_bytes) as zf:
                safe_extract_zip(zf, dest)
            self.assertEqual((dest / "subdir" / "data.bin").read_bytes(), b"trusted")

    def test_rejects_zip_parent_traversal(self):
        with tempfile.TemporaryDirectory() as tmp:
            archive_bytes = io.BytesIO()
            with zipfile.ZipFile(archive_bytes, "w") as zf:
                zf.writestr("../escape.txt", b"escape")
            archive_bytes.seek(0)
            with zipfile.ZipFile(archive_bytes) as zf:
                with self.assertRaisesRegex(ValueError, "escapes extraction root"):
                    safe_extract_zip(zf, Path(tmp) / "out")
            self.assertFalse((Path(tmp) / "escape.txt").exists())

    def test_rejects_zip_symlink(self):
        with tempfile.TemporaryDirectory() as tmp:
            archive_bytes = io.BytesIO()
            info = zipfile.ZipInfo("link")
            info.create_system = 3
            info.external_attr = (stat.S_IFLNK | 0o777) << 16
            with zipfile.ZipFile(archive_bytes, "w") as zf:
                zf.writestr(info, "../outside")
            archive_bytes.seek(0)
            with zipfile.ZipFile(archive_bytes) as zf:
                with self.assertRaisesRegex(ValueError, "non-regular ZIP member"):
                    safe_extract_zip(zf, Path(tmp) / "out")

    def test_rejects_nonempty_destination(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp) / "out"
            dest.mkdir()
            (dest / "stale.txt").write_text("stale")
            with zipfile.ZipFile(io.BytesIO(), "w") as zf:
                with self.assertRaisesRegex(ValueError, "not empty"):
                    safe_extract_zip(zf, dest)

    def test_artifact_name_must_be_single_safe_component(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(
                safe_artifact_directory(tmp, "bqg-proof-shard-01").name,
                "bqg-proof-shard-01",
            )
            for name in ("../escape", "/absolute", "a\\..\\b", ".", "..", ""):
                with self.subTest(name=name):
                    with self.assertRaisesRegex(ValueError, "unsafe artifact name"):
                        safe_artifact_directory(tmp, name)


if __name__ == "__main__":
    unittest.main()
