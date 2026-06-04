from __future__ import annotations

import io
import os
import tempfile
import unittest
import sys
import types
from pathlib import Path
from unittest import mock

from downstream_tasks.downloader import DatasetDownloader
from downstream_tasks.utils.files import sha256sum, sidecar_metadata_path, read_json


class _FakeResponse:
    def __init__(self, *, status_code: int, body: bytes, headers: dict | None = None):
        self.status_code = status_code
        self._body = body
        self.headers = headers or {}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def iter_content(self, chunk_size: int = 1024 * 1024):
        for i in range(0, len(self._body), chunk_size):
            yield self._body[i : i + chunk_size]


class TestDatasetDownloader(unittest.TestCase):
    def test_download_writes_file_and_metadata(self):
        d = DatasetDownloader()
        payload = b"hello-world"

        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "file.bin"

            def _get(url, headers=None, stream=True, timeout=1):
                return _FakeResponse(status_code=200, body=payload, headers={})

            fake_requests = types.SimpleNamespace(get=_get)
            with mock.patch.dict(sys.modules, {"requests": fake_requests}):
                res = d.download(url="https://example.com/file.bin", output_path=out, timeout_s=1, max_retries=1)

            self.assertTrue(out.exists())
            self.assertEqual(out.read_bytes(), payload)
            self.assertEqual(res.bytes_downloaded, len(payload))
            self.assertEqual(res.sha256, sha256sum(out))

            meta_path = sidecar_metadata_path(out)
            self.assertTrue(meta_path.exists())
            meta = read_json(meta_path)
            self.assertEqual(meta["url"], "https://example.com/file.bin")
            self.assertEqual(meta["sha256"], res.sha256)

    def test_download_checksum_mismatch_raises(self):
        d = DatasetDownloader()
        payload = b"abc"

        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "file.bin"
            def _get(url, headers=None, stream=True, timeout=1):
                return _FakeResponse(status_code=200, body=payload, headers={})

            fake_requests = types.SimpleNamespace(get=_get)
            with mock.patch.dict(sys.modules, {"requests": fake_requests}):
                with self.assertRaises(Exception):
                    d.download(url="https://example.com/x", output_path=out, expected_sha256="0" * 64, max_retries=1)

    def test_download_resume_appends_when_206(self):
        d = DatasetDownloader()
        initial = b"hello-"
        rest = b"world"
        payload = rest

        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "file.bin"
            out.write_bytes(initial)

            called = {}

            def _fake_get(url, headers=None, stream=True, timeout=1):
                called["range"] = headers.get("Range")
                return _FakeResponse(status_code=206, body=payload, headers={"Content-Length": str(len(payload))})

            fake_requests = types.SimpleNamespace(get=_fake_get)
            with mock.patch.dict(sys.modules, {"requests": fake_requests}):
                d.download(url="https://example.com/x", output_path=out, max_retries=1, resume=True)

            self.assertEqual(called["range"], f"bytes={len(initial)}-")
            self.assertEqual(out.read_bytes(), initial + rest)

    def test_extract_archive_zip_and_tar(self):
        d = DatasetDownloader()
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            # ZIP
            import zipfile

            zip_path = td_path / "a.zip"
            with zipfile.ZipFile(zip_path, "w") as zf:
                zf.writestr("x.txt", "hi")
            out_dir = td_path / "zip_out"
            d.extract_archive(archive_path=zip_path, extract_dir=out_dir)
            self.assertTrue((out_dir / "x.txt").exists())

            # TAR.GZ
            import tarfile

            tar_path = td_path / "b.tar.gz"
            src_file = td_path / "y.txt"
            src_file.write_text("yo", encoding="utf-8")
            with tarfile.open(tar_path, "w:gz") as tf:
                tf.add(src_file, arcname="y.txt")
            tar_out = td_path / "tar_out"
            d.extract_archive(archive_path=tar_path, extract_dir=tar_out)
            self.assertTrue((tar_out / "y.txt").exists())

    def test_ftp_scheme_dispatches_to_ftp_downloader(self):
        d = DatasetDownloader()
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "x.bin"
            with mock.patch.object(d, "_download_ftp") as m:
                m.return_value = None  # type: ignore
                d.download(url="ftp://example.com/file.bin", output_path=out, max_retries=1)
                self.assertTrue(m.called)

    def test_extract_archive_recursive_nested_zip(self):
        d = DatasetDownloader()
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            import zipfile

            inner = td_path / "inner.zip"
            with zipfile.ZipFile(inner, "w") as zf:
                zf.writestr("a.txt", "a")

            outer = td_path / "outer.zip"
            with zipfile.ZipFile(outer, "w") as zf:
                zf.write(inner, arcname="inner.zip")

            out_dir = td_path / "out"
            d.extract_archive_recursive(archive_path=outer, extract_dir=out_dir, max_depth=2)
            # outer extracted
            self.assertTrue((out_dir / "inner.zip").exists())
            # nested extracted
            self.assertTrue((out_dir / "inner_extracted" / "a.txt").exists())


if __name__ == "__main__":
    unittest.main()
