from __future__ import annotations

import json
import os
import subprocess
import time
from urllib.parse import urlparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

from .errors import DownloadError
from .utils.files import sha256sum, sidecar_metadata_path, write_json


@dataclass(frozen=True)
class DownloadResult:
    path: Path
    url: str
    sha256: str
    bytes_downloaded: int
    resumed: bool


class DatasetDownloader:
    def __init__(self, *, user_agent: str = "downstream_tasks/0.1"):
        self.user_agent = user_agent

    def download(
        self,
        *,
        url: str,
        output_path: Path,
        expected_sha256: Optional[str] = None,
        timeout_s: int = 60,
        max_retries: int = 3,
        resume: bool = True,
    ) -> DownloadResult:
        parsed = urlparse(url)
        if parsed.scheme.lower() == "ftp":
            return self._download_ftp(
                url=url,
                output_path=output_path,
                expected_sha256=expected_sha256,
                timeout_s=timeout_s,
                max_retries=max_retries,
                resume=resume,
            )

        try:
            import requests  # type: ignore
        except Exception:
            return self._download_http_urllib(
                url=url,
                output_path=output_path,
                expected_sha256=expected_sha256,
                timeout_s=timeout_s,
                max_retries=max_retries,
                resume=resume,
            )

        output_path.parent.mkdir(parents=True, exist_ok=True)

        existing = output_path.exists()
        existing_size = output_path.stat().st_size if existing else 0
        resumed = False

        headers = {"User-Agent": self.user_agent}
        if resume and existing_size > 0:
            headers["Range"] = f"bytes={existing_size}-"

        last_err: Optional[BaseException] = None
        for attempt in range(1, max_retries + 1):
            try:
                with requests.get(url, headers=headers, stream=True, timeout=timeout_s) as r:
                    if r.status_code >= 400:
                        raise DownloadError(f"HTTP {r.status_code} for {url}")

                    # If server ignored Range header, restart the download.
                    if "Range" in headers and r.status_code == 200:
                        existing_size = 0
                        resumed = False
                        try:
                            output_path.unlink()
                        except FileNotFoundError:
                            pass

                    mode = "ab" if ("Range" in headers and existing_size > 0 and r.status_code == 206) else "wb"
                    resumed = mode == "ab"

                    total = r.headers.get("Content-Length")
                    total_int = int(total) if total and total.isdigit() else None
                    if resumed and total_int is not None:
                        total_int += existing_size

                    try:
                        from tqdm import tqdm
                    except Exception:  # pragma: no cover
                        tqdm = None  # type: ignore

                    progress = None
                    if tqdm is not None:
                        progress = tqdm(
                            total=total_int,
                            initial=existing_size if resumed else 0,
                            unit="B",
                            unit_scale=True,
                            desc=output_path.name,
                        )

                    bytes_written = existing_size if resumed else 0
                    with output_path.open(mode) as f:
                        for chunk in r.iter_content(chunk_size=1024 * 1024):
                            if not chunk:
                                continue
                            f.write(chunk)
                            bytes_written += len(chunk)
                            if progress is not None:
                                progress.update(len(chunk))
                    if progress is not None:
                        progress.close()

                actual_sha = sha256sum(output_path)
                if expected_sha256 is not None and actual_sha.lower() != expected_sha256.lower():
                    raise DownloadError(
                        f"Checksum mismatch for {output_path} (expected {expected_sha256}, got {actual_sha})"
                    )

                meta = {
                    "url": url,
                    "downloaded_at_epoch": int(time.time()),
                    "sha256": actual_sha,
                    "bytes": int(output_path.stat().st_size),
                    "resumed": bool(resumed),
                }
                write_json(sidecar_metadata_path(output_path), meta)
                return DownloadResult(
                    path=output_path,
                    url=url,
                    sha256=actual_sha,
                    bytes_downloaded=int(output_path.stat().st_size),
                    resumed=resumed,
                )
            except BaseException as e:
                last_err = e
                if attempt == max_retries:
                    break
                time.sleep(2 ** (attempt - 1))

        raise DownloadError(f"Failed to download {url}") from last_err

    def _download_http_urllib(
        self,
        *,
        url: str,
        output_path: Path,
        expected_sha256: Optional[str],
        timeout_s: int,
        max_retries: int,
        resume: bool,
    ) -> DownloadResult:
        import urllib.request

        output_path.parent.mkdir(parents=True, exist_ok=True)
        existing_size = output_path.stat().st_size if (resume and output_path.exists()) else 0

        headers = {"User-Agent": self.user_agent}
        if resume and existing_size > 0:
            headers["Range"] = f"bytes={existing_size}-"

        last_err: Optional[BaseException] = None
        for attempt in range(1, max_retries + 1):
            try:
                req = urllib.request.Request(url, headers=headers, method="GET")
                with urllib.request.urlopen(req, timeout=timeout_s) as resp:
                    status = getattr(resp, "status", 200)

                    # If server ignored Range header, restart.
                    if "Range" in headers and status != 206:
                        existing_size = 0
                        try:
                            output_path.unlink()
                        except FileNotFoundError:
                            pass

                    mode = "ab" if ("Range" in headers and existing_size > 0 and status == 206) else "wb"
                    resumed = mode == "ab"
                    bytes_written = existing_size if resumed else 0
                    with output_path.open(mode) as f:
                        while True:
                            chunk = resp.read(1024 * 1024)
                            if not chunk:
                                break
                            f.write(chunk)
                            bytes_written += len(chunk)

                actual_sha = sha256sum(output_path)
                if expected_sha256 is not None and actual_sha.lower() != expected_sha256.lower():
                    raise DownloadError(
                        f"Checksum mismatch for {output_path} (expected {expected_sha256}, got {actual_sha})"
                    )

                meta = {
                    "url": url,
                    "downloaded_at_epoch": int(time.time()),
                    "sha256": actual_sha,
                    "bytes": int(output_path.stat().st_size),
                    "resumed": bool(resumed),
                }
                write_json(sidecar_metadata_path(output_path), meta)
                return DownloadResult(
                    path=output_path,
                    url=url,
                    sha256=actual_sha,
                    bytes_downloaded=int(output_path.stat().st_size),
                    resumed=bool(resumed),
                )
            except BaseException as e:
                last_err = e
                if attempt == max_retries:
                    break
                time.sleep(2 ** (attempt - 1))

        raise DownloadError(f"Failed to download {url}") from last_err

    def _download_ftp(
        self,
        *,
        url: str,
        output_path: Path,
        expected_sha256: Optional[str],
        timeout_s: int,
        max_retries: int,
        resume: bool,
    ) -> DownloadResult:
        import ftplib

        parsed = urlparse(url)
        host = parsed.hostname
        if not host:
            raise DownloadError(f"Invalid FTP URL (missing host): {url}")
        port = parsed.port or 21
        path = parsed.path or "/"

        output_path.parent.mkdir(parents=True, exist_ok=True)
        existing_size = output_path.stat().st_size if (resume and output_path.exists()) else 0

        last_err: Optional[BaseException] = None
        for attempt in range(1, max_retries + 1):
            try:
                with ftplib.FTP() as ftp:
                    ftp.connect(host=host, port=port, timeout=timeout_s)
                    ftp.login()  # anonymous by default

                    resumed = False
                    mode = "wb"
                    rest = None
                    if existing_size > 0:
                        try:
                            ftp.sendcmd("TYPE I")  # binary
                            rest = existing_size
                            mode = "ab"
                            resumed = True
                        except Exception:
                            rest = None
                            mode = "wb"
                            resumed = False

                    bytes_written = existing_size if resumed else 0
                    with output_path.open(mode) as f:
                        def _cb(chunk: bytes) -> None:
                            nonlocal bytes_written
                            f.write(chunk)
                            bytes_written += len(chunk)

                        if rest is not None:
                            ftp.retrbinary(f"RETR {path}", _cb, rest=rest)
                        else:
                            ftp.retrbinary(f"RETR {path}", _cb)

                actual_sha = sha256sum(output_path)
                if expected_sha256 is not None and actual_sha.lower() != expected_sha256.lower():
                    raise DownloadError(
                        f"Checksum mismatch for {output_path} (expected {expected_sha256}, got {actual_sha})"
                    )

                meta = {
                    "url": url,
                    "downloaded_at_epoch": int(time.time()),
                    "sha256": actual_sha,
                    "bytes": int(output_path.stat().st_size),
                    "resumed": bool(resumed),
                }
                write_json(sidecar_metadata_path(output_path), meta)
                return DownloadResult(
                    path=output_path,
                    url=url,
                    sha256=actual_sha,
                    bytes_downloaded=int(output_path.stat().st_size),
                    resumed=bool(resumed),
                )
            except BaseException as e:
                last_err = e
                if attempt == max_retries:
                    break
                time.sleep(2 ** (attempt - 1))

        raise DownloadError(f"Failed to download (FTP) {url}") from last_err

    def clone_git_repo(
        self,
        *,
        repo_url: str,
        output_dir: Path,
        depth: Optional[int] = 1,
        branch: Optional[str] = None,
    ) -> Path:
        if output_dir.exists():
            return output_dir
        output_dir.parent.mkdir(parents=True, exist_ok=True)

        cmd: List[str] = ["git", "clone"]
        if depth is not None:
            cmd += ["--depth", str(depth)]
        if branch is not None:
            cmd += ["--branch", branch]
        cmd += [repo_url, str(output_dir)]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise DownloadError(f"git clone failed: {result.stderr.strip()}")
        return output_dir

    def extract_archive(self, *, archive_path: Path, extract_dir: Path) -> List[Path]:
        import tarfile
        import zipfile

        extract_dir.mkdir(parents=True, exist_ok=True)
        out: List[Path] = []

        suffixes = "".join(archive_path.suffixes)
        if suffixes.endswith(".tar.gz") or suffixes.endswith(".tgz") or suffixes.endswith(".tar.bz2"):
            with tarfile.open(archive_path, "r:*") as tf:
                # Python 3.12+ supports `filter=` to reduce tar extraction risks.
                try:
                    tf.extractall(extract_dir, filter="data")  # type: ignore[call-arg]
                except TypeError:
                    tf.extractall(extract_dir)
                out = [extract_dir / m.name for m in tf.getmembers() if m.name]
            return out

        if suffixes.endswith(".zip"):
            with zipfile.ZipFile(archive_path, "r") as zf:
                zf.extractall(extract_dir)
                out = [extract_dir / n for n in zf.namelist() if n]
            return out

        raise DownloadError(f"Unsupported archive format: {archive_path}")

    def extract_archive_recursive(
        self,
        *,
        archive_path: Path,
        extract_dir: Path,
        max_depth: int = 2,
    ) -> List[Path]:
        """
        Extract an archive and, if extracted files contain supported nested archives,
        extract them as well (up to max_depth).
        """
        extracted: List[Path] = self.extract_archive(archive_path=archive_path, extract_dir=extract_dir)
        if max_depth <= 0:
            return extracted

        nested: List[Path] = []
        for p in extracted:
            if not p.is_file():
                continue
            suffixes = "".join(p.suffixes)
            if suffixes.endswith(".tar.gz") or suffixes.endswith(".tgz") or suffixes.endswith(".tar.bz2") or suffixes.endswith(".zip"):
                nested.append(p)

        all_paths = list(extracted)
        for n in nested:
            child_dir = extract_dir / (n.stem + "_extracted")
            child_dir.mkdir(parents=True, exist_ok=True)
            all_paths.extend(self.extract_archive_recursive(archive_path=n, extract_dir=child_dir, max_depth=max_depth - 1))
        return all_paths
