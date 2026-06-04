from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Dict, Optional


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def sha256sum(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def write_json(path: Path, data: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True)


def read_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def sidecar_metadata_path(file_path: Path) -> Path:
    return file_path.with_suffix(file_path.suffix + ".metadata.json")


def load_sidecar_metadata(file_path: Path) -> Optional[Dict[str, Any]]:
    meta_path = sidecar_metadata_path(file_path)
    if not meta_path.exists():
        return None
    return read_json(meta_path)

