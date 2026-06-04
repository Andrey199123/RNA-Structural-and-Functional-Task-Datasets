from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from ..errors import ExternalToolError
from .subprocess_utils import require_executable, run


@dataclass(frozen=True)
class DSSRResult:
    json_path: Path
    n_pairs: int


def run_dssr(
    *,
    pdb_path: Path,
    json_out: Path,
    exe: str = "x3dna-dssr",
    timeout_s: int = 120,
) -> DSSRResult:
    exe_path = require_executable(exe)
    json_out.parent.mkdir(parents=True, exist_ok=True)

    argv = [exe_path, f"-i={pdb_path}", "--json", f"-o={json_out}"]
    res = run(argv, timeout_s=timeout_s)
    if res.returncode != 0:
        raise ExternalToolError(f"DSSR failed: {res.stderr.strip() or res.stdout.strip()}")
    if not json_out.exists():
        raise ExternalToolError(f"DSSR did not write JSON output: {json_out}")

    data = json.loads(json_out.read_text(encoding="utf-8"))
    pairs = data.get("pairs", [])
    n_pairs = len(pairs) if isinstance(pairs, list) else 0
    return DSSRResult(json_path=json_out, n_pairs=n_pairs)


def extract_pairs(dssr_json: Dict[str, Any]) -> List[Tuple[str, str]]:
    """
    Extract base pairs from DSSR JSON as list of (nt1_id, nt2_id).
    DSSR uses nt_id-like strings such as 'A.G1'.
    """
    pairs = dssr_json.get("pairs", [])
    out: List[Tuple[str, str]] = []
    if not isinstance(pairs, list):
        return out

    for p in pairs:
        if not isinstance(p, dict):
            continue
        nt1 = p.get("nt1") or p.get("nt1_id") or p.get("nt1_name") or p.get("nt1_long")
        nt2 = p.get("nt2") or p.get("nt2_id") or p.get("nt2_name") or p.get("nt2_long")
        if isinstance(nt1, str) and isinstance(nt2, str):
            out.append((nt1, nt2))
            continue
        # DSSR sometimes stores the combined string in 'nt1_nt2' or similar.
        combo = p.get("nt1_nt2") or p.get("nts_long")
        if isinstance(combo, str) and "," in combo:
            a, b = combo.split(",", 1)
            out.append((a.strip(), b.strip()))
    return out

