from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from ..errors import ExternalToolError
from .subprocess_utils import require_executable, run


def recommended_word_size(identity: float) -> int:
    """
    Recommended -n word size for cd-hit-est.
    Based on CD-HIT user guide heuristics.
    """
    if identity >= 0.90:
        return 8
    if identity >= 0.88:
        return 7
    if identity >= 0.85:
        return 6
    if identity >= 0.80:
        return 5
    if identity >= 0.75:
        return 4
    raise ValueError("cd-hit-est generally requires identity >= 0.75 for stable word sizes")


@dataclass(frozen=True)
class CDHitEstResult:
    representative_fasta: Path
    clusters_path: Path


def run_cdhit_est(
    *,
    input_fasta: Path,
    output_fasta: Path,
    identity: float = 0.80,
    threads: int = 0,
    memory_mb: int = 0,
    word_size: Optional[int] = None,
    exe: str = "cd-hit-est",
) -> CDHitEstResult:
    exe_path = require_executable(exe)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    n = word_size if word_size is not None else recommended_word_size(identity)
    argv = [
        exe_path,
        "-i",
        str(input_fasta),
        "-o",
        str(output_fasta),
        "-c",
        f"{identity:.3f}",
        "-n",
        str(n),
        "-T",
        str(int(threads)),
        "-M",
        str(int(memory_mb)),
        "-d",
        "0",
    ]
    res = run(argv)
    if res.returncode != 0:
        raise ExternalToolError(f"cd-hit-est failed: {res.stderr.strip() or res.stdout.strip()}")

    clusters = Path(str(output_fasta) + ".clstr")
    if not clusters.exists():
        raise ExternalToolError(f"cd-hit-est did not produce cluster file: {clusters}")
    return CDHitEstResult(representative_fasta=output_fasta, clusters_path=clusters)


def parse_clstr(path: Path) -> List[List[str]]:
    """
    Parse a .clstr file into clusters of sequence IDs (FASTA defline tokens).
    Returns list of clusters, each as list of IDs (representative first).
    """
    clusters: List[List[str]] = []
    cur: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">Cluster"):
                if cur:
                    clusters.append(cur)
                cur = []
                continue
            # Example lines:
            # 0   123nt, >seq_id... *
            # 1   120nt, >seq_id... at 99.2%
            if ">" in line:
                after = line.split(">", 1)[1]
                seq_id = after.split("...", 1)[0].strip()
                cur.append(seq_id)
    if cur:
        clusters.append(cur)
    return clusters

