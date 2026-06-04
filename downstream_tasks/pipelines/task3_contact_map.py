from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np

from .base import PipelineContext, TaskPipeline
from ..downloader import DatasetDownloader
from ..types import TaskState, TaskStatus
from ..utils.files import write_json
from ..validator import DataValidator


def _read_fasta(path: Path) -> str:
    seq_parts: List[str] = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            seq_parts.append(line)
    return "".join(seq_parts).strip().upper().replace("T", "U")


def _pick_rna_chain(structure) -> object:
    # Pick the chain with the most nucleic residues.
    from Bio.PDB.Polypeptide import is_nucleic

    model = next(structure.get_models())
    best_chain = None
    best_count = -1
    for chain in model.get_chains():
        count = 0
        for res in chain.get_residues():
            if is_nucleic(res, standard=False):
                count += 1
        if count > best_count:
            best_count = count
            best_chain = chain
    if best_chain is None:
        raise RuntimeError("No chains found in PDB")
    return best_chain


def _residue_atoms(chain) -> Tuple[List[object], np.ndarray, np.ndarray]:
    from Bio.PDB.Polypeptide import is_nucleic

    residues: List[object] = []
    coords: List[List[float]] = []
    res_idx: List[int] = []
    for res in chain.get_residues():
        if not is_nucleic(res, standard=False):
            continue
        residues.append(res)
        i = len(residues) - 1
        for atom in res.get_atoms():
            # Skip hydrogens if present.
            if atom.element == "H":
                continue
            c = atom.get_coord()
            coords.append([float(c[0]), float(c[1]), float(c[2])])
            res_idx.append(i)
    if not residues:
        raise RuntimeError("No nucleic residues found in selected chain")
    if not coords:
        raise RuntimeError("No atoms found in selected chain residues")
    return residues, np.asarray(coords, dtype=np.float32), np.asarray(res_idx, dtype=np.int32)


def _residues_to_sequence(residues: List[object]) -> str:
    seq: List[str] = []
    for r in residues:
        name = (r.get_resname() or "").strip().upper()
        if not name:
            seq.append("N")
            continue
        # Common cases: "A", "C", "G", "U" or deoxy "DA", "DG", ...
        if name in {"A", "C", "G", "U"}:
            seq.append(name)
            continue
        if name in {"DA", "DC", "DG", "DT"}:
            seq.append("U" if name == "DT" else name[-1])
            continue
        base = name[-1]
        seq.append("U" if base == "T" else (base if base in {"A", "C", "G", "U"} else "N"))
    return "".join(seq)


def compute_contact_map(
    *,
    pdb_path: Path,
    distance_threshold_angstrom: float = 8.0,
    long_range_min_sep: int = 24,
) -> Tuple[str, np.ndarray]:
    from Bio.PDB import PDBParser
    from scipy.spatial import cKDTree

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    chain = _pick_rna_chain(structure)
    residues, coords, res_idx = _residue_atoms(chain)
    seq = _residues_to_sequence(residues)

    l = len(residues)
    contact = np.zeros((l, l), dtype=np.uint8)

    tree = cKDTree(coords)
    for a, b in tree.query_pairs(r=distance_threshold_angstrom):
        i = int(res_idx[a])
        j = int(res_idx[b])
        if i == j:
            continue
        contact[i, j] = 1
        contact[j, i] = 1

    if long_range_min_sep > 0:
        ii, jj = np.indices(contact.shape)
        contact[np.abs(ii - jj) < long_range_min_sep] = 0

    np.fill_diagonal(contact, 0)
    return seq, contact


class Task3ContactMapPipeline(TaskPipeline):
    task_id = "task3_contact_map"
    task_name = "Task 3: Contact Map Prediction (RNAcontact)"

    TRAIN_URL = "https://yanglab.qd.sdu.edu.cn/RNAcontact/benchmark/train_set.tar.bz2"
    TEST_URL = "https://yanglab.qd.sdu.edu.cn/RNAcontact/benchmark/test_set.tar.bz2"

    def __init__(
        self,
        ctx: PipelineContext,
        *,
        output_dir: Optional[Path] = None,
        config: Optional[Dict[str, Any]] = None,
    ):
        super().__init__(ctx, output_dir=output_dir, config=config)
        self.task_dir = ctx.eval_root / "Task3"
        self.data_dir = self.task_dir / "data"
        self.prep_dir = self.task_dir / "preprocessed"

    def _ensure_dataset(self) -> Tuple[Path, Path]:
        d = DatasetDownloader()
        self.data_dir.mkdir(parents=True, exist_ok=True)

        train_archive = self.data_dir / "train_set.tar.bz2"
        test_archive = self.data_dir / "test_set.tar.bz2"
        if not train_archive.exists():
            self.logger.info("Downloading %s", self.TRAIN_URL)
            d.download(url=self.TRAIN_URL, output_path=train_archive, timeout_s=180)
        if not test_archive.exists():
            self.logger.info("Downloading %s", self.TEST_URL)
            d.download(url=self.TEST_URL, output_path=test_archive, timeout_s=180)

        train_root = self.data_dir / "train_set"
        test_root = self.data_dir / "test_set"
        if not train_root.exists():
            self.logger.info("Extracting %s", train_archive)
            d.extract_archive(archive_path=train_archive, extract_dir=train_root)
        if not test_root.exists():
            self.logger.info("Extracting %s", test_archive)
            d.extract_archive(archive_path=test_archive, extract_dir=test_root)

        train_set = train_root / "train_set"
        test_set = test_root / "test_set"
        return train_set, test_set

    def _write_split(
        self,
        *,
        split_name: str,
        ids: List[str],
        fasta_dir: Path,
        pdb_dir: Path,
        long_range_min_sep: int = 24,
    ) -> List[str]:
        out_dir = self.prep_dir / split_name
        out_dir.mkdir(parents=True, exist_ok=True)
        artifacts: List[str] = []

        for sid in ids:
            out_path = out_dir / f"{sid}.npz"
            if out_path.exists():
                artifacts.append(str(out_path))
                continue

            fasta_path = fasta_dir / f"{sid}.fasta"
            pdb_path = pdb_dir / f"{sid}.pdb"
            if not fasta_path.exists() or not pdb_path.exists():
                self.logger.warning("Missing inputs for %s (fasta=%s pdb=%s)", sid, fasta_path.exists(), pdb_path.exists())
                continue

            seq = _read_fasta(fasta_path)
            pdb_seq, contact = compute_contact_map(pdb_path=pdb_path, long_range_min_sep=long_range_min_sep)
            if contact.shape[0] != len(pdb_seq):
                raise RuntimeError(f"Internal error for {sid}: pdb_seq_len != contact_map_len")
            if len(seq) != len(pdb_seq):
                self.logger.warning(
                    "Length mismatch for %s (fasta_len=%d pdb_len=%d) — using PDB-derived sequence",
                    sid,
                    len(seq),
                    len(pdb_seq),
                )

            np.savez_compressed(
                out_path,
                id=sid,
                sequence=np.asarray(pdb_seq.encode("ascii")),
                fasta_sequence=np.asarray(seq.encode("ascii")),
                contact_map=contact,
            )
            artifacts.append(str(out_path))
        return artifacts

    def run(self, status: TaskStatus) -> TaskStatus:
        train_set, test_set = self._ensure_dataset()

        train_fasta = train_set / "train_fasta"
        train_pdb = train_set / "train_pdb"
        test_fasta = test_set / "test_fasta"
        test_pdb = test_set / "test_pdb"

        train_ids = sorted(p.stem for p in train_fasta.glob("*.fasta"))
        test_ids = sorted(p.stem for p in test_fasta.glob("*.fasta"))

        if not train_ids or not test_ids:
            status.state = TaskState.FAILED
            status.message = "No FASTA files found after extraction"
            return status

        seed = int(self.config.get("seed", 42))
        train_ratio = float(self.config.get("train_ratio", 0.8))
        long_range_min_sep = int(self.config.get("long_range_min_sep", 24))

        # Reproduce the ERNIE-RNA convention: split training set ~8:2 into TR168/VL43.
        rng = np.random.default_rng(seed)
        perm = rng.permutation(len(train_ids))
        # Use floor to match the canonical 211 → 168/43 split at train_ratio=0.8.
        n_train = int(np.floor(train_ratio * len(train_ids)))
        train_sel = [train_ids[i] for i in perm[:n_train]]
        val_sel = [train_ids[i] for i in perm[n_train:]]

        self.prep_dir.mkdir(parents=True, exist_ok=True)
        artifacts: List[str] = []
        artifacts += self._write_split(
            split_name="train",
            ids=train_sel,
            fasta_dir=train_fasta,
            pdb_dir=train_pdb,
            long_range_min_sep=long_range_min_sep,
        )
        artifacts += self._write_split(
            split_name="val",
            ids=val_sel,
            fasta_dir=train_fasta,
            pdb_dir=train_pdb,
            long_range_min_sep=long_range_min_sep,
        )
        artifacts += self._write_split(
            split_name="test",
            ids=test_ids,
            fasta_dir=test_fasta,
            pdb_dir=test_pdb,
            long_range_min_sep=long_range_min_sep,
        )

        summary = {
            "source": "RNAcontact benchmark (yanglab.qd.sdu.edu.cn)",
            "distance_threshold_angstrom": 8.0,
            "long_range_min_sep": long_range_min_sep,
            "n_total_train_set": len(train_ids),
            "n_train": len(train_sel),
            "n_val": len(val_sel),
            "n_test": len(test_ids),
            "seed": seed,
            "train_ratio": train_ratio,
        }
        write_json(self.prep_dir / "summary.json", summary)
        artifacts.append(str(self.prep_dir / "summary.json"))

        # Validation report
        v = DataValidator()
        checks = [
            v.validate_counts(
                observed=len(train_ids),
                expected=211,
                tolerance_ratio=0.0,
                name="train_set_total(TR221)",
            ),
            v.validate_counts(
                observed=len(test_ids),
                expected=80,
                tolerance_ratio=0.0,
                name="test_set_total(TS80)",
            ),
            v.validate_no_leakage(train_ids=train_sel, val_ids=val_sel, test_ids=test_ids),
        ]
        report = v.report(checks)
        report_path = self.prep_dir / "validation.json"
        v.write_report(report=report, path=report_path)
        artifacts.append(str(report_path))

        status.output_dir = str(self.prep_dir)
        status.artifacts = artifacts
        status.message = f"Prepared contact maps (train={len(train_sel)} val={len(val_sel)} test={len(test_ids)})"
        return status
