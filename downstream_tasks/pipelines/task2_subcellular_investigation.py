from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np

from .base import PipelineContext, TaskPipeline
from ..types import TaskState, TaskStatus
from ..utils.files import write_json


class Task2SubcellularInvestigationPipeline(TaskPipeline):
    task_id = "task2_subcellular_investigation"
    task_name = "Task 2: Subcellular Localization (Investigation)"

    def __init__(
        self,
        ctx: PipelineContext,
        *,
        output_dir: Optional[Path] = None,
        config: Optional[Dict[str, Any]] = None,
    ):
        super().__init__(ctx, output_dir=output_dir, config=config)
        self.task_dir = ctx.eval_root / "Task2"
        # Preferred manual dataset location is under the evaluation directory:
        #   downstream-tasks-evaluation/data/raw/task2_subcellular/RNALocate_v2.0.txt
        #
        # Some users may place the file under repo_root/data/raw/... by mistake;
        # we check that as a secondary fallback.
        self.raw_dir = ctx.eval_root / "data" / "raw" / "task2_subcellular"
        self.rnalocate_v2_file = self.raw_dir / "RNALocate_v2.0.txt"
        self.repo_raw_dir = ctx.repo_root / "data" / "raw" / "task2_subcellular"
        self.repo_rnalocate_v2_file = self.repo_raw_dir / "RNALocate_v2.0.txt"
        self.legacy_data_file = self.task_dir / "data" / "All RNA subcellular localization information.txt"

        # Prefer the manually downloaded RNALocate v2.0 file placed under data/raw/.
        # Fall back to the legacy Task2/data export (annotations-only; typically no sequences).
        override = self.config.get("task2_data_file")
        if override:
            p = Path(str(override))
            self.data_file = p if p.is_absolute() else (ctx.eval_root / p)
        elif self.rnalocate_v2_file.exists():
            self.data_file = self.rnalocate_v2_file
        elif self.repo_rnalocate_v2_file.exists():
            self.data_file = self.repo_rnalocate_v2_file
        else:
            self.data_file = self.legacy_data_file
        self.dm3loc_dir = self.task_dir / "data" / "dm3loc"
        # DM3Loc provides FASTA files with labels embedded in the header.
        self.dm3loc_nonredundant_fasta = self.dm3loc_dir / "modified_multilabel_seq_nonredundent.fasta"
        self.dm3loc_full_fasta = self.dm3loc_dir / "modified_multilabel_seq.fasta"

    @staticmethod
    def _iter_dm3loc_fasta(path: Path):
        """
        Yield (label_bits, header, sequence_length) from DM3Loc-style FASTA.

        Header format begins with a 7-bit one-hot string (0/1) followed by comma-separated metadata:
          >0100000,ACCNUM:NM_...,Gene_ID:...,Gene_Name:...

        Label bit order (from DM3Loc code):
          0 Nucleus
          1 Exosome
          2 Cytosol
          3 Cytoplasm (auxiliary; often skipped in evaluation)
          4 Ribosome
          5 Membrane
          6 Endoplasmic reticulum (ER)
        """
        label = None
        header = None
        seq_len = 0
        with path.open("r", encoding="utf-8", errors="replace") as f:
            for raw in f:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None and label is not None:
                        yield label, header, seq_len
                    header = line[1:]
                    first = header.split(",", 1)[0]
                    if len(first) != 7 or any(c not in "01" for c in first):
                        raise ValueError(f"Unexpected DM3Loc label bits: {first!r}")
                    label = [int(c) for c in first]
                    seq_len = 0
                else:
                    # DM3Loc uses DNA alphabet (T); treat it as a sequence length only for investigation.
                    seq_len += len(line)
            if header is not None and label is not None:
                yield label, header, seq_len

    @staticmethod
    def _summarize_lengths(lengths: np.ndarray) -> Dict[str, float]:
        if lengths.size == 0:
            return {"n": 0.0, "min": 0.0, "p50": 0.0, "p90": 0.0, "p95": 0.0, "p99": 0.0, "max": 0.0, "mean": 0.0}
        return {
            "n": float(lengths.size),
            "min": float(lengths.min()),
            "p50": float(np.percentile(lengths, 50)),
            "p90": float(np.percentile(lengths, 90)),
            "p95": float(np.percentile(lengths, 95)),
            "p99": float(np.percentile(lengths, 99)),
            "max": float(lengths.max()),
            "mean": float(lengths.mean()),
        }

    def run(self, status: TaskStatus) -> TaskStatus:
        report_dir = self.ctx.eval_root / "reports" / self.task_id
        report_dir.mkdir(parents=True, exist_ok=True)

        artifacts = []
        notes = []

        # Preferred: investigate the paper-aligned mRNA benchmark used by DM3Loc/LOCAS.
        if self.dm3loc_nonredundant_fasta.exists():
            classes = ["Nucleus", "Exosome", "Cytosol", "Cytoplasm", "Ribosome", "Membrane", "ER"]

            lengths = []
            labels = []
            max_rec = {"length": 0, "header": None}
            for lb, header, seq_len in self._iter_dm3loc_fasta(self.dm3loc_nonredundant_fasta):
                lengths.append(seq_len)
                labels.append(lb)
                if seq_len > max_rec["length"]:
                    max_rec = {"length": seq_len, "header": header}

            arr = np.asarray(lengths, dtype=int)
            lab = np.asarray(labels, dtype=int)

            # DM3Loc notes truncating to 4000 nt at each end for some training/evaluation comparisons.
            trunc_end = int(self.config.get("truncate_ends_nt", 4000))
            trunc_arr = np.minimum(arr, 2 * trunc_end)

            per_class_counts = {classes[i]: int(lab[:, i].sum()) for i in range(len(classes))}
            # Multi-label statistics excluding the cytoplasm auxiliary label at index 3.
            non_cyto = lab.copy()
            non_cyto[:, 3] = 0
            n_labels = non_cyto.sum(axis=1)

            thresholds = [512, 1024, 2048, 4096, 8192, 12000, 20000]
            coverage = {f"<= {t}": float((arr <= t).mean()) for t in thresholds}
            trunc_coverage = {f"<= {t}": float((trunc_arr <= t).mean()) for t in thresholds}

            report = {
                "dataset": {
                    "name": "DM3Loc nonredundant mRNA benchmark (FASTA)",
                    "file": str(self.dm3loc_nonredundant_fasta),
                    "n_sequences": int(arr.size),
                    "label_encoding": {
                        "format": "7-bit one-hot prefix in FASTA header",
                        "class_order": classes,
                        "note": "Cytoplasm (index 3) is treated as auxiliary in DM3Loc code and often skipped in evaluation.",
                    },
                    "max_length_record": {"length": int(max_rec["length"]), "header": max_rec["header"]},
                },
                "lengths": {
                    "original": self._summarize_lengths(arr),
                    "truncated_ends": {
                        "truncate_ends_nt": trunc_end,
                        "max_after_truncation": 2 * trunc_end,
                        "summary": self._summarize_lengths(trunc_arr),
                    },
                    "coverage_by_max_len": {"original": coverage, "truncated_ends": trunc_coverage},
                },
                "labels": {
                    "per_class_positive_counts": per_class_counts,
                    "num_labels_excluding_cytoplasm": {
                        "min": int(n_labels.min()) if n_labels.size else 0,
                        "max": int(n_labels.max()) if n_labels.size else 0,
                        "mean": float(n_labels.mean()) if n_labels.size else 0.0,
                        "single_label_count": int((n_labels == 1).sum()) if n_labels.size else 0,
                        "multi_label_count": int((n_labels > 1).sum()) if n_labels.size else 0,
                        "zero_label_count": int((n_labels == 0).sum()) if n_labels.size else 0,
                    },
                },
                "notes_for_azam": [
                    "This DM3Loc-provided FASTA includes raw sequences and labels; it is suitable for length feasibility analysis.",
                    "Many sequences exceed common model context lengths (only a small fraction are <=1024 or <=2048).",
                    "DM3Loc mentions sequence lengths spanning ~200 to ~30,000 nt in the paper, but the provided FASTA includes some very long mRNAs (e.g., TTN ~104k nt).",
                    "Next decision: confirm whether we should (a) skip Task 2, (b) adopt DM3Loc's truncation (e.g., 4000 nt from each end), or (c) use a different benchmark with shorter sequences.",
                ],
                "sources": {
                    "dm3loc_paper_pdf": "https://pdfs.semanticscholar.org/fa93/530d6855d02ac18abba70df2cb4ebeeaa01d.pdf",
                    "dm3loc_github_testdata": "https://github.com/duolinwang/DM3Loc/tree/master/testdata",
                },
            }

            report_path = report_dir / "dm3loc_investigation_report.json"
            write_json(report_path, report)

            artifacts.append(str(report_path))
            notes.append(f"DM3Loc n={int(arr.size)}, max_len={int(arr.max())}")

        if not self.data_file.exists():
            if artifacts:
                status.output_dir = str(report_dir)
                status.artifacts = artifacts
                status.message = "; ".join(notes + [f"RNALocate file missing: {self.data_file}"])
                status.state = TaskState.SUCCESS
                return status

            status.state = TaskState.SKIPPED
            status.message = f"Missing dataset file: {self.data_file}"
            return status

        lengths = []
        total = 0
        # The RNALocate export is large; stream line-by-line and extract the sequence column heuristically.
        with self.data_file.open("r", encoding="utf-8", errors="replace") as f:
            header = f.readline()
            sep = "\t" if "\t" in header else ","
            cols = [c.strip().lower() for c in header.strip().split(sep)]
            seq_idx = None
            for i, c in enumerate(cols):
                if c in {"sequence", "seq", "rna_sequence"}:
                    seq_idx = i
                    break
            # Fallback: many exports keep the RNA sequence in a column named "RNA sequence" or similar.
            if seq_idx is None:
                for i, c in enumerate(cols):
                    if "sequence" in c:
                        seq_idx = i
                        break

            if seq_idx is None:
                # Many RNALocate dumps contain localization annotations but no sequences.
                # Record what we can and mark investigation as successful (but blocked for length stats).
                for _ in f:
                    total += 1

                report_path = report_dir / "schema_report.json"
                write_json(
                    report_path,
                    {
                        "file": str(self.data_file),
                        "delimiter": "tab" if sep == "\t" else "comma",
                        "columns": cols,
                        "n_rows_total_excluding_header": int(total),
                        "has_sequence_column": False,
                        "note": "No sequence column detected in this RNALocate export; sequence-length investigation requires a dataset/version that includes sequences.",
                    },
                )

                artifacts.append(str(report_path))
                notes.append("RNALocate: no sequence column (schema_report.json)")

                status.output_dir = str(report_dir)
                status.artifacts = artifacts
                status.message = "; ".join(notes)
                status.state = TaskState.SUCCESS
                return status

            for line in f:
                total += 1
                parts = line.rstrip("\n").split(sep)
                if seq_idx is None or seq_idx >= len(parts):
                    continue
                seq = parts[seq_idx].strip()
                if not seq:
                    continue
                lengths.append(len(seq))

        if not lengths:
            notes.append("RNALocate: sequence column present but no sequences extracted")
            status.output_dir = str(report_dir)
            status.artifacts = artifacts
            status.message = "; ".join(notes)
            status.state = TaskState.SUCCESS
            return status

        arr = np.asarray(lengths, dtype=int)
        stats: Dict[str, float] = {
            "n_rows_total": float(total),
            "n_sequences_parsed": float(arr.size),
            "min_len": float(arr.min()),
            "p50_len": float(np.percentile(arr, 50)),
            "p90_len": float(np.percentile(arr, 90)),
            "p95_len": float(np.percentile(arr, 95)),
            "p99_len": float(np.percentile(arr, 99)),
            "max_len": float(arr.max()),
            "mean_len": float(arr.mean()),
        }

        report_path = report_dir / "length_report.json"
        write_json(report_path, {"file": str(self.data_file), **stats})

        artifacts.append(str(report_path))
        notes.append(f"RNALocate: parsed {int(stats['n_sequences_parsed'])} sequences; max_len={int(stats['max_len'])}")

        status.output_dir = str(report_dir)
        status.artifacts = artifacts
        status.message = "; ".join(notes)
        status.state = TaskState.SUCCESS
        return status
