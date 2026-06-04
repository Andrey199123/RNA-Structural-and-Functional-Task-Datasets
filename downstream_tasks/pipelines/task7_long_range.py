from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np

from .base import PipelineContext, TaskPipeline
from ..types import TaskState, TaskStatus
from ..utils.files import write_json
from ..validator import DataValidator


BRACKET_PAIRS = {"(": ")", "[": "]", "{": "}", "<": ">"}
BRACKET_CLOSE_TO_OPEN = {v: k for k, v in BRACKET_PAIRS.items()}


def parse_dot_bracket(struct: str) -> List[Tuple[int, int]]:
    stacks: Dict[str, List[int]] = {op: [] for op in BRACKET_PAIRS}
    pairs: List[Tuple[int, int]] = []
    for i, ch in enumerate(struct):
        if ch in BRACKET_PAIRS:
            stacks[ch].append(i)
            continue
        if ch in BRACKET_CLOSE_TO_OPEN:
            op = BRACKET_CLOSE_TO_OPEN[ch]
            if not stacks[op]:
                continue
            j = stacks[op].pop()
            pairs.append((j, i))
    return pairs


def has_pseudoknot(struct: str) -> bool:
    return any(ch in "[]{}<>" for ch in struct)


class Task7LongRangePipeline(TaskPipeline):
    task_id = "task7_long_range"
    task_name = "Task 7: Long-Range Base Pairs (bpRNA-1m subset)"

    def __init__(
        self,
        ctx: PipelineContext,
        *,
        output_dir: Optional[Path] = None,
        config: Optional[Dict[str, Any]] = None,
    ):
        super().__init__(ctx, output_dir=output_dir, config=config)
        self.task_dir = ctx.eval_root / "Task7"
        self.prep_dir = self.task_dir / "preprocessed"

    def run(self, status: TaskStatus) -> TaskStatus:
        # Use the official TR0 split from `multimolecule/bprna-spot-0` if available locally.
        data_dir = self.task_dir / "data" / "bprna_spot0"
        tr0_jsonl = data_dir / "TR0.jsonl"
        if not tr0_jsonl.exists():
            status.state = TaskState.SKIPPED
            status.message = f"Missing local TR0 JSONL: {tr0_jsonl} (run Task7/download_bprna_tr0.py)"
            return status

        max_len = int(self.config.get("max_length", 1024))
        seed = int(self.config.get("seed", 42))
        long_range_min_sep = int(self.config.get("long_range_min_sep", 101))
        import json

        ids: List[str] = []
        seqs: List[str] = []
        lengths: List[int] = []
        offsets: List[int] = [0]
        pairs_all: List[Tuple[int, int]] = []

        n_pk = 0
        with tr0_jsonl.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                row = json.loads(line)
                seq = str(row.get("sequence", "")).upper().replace("T", "U")
                st = str(row.get("secondary_structure", ""))
                if not seq or not st:
                    continue
                if len(seq) > max_len:
                    continue
                if has_pseudoknot(st):
                    n_pk += 1

                pairs = [(i, j) for (i, j) in parse_dot_bracket(st) if abs(i - j) >= long_range_min_sep]
                pairs_all.extend(pairs)
                offsets.append(offsets[-1] + len(pairs))

                ids.append(str(row.get("id", "")))
                seqs.append(seq)
                lengths.append(len(seq))

        self.prep_dir.mkdir(parents=True, exist_ok=True)
        out_tsv = self.prep_dir / "TR0.tsv"
        with out_tsv.open("w", encoding="utf-8") as f:
            for i, s, line in zip(ids, seqs, lengths):
                # Note: structure is not written here to keep file smaller; use TR0.jsonl for full records.
                f.write(f"{i}\t{s}\t{line}\n")

        pairs_arr = np.asarray(pairs_all, dtype=np.int32) if pairs_all else np.zeros((0, 2), dtype=np.int32)
        offsets_arr = np.asarray(offsets, dtype=np.int64)

        out_npz = self.prep_dir / "TR0_long_range_pairs.npz"
        np.savez_compressed(
            out_npz,
            ids=np.asarray(ids, dtype=object),
            sequences=np.asarray(seqs, dtype=object),
            lengths=np.asarray(lengths, dtype=np.int32),
            pair_offsets=offsets_arr,
            pairs=pairs_arr,
            min_sep=np.asarray(long_range_min_sep, dtype=np.int32),
        )

        long_range_pairs_total = int(pairs_arr.shape[0])

        summary = {
            "source": str(tr0_jsonl),
            "split": "TR0 (bpRNA-spot-0)",
            "filters": {"max_length": max_len},
            "n_sequences_used": int(len(ids)),
            "seed": seed,
            "pseudoknots": {"n_with_pseudoknot": int(n_pk)},
            "long_range": {"min_sep": long_range_min_sep, "total_pairs_in_subset": int(long_range_pairs_total)},
        }
        write_json(self.prep_dir / "summary.json", summary)

        v = DataValidator()
        checks = [
            v.validate_counts(observed=len(ids), expected=10814, tolerance_ratio=0.0, name="TR0_count"),
            v.validate_counts(observed=n_pk, expected=999, tolerance_ratio=0.0, name="TR0_pseudoknot_count"),
        ]
        report_path = self.prep_dir / "validation.json"
        v.write_report(report=v.report(checks), path=report_path)

        status.output_dir = str(self.prep_dir)
        status.artifacts = [str(out_tsv), str(out_npz), str(self.prep_dir / "summary.json"), str(report_path)]
        status.message = f"Prepared TR0 long-range pairs (n={len(ids)}; pseudoknots={n_pk})"
        return status
