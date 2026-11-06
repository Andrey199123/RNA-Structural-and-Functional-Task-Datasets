import json

import numpy as np
import argparse
import atexit
import tempfile
import zipfile
from pathlib import Path


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path, nargs="?", default=Path(__file__).resolve().parent / "inputs.zip")
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    if args.input.is_file():
        temporary = tempfile.TemporaryDirectory()
        atexit.register(temporary.cleanup)
        destination = Path(temporary.name)
        with zipfile.ZipFile(args.input) as archive:
            for member in archive.infolist():
                path = Path(member.filename)
                if path.is_absolute() or ".." in path.parts or member.external_attr >> 16 & 0o170000 == 0o120000:
                    raise ValueError("Unsafe input archive member")
            archive.extractall(destination)
        args.input = destination
    args.output.mkdir(parents=True, exist_ok=True)
    return args


def main():
    args = arguments()
    for split, count in (("train", 168), ("val", 43), ("test", 80)):
        rows = [json.loads(line) for line in (args.input / f"{split}.jsonl").read_text().splitlines()]
        if len(rows) != count:
            raise ValueError(f"Unexpected {split} size")
        matrices = {}
        with (args.output / f"{split}.fasta").open("w") as fasta:
            for index, row in enumerate(rows):
                length = len(row["sequence"])
                target = np.asarray(row["target"], dtype=np.uint8)
                mask = np.asarray(row["pair_mask"], dtype=bool)
                if target.shape != (length, length) or mask.shape != target.shape:
                    raise ValueError(row["id"])
                fasta.write(f">{row['id']} matrix={index}\n{row['sequence']}\n")
                matrices[f"contact_{index}"] = target
                matrices[f"valid_{index}"] = mask
        np.savez_compressed(args.output / f"{split}_contacts.npz", **matrices)


if __name__ == "__main__":
    main()
