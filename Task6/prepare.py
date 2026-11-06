import csv
import math
import argparse
from pathlib import Path


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    return args


def main():
    args = arguments()
    for split in ("train", "Random7600", "Human7600"):
        with (args.input / f"{split}.csv").open(newline="") as source:
            rows = list(csv.DictReader(source))
        if not rows or (split != "train" and len(rows) != 7600):
            raise ValueError(f"Unexpected number of records in {split}")
        with (args.output / f"{split}.fasta").open("w") as fasta:
            for index, row in enumerate(rows, 1):
                sequence = row["sequence"].strip().upper()
                value = float(row["mrl"])
                if not sequence or set(sequence) - set("ACGTUN") or not math.isfinite(value):
                    raise ValueError(f"Invalid record {index} in {split}")
                fasta.write(f">{split}_{index:05d} mrl={value}\n{sequence}\n")


if __name__ == "__main__":
    main()
