import csv
from urllib.parse import quote
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
    for split in ("train", "val", "test"):
        with (args.input / f"{split}.csv").open(newline="") as source, (args.output / f"{split}.fasta").open("w") as fasta:
            for index, row in enumerate(csv.DictReader(source), 1):
                for role, key in (("mirna", "a"), ("target", "b")):
                    identifier = quote(row[f"{key}_name"], safe="")
                    fasta.write(f">{split}_{index:05d}/{role} label={row['label']} source_id={identifier}\n{row[f'{key}_seq']}\n")


if __name__ == "__main__":
    main()
