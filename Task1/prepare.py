import gzip
import json
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
    for split, count in (("train", 95), ("val", 24), ("test", 44)):
        rows = [json.loads(line) for line in (args.input / f"{split}.jsonl").read_text().splitlines()]
        if len(rows) != count:
            raise ValueError(f"Unexpected {split} size")
        with (args.output / f"{split}.fasta").open("w") as fasta, (args.output / f"{split}_rsa.tsv.gz").open("wb") as compressed:
            with gzip.GzipFile(fileobj=compressed, mode="wb", mtime=0) as labels:
                labels.write(b"id\tposition\trsa\n")
                for row in rows:
                    sequence = row["sequence"]
                    if not len(sequence) == len(row["target"]) == len(row["mask"]):
                        raise ValueError(row["id"])
                    fasta.write(f">{row['id']}\n{sequence}\n")
                    for position, (value, valid) in enumerate(zip(row["target"], row["mask"]), 1):
                        if valid:
                            labels.write(f"{row['id']}\t{position}\t{value}\n".encode())


if __name__ == "__main__":
    main()
