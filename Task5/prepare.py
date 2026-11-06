import csv
from urllib.parse import quote
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
    for split in ("train", "val", "test"):
        with (args.input / f"{split}.csv").open(newline="") as source, (args.output / f"{split}.fasta").open("w") as fasta:
            for index, row in enumerate(csv.DictReader(source), 1):
                for role, key in (("mirna", "a"), ("target", "b")):
                    identifier = quote(row[f"{key}_name"], safe="")
                    fasta.write(f">{split}_{index:05d}/{role} label={row['label']} source_id={identifier}\n{row[f'{key}_seq']}\n")


if __name__ == "__main__":
    main()
