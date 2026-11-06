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
    for split, count in (("TR0", 10814), ("VL0", 1300), ("TS0", 1305)):
        rows = [json.loads(line) for line in (args.input / f"{split}.jsonl").read_text().splitlines()]
        if len(rows) != count:
            raise ValueError(f"Unexpected {split} size")
        with (args.output / f"{split}.fasta").open("w") as fasta, (args.output / f"{split}.dbn").open("w") as structures:
            for row in rows:
                sequence = row["sequence"]
                structure = row["secondary_structure"]
                if len(sequence) != len(structure):
                    raise ValueError(row["id"])
                fasta.write(f">{row['id']}\n{sequence}\n")
                structures.write(f">{row['id']}\n{sequence}\n{structure}\n")


if __name__ == "__main__":
    main()
