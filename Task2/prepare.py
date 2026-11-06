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
    for split, count in (("train", 12983), ("test", 1164)):
        records = []
        for label in ("positive", "negative"):
            path = args.input / split / label / f"{label}_{split}_data.fa"
            for record in path.read_text().split(">")[1:]:
                identifier, *lines = record.splitlines()
                sequence = "".join(lines).strip().upper()
                if len(sequence) != 174:
                    raise ValueError(identifier)
                records.append((identifier.split()[0], sequence, label))
        if len(records) != count:
            raise ValueError(f"Unexpected {split} size: {len(records)}")
        with (args.output / f"{split}.fasta").open("w") as handle:
            for index, (identifier, sequence, label) in enumerate(records, 1):
                handle.write(f">{split}_{index:05d} label={int(label == 'positive')} source_id={identifier}\n{sequence}\n")


if __name__ == "__main__":
    main()
