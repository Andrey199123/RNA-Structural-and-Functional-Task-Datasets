import h5py
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
    inputs = sorted(args.input.glob("*.h5"))
    if len(inputs) != 17:
        raise ValueError("Expected the 17 HeLa HDF5 exports")
    for split in ("train", "val", "test"):
        tracks = {}
        with (args.output / f"{split}.fasta").open("w") as fasta:
            for path in inputs:
                with h5py.File(path) as source:
                    if "72/8/20" not in source.attrs["split_method"]:
                        raise ValueError(f"Unexpected split for {path.name}")
                    group = source[split]
                    tracks[path.stem] = group["icshape"][:]
                    for index, (sequence, label) in enumerate(zip(group["sequences"], group["labels"])):
                        fasta.write(f">{path.stem}:{index} label={int(label)}\n{sequence.decode('ascii')}\n")
        np.savez_compressed(args.output / f"{split}_icshape.npz", **tracks)


if __name__ == "__main__":
    main()
