#!/usr/bin/env python3
"""
Download bpRNA-1m TR0 split for Task 7 (Long-Range Base Pair Prediction)

This script downloads the **official TR0/VL0/TS0 split** (CD-HIT-EST 80% + ~8:1:1)
as released by the SPOT-RNA team on Hugging Face:

  `multimolecule/bprna-spot-0`

Counts:
  TR0 = 10,814
  VL0 = 1,300
  TS0 = 1,305
"""

import os
import json
from pathlib import Path

def download_bprna_tr0():
    """
    Download the bpRNA TR0/VL0/TS0 split using Hugging Face datasets.
    """
    try:
        from datasets import load_dataset
    except ImportError:
        print("ERROR: datasets library not installed")
        print("Install with: pip install datasets")
        return False
    
    print("=" * 80)
    print("Downloading bpRNA TR0/VL0/TS0 for Task 7")
    print("=" * 80)
    print()
    print("Dataset: Hugging Face (multimolecule/bprna-spot-0)")
    print("Splits:  TR0 (train)=10814, VL0 (validation)=1300, TS0 (test)=1305")
    print()
    print("NOTE: These are the standard benchmark splits used by SPOT-RNA / UFold.")
    print()
    
    # Create output directory
    output_dir = Path(__file__).parent / "data" / "bprna_spot0"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Output directory: {output_dir}")
    print()
    print("Downloading dataset... (this may take a few minutes)")
    print()
    
    try:
        ds = load_dataset("multimolecule/bprna-spot-0")

        print("✓ Downloaded dataset")
        for split in ["train", "validation", "test"]:
            print(f"  - {split}: {len(ds[split])}")
        print()
        
        # Save dataset information
        info = {
            "dataset_name": "bpRNA-spot-0 (TR0/VL0/TS0)",
            "source": "Hugging Face (multimolecule/bprna-spot-0)",
            "splits": {
                "TR0": len(ds["train"]),
                "VL0": len(ds["validation"]),
                "TS0": len(ds["test"]),
            },
            "task": "Task 7: Long-Range Base Pair Prediction",
            "usage": "Zero-shot evaluation only (no fine-tuning)",
            "columns": list(ds["train"].column_names) if hasattr(ds["train"], "column_names") else []
        }
        
        info_file = output_dir / "dataset_info.json"
        with open(info_file, 'w') as f:
            json.dump(info, f, indent=2)
        
        print(f"✓ Saved dataset info to: {info_file}")
        print()
        
        # Save splits
        tr0_path = output_dir / "TR0.jsonl"
        vl0_path = output_dir / "VL0.jsonl"
        ts0_path = output_dir / "TS0.jsonl"

        ds["train"].to_json(tr0_path)
        ds["validation"].to_json(vl0_path)
        ds["test"].to_json(ts0_path)

        print(f"✓ Saved TR0 → {tr0_path}")
        print(f"✓ Saved VL0 → {vl0_path}")
        print(f"✓ Saved TS0 → {ts0_path}")
        print()
        
        # Print sample
        print("Sample sequence:")
        print(f"  ID: {ds['train'][0].get('id', 'N/A')}")
        print(f"  Sequence length: {len(ds['train'][0].get('sequence', ''))}")
        print(f"  Structure: {ds['train'][0].get('secondary_structure', 'N/A')[:50]}...")
        print()
        
        print("=" * 80)
        print("✓ Download complete!")
        print("=" * 80)
        print()
        print("Next steps:")
        print("1. Filter for sequences with pseudoknots (paper reports ~999 for TR0)")
        print("2. Extract long-range base pairs (|i-j| ≥ 101)")
        print("3. Prepare for zero-shot evaluation (no fine-tuning!)")
        print()
        
        return True
        
    except Exception as e:
        print(f"✗ Error downloading dataset: {e}")
        print()
        print("Alternative options:")
        print("1. Visit https://bprna.cgrb.oregonstate.edu/ and download TR0 split manually")
        print("2. Contact the bpRNA authors for the TR0 split")
        print("3. Check if the dataset is available on Zenodo or Figshare")
        return False

if __name__ == "__main__":
    success = download_bprna_tr0()
    exit(0 if success else 1)
