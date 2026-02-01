# RINmaker - Residue Interaction Network Maker

## Overview
RINmaker is a fast, versatile tool to determine residue interaction networks (RINs) in proteins.

## Installation
Clone from: https://github.com/RINmaker/RINmaker

## Features
- Hydrogen bonds
- Van der Waals interactions
- Ionic bridges
- π-π stacking
- π-cation interactions
- Hydrophobic interactions
- Interaction energies (not just distances)

## Benchmark Tests (from paper - Table 2)
Performance comparison with RING 3.0 and PDBe-Arpeggio on proteins of increasing size:

| PDB ID | Residues | Description |
|--------|----------|-------------|
| 1al1   | 12       | Very small protein |
| 1l2y   | 20       | Trp-cage miniprotein |
| 1bk8   | 50       | Small protein |
| 1eod   | 100      | Medium protein |
| 1hhq   | 150      | Medium protein |
| 1ab5   | 250      | Large protein |
| 1bge   | 318      | Large protein |
| 1b1y   | 500      | Very large protein |
| 2vji   | 598      | Very large protein |
| 1miq   | 750      | Very large protein |

## Expected Performance
- RINmaker should be ~10-100x faster than RING 3.0
- RINmaker should be ~1000x faster than PDBe-Arpeggio
- Linear time complexity with number of residues
