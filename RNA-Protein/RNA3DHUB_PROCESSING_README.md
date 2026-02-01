# RNA 3D Hub Dataset Processing

## Overview

This document describes the fixes implemented to generate trRosettaRNA2-compatible A3M files from the RNA 3D Hub dataset (13,782 RNAs).

## Critical Fixes Implemented

### 1. BLAST + MAFFT Branch

**Problem:** 
- MAFFT produced correct aligned FASTA with gaps
- `fasta_to_a3m()` function was removing gap columns
- Query sequence was not included in MAFFT input
- Result: Sequences had different lengths, no alignment information

**Solution:**
- Query sequence is now prepended to homologs file before MAFFT
- `fasta_to_a3m()` now preserves alignment columns
- Query is written first with NO gaps
- All other sequences have SAME length as query (gaps preserved as `-`)
- Each sequence written on ONE line

**Code changes:**
```python
# Before MAFFT, add query to homologs
with open(homologs_file, 'w') as f:
    f.write(f">{rna['id']}\n{rna['seq']}\n")  # Query first
    f.write(homolog_content)  # Then homologs

# After MAFFT, preserve alignment
fasta_to_a3m(aligned_file, a3m_file, query_id=rna['id'])
# Now keeps gaps and alignment structure
```

### 2. Infernal + Rfam Branch

**Problem:**
- `stockholm_to_a3m()` was removing alignment columns from Stockholm file
- Sequences had different lengths
- No gap characters preserved

**Solution:**
- Stockholm alignment columns now preserved during conversion
- Query sequence provided separately (not in seed alignment)
- All sequences aligned to query positions
- Gaps preserved as `-`

**Code changes:**
```python
stockholm_to_a3m(sto_file, a3m_file, 
                query_id=rna['id'], 
                query_sequence=rna['seq'])
# Now preserves alignment structure
```

### 3. Ambiguity Code Normalization

**Problem:**
- IUPAC ambiguity codes (W, R, Y, K, M, S, B, D, H, V) in sequences
- trRosettaRNA2 may not support these

**Solution:**
- All ambiguity codes normalized to 'N'
- Applied in both `fasta_to_a3m()` and `stockholm_to_a3m()`

**Code:**
```python
def normalize_sequence(seq):
    """Replace IUPAC ambiguity codes with N"""
    ambiguity_codes = 'WRYKMSBDHVN'
    normalized = []
    for char in seq.upper():
        if char in ambiguity_codes:
            normalized.append('N')
        elif char in 'ACGU-':
            normalized.append(char)
        elif char == '.':
            normalized.append('-')  # Stockholm gaps
        else:
            normalized.append(char)
    return ''.join(normalized)
```

### 4. Single-Line Sequences

**Problem:**
- Sequences were wrapped at 80 characters
- trRosettaRNA2 expects one line per sequence

**Solution:**
- Each sequence now written on ONE continuous line
- No line wrapping

### 5. Query Sequence First

**Problem:**
- Query sequence position was not guaranteed

**Solution:**
- Query sequence always written first in A3M file
- Query has NO gaps (ungapped)
- All other sequences aligned to query

## Expected A3M Format

Based on trRosettaRNA2 example (`trRosettaRNA2/example/seq.a3m`):

```
>query
CCGCCGCGCCAUGCCUGUGGCGGAAACCGCCGCGCCAUGCCUGUGGCGG
>homolog1
------------GCCCGTGGCGGAAACCGCCGCGCC-------------
>homolog2
CCGCCGCGCCGTGCCTGTG-----------CGCGCCGTGCCTGTGGCGC
>homolog3
CCGCCGCGCTGTGGCTGTGGCGGAA-CCGCCGCGGCATGCTCACG-CGG
```

**Key properties:**
- Query (first sequence): NO gaps, length N
- All other sequences: length N (same as query)
- Gaps represented as `-`
- Each sequence on ONE line
- No ambiguity codes (normalized to N)

## Validation

The `validate_a3m()` function now checks:
1. Query sequence (first) has no gaps
2. All sequences have same length
3. No ambiguity codes present
4. Proper format

## Usage

### Basic Usage

```bash
python RNA-Protein/process_rna3dhub_dataset.py \
  --input rna3dhub_sequences.fasta \
  --output output_a3m_files/ \
  --rfam-cm data/rfam/Rfam.cm \
  --rfam-seed data/rfam/Rfam.seed \
  --rnacentral-db data/rnacentral/rnacentral_active
```

### With Custom Tool Paths

```bash
python RNA-Protein/process_rna3dhub_dataset.py \
  --input rna3dhub_sequences.fasta \
  --output output_a3m_files/ \
  --rfam-cm data/rfam/Rfam.cm \
  --rfam-seed data/rfam/Rfam.seed \
  --rnacentral-db data/rnacentral/rnacentral_active \
  --cmscan /path/to/cmscan \
  --blastn /path/to/blastn \
  --blastdbcmd /path/to/blastdbcmd \
  --mafft /path/to/mafft
```

### Parallel Processing

For large datasets, process in chunks:

```bash
# Process first 5000 RNAs
python RNA-Protein/process_rna3dhub_dataset.py \
  --input rna3dhub_sequences.fasta \
  --output output_a3m_files/ \
  --rfam-cm data/rfam/Rfam.cm \
  --rfam-seed data/rfam/Rfam.seed \
  --rnacentral-db data/rnacentral/rnacentral_active \
  --start 0 \
  --end 5000

# Process next 5000 RNAs
python RNA-Protein/process_rna3dhub_dataset.py \
  --input rna3dhub_sequences.fasta \
  --output output_a3m_files/ \
  --rfam-cm data/rfam/Rfam.cm \
  --rfam-seed data/rfam/Rfam.seed \
  --rnacentral-db data/rnacentral/rnacentral_active \
  --start 5000 \
  --end 10000
```

## Output

### A3M Files

For each RNA, an A3M file is generated:
- `{RNA_ID}_blast.a3m` - if BLAST+MAFFT succeeded
- `{RNA_ID}_{RFAM_ACC}.a3m` - if Infernal+Rfam succeeded

### Results TSV

Tab-separated file with columns:
- `rna_id`: RNA identifier
- `length`: RNA sequence length
- `infernal_success`: Boolean, Infernal succeeded
- `blast_success`: Boolean, BLAST succeeded
- `hybrid_success`: Boolean, either method succeeded
- `final_method`: 'Infernal', 'BLAST', or None
- `a3m_file`: Name of generated A3M file
- `n_sequences`: Number of sequences in MSA
- `time`: Processing time in seconds

Example:
```
rna_id	length	infernal_success	blast_success	hybrid_success	final_method	a3m_file	n_sequences	time
RNA001	150	True	False	True	Infernal	RNA001_RF00001.a3m	45	12.3
RNA002	200	False	True	True	BLAST	RNA002_blast.a3m	78	25.6
RNA003	100	False	False	False	None	None	0	8.2
```

## Hybrid Strategy

The pipeline uses a two-step hybrid approach:

1. **Step 1: Infernal + Rfam**
   - Fast, structure-aware
   - Works for RNAs with known families
   - Uses curated seed alignments

2. **Step 2: BLAST + MAFFT** (fallback)
   - Slower but more general
   - Works for RNAs without Rfam families
   - Finds homologs by sequence similarity

## Expected Success Rates

Based on previous tests:
- Infernal-only: ~63%
- BLAST-only: ~74%
- **Hybrid: ~85-90%** (expected)

## Troubleshooting

### Issue: "Query sequence contains gaps"

**Cause:** Query was not properly identified or extracted
**Solution:** Check that query_id matches a sequence in the alignment

### Issue: "Sequences have different lengths"

**Cause:** Alignment columns were removed during conversion
**Solution:** Ensure using updated `fasta_to_a3m()` or `stockholm_to_a3m()`

### Issue: "Found ambiguity codes"

**Cause:** Normalization not applied
**Solution:** Ensure using updated conversion functions with normalization

### Issue: MAFFT timeout

**Cause:** Too many homologs or very long sequences
**Solution:** Increase timeout or reduce max_target_seqs in BLAST

## Files Modified

1. `RNA-Protein/src/rna_msa/alignment_formats.py`
   - `fasta_to_a3m()`: Preserves alignment, normalizes ambiguity codes
   - `stockholm_to_a3m()`: Preserves alignment, normalizes ambiguity codes
   - `validate_a3m()`: Enhanced validation checks

2. `RNA-Protein/test_hybrid_diverse.py`
   - Added query to MAFFT input
   - Uses updated conversion functions

3. `RNA-Protein/process_rna3dhub_dataset.py` (NEW)
   - Production script for RNA 3D Hub dataset
   - Supports parallelization
   - Comprehensive error handling

## Next Steps

After generating A3M files:

1. **Validate A3M files:**
   ```python
   from rna_msa.alignment_formats import validate_a3m
   stats = validate_a3m('output.a3m')
   print(stats)
   ```

2. **Generate trRosettaRNA2 features:**
   ```bash
   python trRosettaRNA2/trRNA2/predict.py \
     --input output.a3m \
     --output features.npz
   ```

3. **Extract DSSR features:**
   ```bash
   x3dna-dssr -i structure.pdb -o=dssr_output.json
   ```

4. **Combine features for training:**
   - 2D features from trRosettaRNA2 (distance/orientation maps)
   - 3D features from DSSR (structural annotations)
   - Sequence embeddings from RNA-FM

## Contact

For questions or issues, contact the research team.

## References

- trRosettaRNA2: https://github.com/uw-ipd/trRosettaRNA2
- Rfam: https://rfam.org/
- RNAcentral: https://rnacentral.org/
- RNA 3D Hub: https://rna.bgsu.edu/rna3dhub/
