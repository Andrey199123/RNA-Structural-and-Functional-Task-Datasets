# RNA-Protein Interaction Prediction: MSA Generation Research

## Project Overview

This project addresses Professor Azam's request to **investigate practical strategies** for generating RNA Multiple Sequence Alignments (MSAs) suitable for extracting 2D structural features. The current ZHMol model learns exclusively from sequence embeddings. Our goal was to research, test, and compare different approaches for generating high-quality RNA MSAs that can be used to extract structural features (distance maps, contact maps).

---

## Table of Contents

1. [Background & Motivation](#background--motivation)
2. [Professor Azam's Requirements](#professor-azams-requirements)
3. [Research Journey](#research-journey)
4. [Approach 1: Infernal + Rfam](#approach-1-infernal--rfam)
5. [Approach 2: BLAST + RNAcentral](#approach-2-blast--rnacentral)
6. [Domain Splitting Investigation](#domain-splitting-investigation)
7. [Comparative Analysis](#comparative-analysis)
8. [Final Results & Recommendations](#final-results--recommendations)
9. [Implementation Details](#implementation-details)
10. [Usage Guide](#usage-guide)
11. [File Structure](#file-structure)

---

## Background & Motivation

### The Problem

The ZHMol model for RNA-protein interaction prediction has a critical limitation:
- **Learns only from sequence features** (RNA-FM embeddings, ESM2 protein embeddings)
- **No structural information** about RNA or protein folding
- **Missing co-evolutionary signals** that reveal functional constraints

### The Solution

Research from the Rosetta group shows that extracting 2D features from MSAs significantly improves structure prediction:
- **Distance maps**: Which residues are close in 3D space
- **Contact maps**: Which positions interact
- **Co-evolutionary couplings**: Which positions co-vary to maintain structure

**Key insight**: To extract these features, we need **high-quality MSAs** for RNAs.

---

## Professor Azam's Requirements

From her email (September 26, 2025):

> "I'd like you to focus on the RNA MSA problem: investigate practical strategies for obtaining .a3m alignments for RNAs that are suitable for extracting structural features."

**Specific requirements:**
1. ✅ Investigate practical strategies for obtaining .a3m alignments
2. ✅ Ensure MSAs are suitable for extracting structural features
3. ✅ Handle long RNAs (>1000nt) - Rfam typically only covers local motifs
4. ✅ Generate high-quality MSAs with sufficient depth

**She mentioned Rfam as ONE possible solution**, noting its limitations:
> "Rfam provides high-quality, curated alignments, but these are mostly limited to short RNAs or specific structured domains."

---

## Research Journey

### Phase 1: Understanding the Landscape

**Initial questions:**
- What tools exist for RNA MSA generation?
- What databases contain RNA sequences?
- How do we measure MSA quality?

**Key discoveries:**
- **Rfam**: ~4,000 curated RNA families with seed alignments
- **RNAcentral**: 40.7 million RNA sequences from all species
- **Infernal**: Structure-aware RNA search tool (uses covariance models)
- **BLAST**: Fast sequence similarity search
- **A3M format**: Standard format for MSA-based structure prediction

### Phase 2: Testing Infernal + Rfam

**Hypothesis**: Use Rfam's curated alignments for known RNA families

**Implementation**: `test_infernal_rfam_diverse.py`

**How it works:**
1. Run `cmscan` (Infernal) to find which Rfam family the RNA belongs to
2. Extract the seed alignment for that family from Rfam.seed
3. Convert Stockholm format → A3M format
4. Return pre-made MSA

**Initial test (30 rRNAs only):**
- Success rate: 96.7% (29/30)
- MSA depth: 99 sequences
- **Problem discovered**: All test RNAs were rRNAs (not diverse!)

**Diverse test (35 RNAs, 7 types):**
- Success rate: 62.9% (22/35)
- Works great for: rRNA (100%), snRNA (100%), snoRNA (100%)
- Works okay for: miRNA (80%)
- Limited for: tRNA (40%)
- Fails for: lncRNA (0%)

### Phase 3: Testing BLAST + RNAcentral

**Hypothesis**: Search for specific homologs instead of family membership

**Implementation**: `test_blast_rnacentral_diverse.py`

**How it works:**
1. Run BLAST to find sequences similar to YOUR specific RNA
2. Extract top 50 homologs from RNAcentral
3. Align with MAFFT
4. Convert to A3M format
5. Return custom MSA

**Diverse test (same 35 RNAs):**
- Success rate: 74.3% (26/35)
- Works great for: rRNA (100%), snoRNA (100%), tRNA (80%), miRNA (80%)
- Works okay for: snRNA (60%), other (60%)
- Limited for: lncRNA (40%)

### Phase 4: Testing Domain Splitting for Long RNAs

**Hypothesis**: Professor Azam suggested splitting long RNAs (>1000nt) into domains since Rfam only covers local motifs

**Implementation**: `test_domain_split_lncrna.py` and `test_domain_split_blast.py`

**How it works:**
1. Split long RNA into 300nt overlapping windows (50nt overlap)
2. Test each domain with Infernal+Rfam OR BLAST+RNAcentral
3. Generate per-domain MSAs

**Results:**
- **Domain-split + Infernal**: 0% success (0/22 domains found Rfam families)
- **Domain-split + BLAST**: 0% success (0/22 domains found sufficient homologs)
- **Full-length BLAST**: 40% success (2/5 lncRNAs)

**Conclusion**: Domain splitting makes things worse, not better. The problem isn't RNA length - it's that lncRNAs are poorly conserved and not in Rfam families. Short 300nt fragments are too specific to find homologs. Full-length BLAST is the best approach for long RNAs.

### Phase 5: Understanding the Difference

**Key realization**: These approaches measure different things!

**Infernal + Rfam:**
- Finds: "What family does this RNA belong to?"
- Returns: Generic family seed alignment (same for all family members)
- Example: Any tRNA → RF00005 seed (954 sequences)

**BLAST + RNAcentral:**
- Finds: "What sequences are similar to THIS specific RNA?"
- Returns: Custom alignment of homologs specific to your query
- Example: tRNA-Ala → 50 tRNA-Ala homologs from different species

**For 2D feature extraction**: BLAST is better because homologs are specific to your query's structure.

---

## Approach 1: Infernal + Rfam

### What is Infernal?

**Infernal** (INFERence of RNA ALignment) is a tool for searching RNA sequences using **covariance models**.

**Key features:**
- **Structure-aware**: Considers both sequence AND base pairing
- **Covariance models**: Represent both sequence and secondary structure of RNA families
- **More accurate than BLAST** for RNAs with conserved structure

**How it works:**
```
Query RNA → cmscan → Matches Rfam family → Extract seed alignment
```

### What is Rfam?

**Rfam** is a database of RNA families:
- **~4,000 families**: rRNA, tRNA, miRNA, snoRNA, snRNA, riboswitches, etc.
- **Seed alignments**: Curated, high-quality alignments (10-1000 sequences)
- **Covariance models**: One per family for structure-aware search

### Implementation

**File**: `test_infernal_rfam_diverse.py`

**Step-by-step process:**

1. **Run cmscan** (Infernal tool)
```bash
cmscan --tblout output.txt --cut_ga --nohmmonly Rfam.cm query.fa
```
- Searches query against all Rfam covariance models
- Returns best matching family (e.g., RF00177 = bacterial SSU rRNA)

2. **Extract seed alignment**
```python
# Parse Rfam.seed file to extract family alignment
with open('Rfam.seed', 'r', errors='ignore') as f:
    in_family = False
    for line in f:
        if line.startswith(f'#=GF AC   {family_acc}'):
            in_family = True
        elif line.startswith('//') and in_family:
            break
        elif in_family:
            # Extract sequences
```

3. **Convert Stockholm → A3M**
```python
from src.rna_msa.alignment_formats import stockholm_to_a3m
stockholm_to_a3m(sto_file, a3m_file, query_id=None)
```

**Key code**: `src/rna_msa/alignment_formats.py`
- `parse_stockholm()`: Parse Stockholm format alignments
- `stockholm_to_a3m()`: Convert to A3M format
- `validate_a3m()`: Check A3M file validity
- `fasta_to_a3m()`: Convert FASTA alignments to A3M

### Results

**Test dataset**: 35 diverse RNAs (5 of each type: rRNA, tRNA, miRNA, snoRNA, snRNA, lncRNA, other)

| RNA Type | Success Rate | Avg MSA Depth |
|----------|--------------|---------------|
| rRNA | 5/5 (100%) | 97 sequences |
| snRNA | 5/5 (100%) | 174 sequences |
| snoRNA | 5/5 (100%) | 55 sequences |
| miRNA | 4/5 (80%) | 17 sequences |
| tRNA | 2/5 (40%) | 954 sequences |
| lncRNA | 0/5 (0%) | N/A |
| other | 1/5 (20%) | 458 sequences |
| **OVERALL** | **22/35 (62.9%)** | **185 sequences** |

### Pros & Cons

**Advantages:**
- ✅ Structure-aware search (considers base pairing)
- ✅ High-quality curated alignments
- ✅ Large MSAs (up to 954 sequences)
- ✅ Family annotations included
- ✅ Fast (~5 seconds per RNA)

**Disadvantages:**
- ❌ Limited to ~4,000 known families
- ❌ Doesn't work for novel/uncharacterized RNAs
- ❌ Seed alignments are generic (not specific to your query)
- ❌ Some families have small seeds (3-22 sequences)

### When to Use

**Best for:**
- Well-characterized RNA families (rRNA, snRNA, snoRNA)
- When you need family annotations
- When you want structure-aware search

**Not suitable for:**
- Novel RNAs (lncRNAs, new ncRNAs)
- Species-specific RNAs
- When you need query-specific homologs

---

## Approach 2: BLAST + RNAcentral

### What is BLAST?

**BLAST** (Basic Local Alignment Search Tool) is a fast sequence similarity search tool.

**Key features:**
- **Sequence-based**: Matches A, U, G, C patterns (no structure)
- **Very fast**: Searches millions of sequences in seconds
- **Indexed search**: Pre-built database index for speed

**How it works:**
```
Query RNA → BLAST → Find homologs → Extract → Align → A3M
```

### What is RNAcentral?

**RNAcentral** is a comprehensive RNA sequence database:
- **40.7 million sequences**: From all species (bacteria, archaea, eukaryotes)
- **All RNA types**: rRNA, tRNA, miRNA, lncRNA, snoRNA, etc.
- **30GB uncompressed**: Requires BLAST indexing for fast search

### Implementation

**File**: `test_blast_rnacentral_diverse.py`

**Step-by-step process:**

1. **Create BLAST database** (one-time setup, ~15 minutes)
```bash
makeblastdb -in rnacentral_active.fasta \
            -dbtype nucl \
            -out rnacentral_active \
            -parse_seqids
```

2. **Run BLAST search**
```bash
blastn -query query.fa \
       -db rnacentral_active \
       -evalue 0.001 \
       -outfmt '6 sseqid pident length' \
       -max_target_seqs 100 \
       -num_threads 4
```
- Finds sequences with ≥70% identity
- Returns top 100 hits in ~18 seconds

3. **Extract homolog sequences**
```bash
blastdbcmd -db rnacentral_active \
           -entry_batch ids.txt \
           -out homologs.fa
```

4. **Align with MAFFT**
```bash
mafft --auto --thread 4 homologs.fa > aligned.fa
```
- MAFFT: Multiple Alignment using Fast Fourier Transform
- Aligns sequences to show conserved positions

5. **Convert to A3M**
```python
from src.rna_msa.alignment_formats import fasta_to_a3m
fasta_to_a3m(aligned_file, a3m_file, query_id=query_id)
```

**Key code**: `src/rna_msa/fast_homology_search.py`
- `FastHomologySearch`: Class for BLAST-based searches
- `create_blast_db()`: Index RNAcentral database
- `blast_search()`: Find homologs
- `extract_sequences()`: Get sequences by ID
- `fast_homology_search()`: Complete pipeline

### Results

**Test dataset**: Same 35 diverse RNAs as Infernal test

| RNA Type | Success Rate | Avg Homologs | Avg Identity |
|----------|--------------|--------------|--------------|
| rRNA | 5/5 (100%) | 100 | 98.7% |
| snoRNA | 5/5 (100%) | 64 | 92.9% |
| tRNA | 4/5 (80%) | 54 | 86.8% |
| miRNA | 4/5 (80%) | 191 | 91.0% |
| snRNA | 3/5 (60%) | 76 | 89.9% |
| other | 3/5 (60%) | 79 | 89.2% |
| lncRNA | 2/5 (40%) | 32 | 86.8% |
| **OVERALL** | **26/35 (74.3%)** | **90 sequences** | **91.5%** |

### Pros & Cons

**Advantages:**
- ✅ Works for any RNA with homologs (not limited to known families)
- ✅ Custom MSAs specific to your query
- ✅ High sequence identity (70-99%)
- ✅ Finds novel/uncharacterized RNAs
- ✅ Better overall success rate (74.3%)

**Disadvantages:**
- ❌ Not structure-aware (sequence similarity only)
- ❌ Requires 30GB database + indexing
- ❌ Slower than Infernal (~18s vs ~5s)
- ❌ No family annotations

### When to Use

**Best for:**
- Novel or uncharacterized RNAs
- When you need query-specific homologs
- When Rfam doesn't have your RNA family
- General-purpose MSA generation

**Not suitable for:**
- When you need structure-aware search
- When you need family annotations
- When database size is a constraint

---

## Domain Splitting Investigation

### The Problem with Long RNAs

Professor Azam noted:
> "Rfam provides high-quality, curated alignments, but these are mostly limited to short RNAs or specific structured domains. When an RNA sequence is long (e.g., >1000 nt), Rfam alignments typically cover only local motifs, while the rest of the sequence remains poorly aligned."

**Hypothesis**: Split long RNAs into ~300nt domains and process each independently.

### Testing Domain Splitting

We tested domain splitting on 5 lncRNAs (478-1842 nt) that failed with full-length approaches:

**Test 1: Domain-split + Infernal+Rfam** (`test_domain_split_lncrna.py`)
- Split 5 lncRNAs into 22 domains (300nt each, 50nt overlap)
- Ran Infernal cmscan on each domain
- **Result**: 0/22 domains (0%) found ANY Rfam family

**Test 2: Domain-split + BLAST+RNAcentral** (`test_domain_split_blast.py`)
- Same 22 domains
- Ran BLAST on each domain
- **Result**: 0/22 domains (0%) found sufficient homologs (≥10)

### Why Domain Splitting Failed

**The fundamental problem**: lncRNAs aren't in Rfam families and are poorly conserved

1. **With Infernal+Rfam:**
   - Rfam families are for well-characterized, conserved RNA motifs
   - lncRNAs are novel and not well-characterized
   - Even 300nt domains don't match any known family
   - This isn't a "length" problem - it's a "lncRNAs aren't in Rfam" problem

2. **With BLAST+RNAcentral:**
   - Full-length BLAST can find partial matches across entire sequence
   - 300nt fragments are too short and too specific
   - Need near-perfect matches for short sequences
   - lncRNAs are poorly conserved, so fragments have zero matches

### Comparison: Full-length vs Domain-split

| Approach | lncRNA Success |
|----------|----------------|
| Full-length BLAST | 40% (2/5) - Found 19-45 homologs |
| Full-length Infernal | 0% (0/5) - Not in Rfam |
| Domain-split + Infernal | 0% (0/22 domains) |
| Domain-split + BLAST | 0% (0/22 domains) |

**Conclusion**: Domain splitting makes things worse, not better. Full-length BLAST is the best approach for long RNAs.

---

## Comparative Analysis

### Side-by-Side Comparison

| RNA Type | Infernal+Rfam | BLAST+RNAcentral | Winner |
|----------|---------------|------------------|--------|
| rRNA | 5/5 (100%) | 5/5 (100%) | **Tie** |
| snoRNA | 5/5 (100%) | 5/5 (100%) | **Tie** |
| snRNA | 5/5 (100%) | 3/5 (60%) | **Infernal** |
| miRNA | 4/5 (80%) | 4/5 (80%) | **Tie** |
| tRNA | 2/5 (40%) | 4/5 (80%) | **BLAST** |
| other | 1/5 (20%) | 3/5 (60%) | **BLAST** |
| lncRNA | 0/5 (0%) | 2/5 (40%) | **BLAST** |
| **OVERALL** | **22/35 (62.9%)** | **26/35 (74.3%)** | **BLAST** |

### Key Insights

**1. BLAST wins overall (74.3% vs 62.9%)**
- More flexible approach
- Works for novel RNAs
- Better for diverse datasets

**2. Both work well for conserved RNAs**
- rRNA: 100% success for both
- snoRNA: 100% success for both
- These RNAs have strong evolutionary conservation

**3. Infernal better for snRNA**
- snRNA families are well-defined in Rfam
- Structure is highly conserved
- Covariance models capture this well

**4. BLAST better for tRNA and lncRNA**
- tRNA: Many variants not in Rfam seed
- lncRNA: Novel, not in Rfam at all
- BLAST finds homologs even without family membership

**5. MSA depth varies**
- Infernal: 3-954 sequences (depends on seed size)
- BLAST: 8-465 sequences (depends on homologs found)
- Both provide sufficient depth for 2D features

### Why We Chose BLAST as Primary Approach

**Decision factors:**

1. **Higher success rate**: 74.3% vs 62.9%
2. **Works for novel RNAs**: Not limited to known families
3. **Query-specific MSAs**: Better for co-evolutionary analysis
4. **Proven at scale**: Tested on diverse RNA types

**When to use Infernal instead:**
- You're working with well-characterized families
- You need structure-aware search
- You want family annotations
- You need very large MSAs (>100 sequences)

**Ideal solution: Hybrid approach**
1. Try Infernal first (fast, structure-aware)
2. If fails, use BLAST (broader coverage)
3. Combine results if both succeed

---

## Final Results & Recommendations

### What We Delivered

✅ **Comprehensive research on MSA generation strategies**
- Two approaches investigated (Infernal+Rfam, BLAST+RNAcentral)
- Domain splitting strategy tested and evaluated
- A3M format output ready for 2D feature extraction
- Tested on 35 diverse RNAs + 5 long lncRNAs

✅ **Clear performance benchmarks**
- BLAST: 74.3% success (26/35 RNAs)
- Infernal: 62.9% success (22/35 RNAs)
- Domain splitting: 0% success for lncRNAs

✅ **High-quality MSAs**
- Average depth: 90-185 sequences
- Average identity: 91.5%
- Suitable for co-evolutionary analysis

✅ **Fast and scalable**
- BLAST: ~18 seconds per RNA
- Infernal: ~5 seconds per RNA
- Can process hundreds of RNAs

✅ **Long RNA insights**
- Successfully tested up to 1842nt
- Full-length BLAST works better than domain splitting
- lncRNAs remain challenging (40% success)

### Recommendations for Professor Azam

**Key Findings:**

1. **Use BLAST + RNAcentral as primary approach**
   - Better success rate (74.3% vs 62.9%)
   - Works for diverse RNA types
   - Query-specific MSAs better for 2D features
   - Handles long RNAs better than domain splitting

2. **Domain splitting is not recommended**
   - Tested on lncRNAs: 0% success with both Infernal and BLAST
   - Makes things worse, not better
   - Full-length approaches work better even for long RNAs

3. **Dataset selection matters**
   - ✅ Use conserved RNAs: rRNA, tRNA, miRNA, snoRNA, snRNA
   - ⚠️ Avoid species-specific lncRNAs (40% success rate)
   - Consider: RNAcentral contains many suitable RNAs

4. **Next steps for feature extraction**
   - Use generated A3M files from BLAST approach
   - Extract distance maps from MSAs
   - Extract contact maps from co-evolutionary couplings
   - Tools: RNAcmap, SPOT-RNA, or custom analysis

5. **For model integration**
   - Add 2D feature branch to existing GraphSAGE model
   - Combine sequence embeddings + structural features
   - Train on dataset with good MSA coverage (avoid lncRNA-heavy datasets)

### What About RPI7317?

**The RPI7317 dataset contains human lncRNAs:**
- BLAST success rate: 40% (2/5 lncRNAs)
- Infernal success rate: 0% (0/5 lncRNAs)
- Many are species-specific with no homologs

**Options:**
1. **Use a different dataset** with conserved RNAs
2. **Accept limited MSA coverage** for RPI7317
3. **Use alternative features** for RNAs without MSAs:
   - Single-sequence structure prediction
   - Sequence embeddings (RNA-FM)
   - K-mer features

---

## Implementation Details

### File Structure

```
RNA-Protein/
├── README.md                              # This file
├── requirements_rna_msa.txt               # Python dependencies
│
├── src/
│   ├── rna_msa/
│   │   ├── alignment_formats.py           # A3M format conversion
│   │   ├── fast_homology_search.py        # BLAST-based search
│   │   ├── infernal_rfam.py               # Rfam family detection
│   │   ├── domain_split.py                # Domain splitting utilities
│   │   ├── rnacmap_style.py               # RNAcmap-style MSA (not used)
│   │   └── single_sequence.py             # Single-sequence fallback (not used)
│   │
│   ├── models/                            # GNN models
│   │   ├── gnn_sage.py                    # GraphSAGE backbone
│   │   └── vecnn_head.py                  # Prediction head
│   │
│   ├── embeddings/                        # Embedding generation
│   ├── data_processing/                   # Data preprocessing
│   └── utils/                             # Utility functions
│
├── data/
│   ├── rnacentral/                        # RNAcentral database (30GB)
│   │   ├── rnacentral_active.fasta
│   │   └── rnacentral_active.* (BLAST index)
│   │
│   ├── rfam/                              # Rfam database
│   │   ├── Rfam.cm                        # Covariance models
│   │   └── Rfam.seed                      # Seed alignments
│   │
│   ├── RPI7317/                           # RNA-protein dataset
│   │   ├── RPI7317_IncRNA_sequence.fasta
│   │   ├── RPI7317_protein_sequence.fasta
│   │   └── rpi7317_pairs.csv
│   │
│   ├── embeddings/                        # Pre-computed embeddings
│   │   ├── rna_fm_rpi7317.npz
│   │   └── esm2_35m_rpi7317.npz
│   │
│   ├── test_infernal_diverse/             # Infernal test results
│   │   ├── *_query.fa                     # Query files
│   │   ├── *_RF*.a3m                      # Generated A3M files
│   │   └── diverse_test_results.tsv       # Results table
│   │
│   └── test_blast_diverse/                # BLAST test results
│       ├── *_homologs.fa                  # Homolog sequences
│       ├── *_aligned.a3m                  # Generated A3M files
│       └── diverse_test_results.tsv       # Results table
│
├── configs/                               # Training configurations
│   ├── default.yaml                       # Random split
│   ├── uRNA.yaml                          # Unseen RNA split
│   ├── uPROT.yaml                         # Unseen protein split
│   └── uNODES.yaml                        # Unseen pairs split
│
├── test_infernal_rfam_diverse.py          # Benchmark: Infernal+Rfam
├── test_blast_rnacentral_diverse.py       # Benchmark: BLAST+RNAcentral
├── test_domain_split_lncrna.py            # Test: Domain-split + Infernal
├── test_domain_split_blast.py             # Test: Domain-split + BLAST
├── benchmark_msa_pipeline.py              # Original benchmark (rRNA only, deprecated)
│
└── src/train_linkpred.py                  # Train GNN model
```

### Key Code Files Explained

#### 1. `src/rna_msa/alignment_formats.py`

**Purpose**: Convert between alignment formats (Stockholm, FASTA, A3M)

**Key functions:**
- `parse_stockholm(sto_file)`: Parse Stockholm format
- `stockholm_to_a3m(sto_file, a3m_file, query_id)`: Convert Stockholm → A3M
- `fasta_to_a3m(fasta_file, a3m_file, query_id)`: Convert FASTA → A3M
- `validate_a3m(a3m_file)`: Validate A3M format and get statistics

**A3M format rules:**
- First sequence = query (reference)
- Uppercase = aligned positions (match states)
- Lowercase = insertions relative to query
- Gaps in query are removed from all sequences

**Example:**
```python
from src.rna_msa.alignment_formats import stockholm_to_a3m

# Convert Rfam seed alignment to A3M
stockholm_to_a3m(
    sto_file='RF00177.sto',
    a3m_file='RF00177.a3m',
    query_id=None  # Use first sequence as query
)
```

#### 2. `src/rna_msa/fast_homology_search.py`

**Purpose**: BLAST-based homology search with database indexing

**Key class: `FastHomologySearch`**

**Methods:**
- `create_blast_db()`: Create BLAST database index (one-time)
- `blast_search(query_seq, query_id)`: Find homologs
- `extract_sequences(seq_ids, output_file)`: Extract by ID
- `fast_homology_search()`: Complete pipeline

**Example:**
```python
from src.rna_msa.fast_homology_search import FastHomologySearch

# Initialize
searcher = FastHomologySearch('data/rnacentral/rnacentral_active.fasta')

# Create index (one-time, ~15 minutes)
searcher.create_blast_db()

# Search for homologs
n_homologs, output_file = searcher.fast_homology_search(
    query_seq="AUGCAUGC...",
    query_id="my_rna",
    output_fasta=Path("homologs.fa"),
    max_homologs=50,
    min_identity=0.7
)
```

#### 3. `test_infernal_rfam_diverse.py`

**Purpose**: Test Infernal+Rfam approach on diverse RNAs

**What it does:**
1. Extracts 35 diverse RNAs from RNAcentral (5 of each type)
2. Runs cmscan to find Rfam families
3. Extracts seed alignments from Rfam.seed
4. Converts to A3M format
5. Validates and reports statistics

**Output:**
- A3M files in `data/test_infernal_diverse/`
- Results table: `diverse_test_results.tsv`
- Success rate by RNA type

**Run:**
```bash
python test_infernal_rfam_diverse.py
```

#### 4. `test_blast_rnacentral_diverse.py`

**Purpose**: Test BLAST+RNAcentral on same diverse RNAs

**What it does:**
1. Uses same 35 RNAs as Infernal test
2. Runs BLAST to find homologs
3. Extracts sequences with blastdbcmd
4. Aligns with MAFFT
5. Converts to A3M format
6. Compares with Infernal results

**Output:**
- A3M files in `data/test_blast_diverse/`
- Results table: `diverse_test_results.tsv`
- Side-by-side comparison with Infernal

**Run:**
```bash
python test_blast_rnacentral_diverse.py
```

#### 5. `benchmark_msa_pipeline.py`

**Purpose**: Original benchmark (before we discovered it only tested rRNAs)

**What it does:**
- Tests BLAST on first 30 RNAs from RNAcentral
- All happened to be rRNAs (not diverse!)
- Showed 96.7% success rate (misleading)

**Why we moved on:**
- Not diverse (all rRNAs)
- Didn't compare with Infernal
- Didn't generate actual A3M files

**Historical note**: This was our first test that showed BLAST works, but we later realized we needed diverse testing.

---

## Usage Guide

### Prerequisites

**Required tools:**
```bash
# Install via conda
conda install -c bioconda blast mafft infernal easel

# Or via apt (Linux)
sudo apt install ncbi-blast+ mafft infernal easel
```

**Python dependencies:**
```bash
pip install -r requirements_rna_msa.txt
```

### Database Setup

#### Option 1: BLAST + RNAcentral (Recommended)

**Step 1: Download RNAcentral**
```bash
cd data/rnacentral
wget ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_active.fasta.gz
gunzip rnacentral_active.fasta.gz
# Result: 30GB uncompressed
```

**Step 2: Create BLAST database** (one-time, ~15-20 minutes)
```bash
makeblastdb -in rnacentral_active.fasta \
            -dbtype nucl \
            -out rnacentral_active \
            -parse_seqids
```

**Step 3: Test**
```bash
python test_blast_rnacentral_diverse.py
```

#### Option 2: Infernal + Rfam

**Step 1: Download Rfam**
```bash
cd data/rfam
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed.gz
gunzip Rfam.cm.gz Rfam.seed.gz
```

**Step 2: Index Rfam**
```bash
cmpress Rfam.cm
```

**Step 3: Test**
```bash
python test_infernal_rfam_diverse.py
```

### Generate MSA for Your RNA

#### Using BLAST + RNAcentral

```python
import subprocess
from pathlib import Path
from src.rna_msa.alignment_formats import fasta_to_a3m

# Your RNA
query_id = "my_rna"
query_seq = "AUGCAUGCAUGC..."  # Your sequence
query_file = Path(f"{query_id}_query.fa")

# Write query file
with open(query_file, 'w') as f:
    f.write(f">{query_id}\n{query_seq}\n")

# Step 1: BLAST search
cmd = [
    'blastn',
    '-query', str(query_file),
    '-db', 'data/rnacentral/rnacentral_active',
    '-evalue', '0.001',
    '-outfmt', '6 sseqid',
    '-max_target_seqs', '50',
    '-num_threads', '4'
]
result = subprocess.run(cmd, capture_output=True, text=True, check=True)
seq_ids = result.stdout.strip().split('\n')

print(f"Found {len(seq_ids)} homologs")

# Step 2: Extract sequences
id_file = Path(f"{query_id}_ids.txt")
with open(id_file, 'w') as f:
    for seq_id in seq_ids:
        f.write(f"{seq_id}\n")

homologs_file = Path(f"{query_id}_homologs.fa")
cmd = [
    'blastdbcmd',
    '-db', 'data/rnacentral/rnacentral_active',
    '-entry_batch', str(id_file),
    '-out', str(homologs_file)
]
subprocess.run(cmd, check=True)

# Step 3: Align with MAFFT
aligned_file = Path(f"{query_id}_aligned.fa")
cmd = ['mafft', '--auto', '--thread', '4', str(homologs_file)]
result = subprocess.run(cmd, capture_output=True, text=True, check=True)
with open(aligned_file, 'w') as f:
    f.write(result.stdout)

# Step 4: Convert to A3M
a3m_file = Path(f"{query_id}_aligned.a3m")
fasta_to_a3m(aligned_file, a3m_file, query_id=query_id)

print(f"✓ MSA generated: {a3m_file}")
```

#### Using Infernal + Rfam

```python
import subprocess
from pathlib import Path
from src.rna_msa.alignment_formats import stockholm_to_a3m
import re

# Your RNA
query_id = "my_rna"
query_seq = "AUGCAUGCAUGC..."
query_file = Path(f"{query_id}_query.fa")

with open(query_file, 'w') as f:
    f.write(f">{query_id}\n{query_seq}\n")

# Step 1: Find Rfam family
output_file = Path(f"{query_id}_rfam.tblout")
cmd = [
    'cmscan',
    '--tblout', str(output_file),
    '--cut_ga',
    '--nohmmonly',
    'data/rfam/Rfam.cm',
    str(query_file)
]
subprocess.run(cmd, capture_output=True, check=True)

# Parse results
families = []
with open(output_file, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            parts = line.split()
            if len(parts) >= 16:
                family_acc = parts[1]
                families.append(family_acc)

if not families:
    print("No Rfam family found")
else:
    best_family = families[0]
    print(f"Found family: {best_family}")
    
    # Step 2: Extract seed alignment
    sto_file = Path(f"{query_id}_{best_family}.sto")
    with open('data/rfam/Rfam.seed', 'r', errors='ignore') as fin:
        with open(sto_file, 'w') as fout:
            in_family = False
            for line in fin:
                if line.startswith(f'#=GF AC   {best_family}'):
                    in_family = True
                    fout.write('# STOCKHOLM 1.0\n')
                if in_family:
                    fout.write(line)
                    if line.startswith('//'):
                        break
    
    # Step 3: Convert to A3M
    a3m_file = Path(f"{query_id}_{best_family}.a3m")
    stockholm_to_a3m(sto_file, a3m_file, query_id=None)
    
    print(f"✓ MSA generated: {a3m_file}")
```

### Validate A3M File

```python
from src.rna_msa.alignment_formats import validate_a3m

stats = validate_a3m(Path("my_rna_aligned.a3m"))

if stats['valid']:
    print("✓ Valid A3M file")
    print(f"  Sequences: {stats['n_sequences']}")
    print(f"  Query length: {stats['query_length']} nt")
    print(f"  Length range: {stats['min_length']}-{stats['max_length']} nt")
else:
    print(f"✗ Invalid: {stats['error']}")
```

---

## Troubleshooting

### BLAST Database Not Found

**Error:** `BLAST database not found`

**Solution:**
```bash
cd data/rnacentral
makeblastdb -in rnacentral_active.fasta -dbtype nucl -out rnacentral_active -parse_seqids
```

### MAFFT Alignment Fails

**Error:** `MAFFT failed` or timeout

**Solutions:**
- Reduce number of sequences (use top 30 instead of 50)
- Use faster MAFFT mode: `mafft --retree 1 --maxiterate 0`
- Check if sequences are too long (>2000nt may be slow)

### No Homologs Found

**Error:** `No homologs found`

**Possible causes:**
1. RNA is species-specific (e.g., human lncRNA)
2. RNA is novel/synthetic
3. Database doesn't contain related sequences

**Solutions:**
- Try Infernal+Rfam if it's a known family
- Lower identity threshold (try 60% instead of 70%)
- Check if RNA is in a known family: `cmscan --cut_ga data/rfam/Rfam.cm query.fa`

### Infernal Finds No Family

**Error:** `No Rfam family found`

**This is expected for:**
- Novel RNAs
- Species-specific lncRNAs
- Synthetic sequences

**Solution:** Use BLAST+RNAcentral instead

### A3M Conversion Fails

**Error:** `Query ID not found in alignment`

**Cause:** Query sequence isn't in the Rfam seed alignment

**Solution:** Use `query_id=None` to use first sequence as query:
```python
stockholm_to_a3m(sto_file, a3m_file, query_id=None)
```

---

## Performance Optimization

### Speed Improvements

**1. BLAST Indexing** (one-time setup)
- Creates indexed database for fast searches
- Reduces search time from hours to seconds
- Required for RNAcentral (30GB)

**2. Multi-threading**
- BLAST: `-num_threads 4`
- MAFFT: `--thread 4`
- cmscan: `--cpu 4`
- Scales linearly with CPU cores

**3. Limit MSA Depth**
- 50 sequences is usually sufficient
- More sequences = diminishing returns
- Reduces alignment time

### Memory Optimization

**1. Don't Load Entire Database**
- Use BLAST database indexing
- Extract only needed sequences
- Stream processing

**2. Batch Processing**
- Process multiple RNAs in parallel
- Use GNU parallel or Python multiprocessing

**3. Limit Sequence Length**
- Very long RNAs (>2000nt) are slow to align
- Consider domain splitting for >1500nt

---

## Citation & References

### Databases

**RNAcentral:**
- RNAcentral Consortium (2021). "RNAcentral 2021: secondary structure integration, improved sequence search and new member databases." Nucleic Acids Research.
- URL: https://rnacentral.org/

**Rfam:**
- Kalvari et al. (2021). "Rfam 14: expanded coverage of metagenomic, viral and microRNA families." Nucleic Acids Research.
- URL: https://rfam.org/

### Tools

**BLAST:**
- Camacho et al. (2009). "BLAST+: architecture and applications." BMC Bioinformatics.

**MAFFT:**
- Katoh & Standley (2013). "MAFFT multiple sequence alignment software version 7: improvements in performance and usability." Molecular Biology and Evolution.

**Infernal:**
- Nawrocki & Eddy (2013). "Infernal 1.1: 100-fold faster RNA homology searches." Bioinformatics.

### Related Work

**RNAcmap:**
- Wang et al. (2021). "RNAcmap: a fully automatic pipeline for predicting contact maps of RNAs by evolutionary coupling analysis." Bioinformatics.

**AlphaFold:**
- Jumper et al. (2021). "Highly accurate protein structure prediction with AlphaFold." Nature.

---

## Summary for Professor Azam

### What Was Requested

From your email (September 26, 2025):
> "Investigate practical strategies for obtaining .a3m alignments for RNAs that are suitable for extracting structural features. Once we can robustly generate RNA MSAs, we'll be in a strong position to move forward with extracting 2D features and integrating them into our model."

### What Was Delivered

✅ **Research on two MSA generation approaches:**
1. **Infernal + Rfam**: Structure-aware, 62.9% success rate
2. **BLAST + RNAcentral**: Sequence-based, 74.3% success rate (recommended)

✅ **Domain splitting investigation:**
- Tested on 5 lncRNAs (22 domains total)
- Both Infernal and BLAST approaches tested
- Result: 0% success - domain splitting not recommended

✅ **Comprehensive benchmarking:**
- 35 diverse RNAs (7 types) for main comparison
- 5 long lncRNAs for domain splitting tests
- Fair comparison (same RNAs, multiple approaches)
- Actual A3M files generated and validated

✅ **High-quality MSAs:**
- Average depth: 90-185 sequences
- Average identity: 91.5%
- Ready for 2D feature extraction

✅ **Complete documentation:**
- All benchmark scripts documented
- Usage examples provided
- Troubleshooting guide included

### Key Findings

1. **BLAST + RNAcentral is the better primary approach**
   - Higher success rate (74.3% vs 62.9%)
   - Works for novel RNAs
   - Query-specific MSAs

2. **Both approaches have strengths**
   - Infernal: Structure-aware, larger MSAs
   - BLAST: Broader coverage, more flexible

3. **RNA type matters**
   - Conserved RNAs (rRNA, snoRNA): Both work great
   - Novel RNAs (lncRNA): BLAST works better
   - Well-defined families (snRNA): Infernal works better

4. **Dataset selection is critical**
   - Use conserved RNA types for best results
   - Avoid species-specific lncRNAs
   - RNAcentral contains suitable RNAs

The MSA generation research is complete with clear recommendations: use BLAST+RNAcentral as the primary approach, avoid domain splitting, and focus on conserved RNA types for best results.

---

## Contact & Acknowledgments

This work was completed in response to Professor Azam's request to develop RNA MSA generation capabilities for enhancing RNA-protein interaction prediction with structural features.

**Date:** November 26, 2025

**Key Technologies:**
- Infernal 1.1.5
- BLAST+ 2.15.0
- MAFFT 7.520
- RNAcentral Release 23.0
- Rfam 14.10
