# RNA-Protein Databases

This directory contains database files needed for RNA MSA generation.

## Required Files

### Rfam Database
Download from: ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/

Required files:
- `Rfam.cm` - Covariance models for all Rfam families (~311 MB)
- `Rfam.seed` - Seed alignments for all families (~1.4 GB)

After downloading `Rfam.cm`, press and index it:
```bash
cmpress Rfam.cm
```

This will create index files:
- `Rfam.cm.i1f`
- `Rfam.cm.i1i`
- `Rfam.cm.i1m`
- `Rfam.cm.i1p`

### Optional: Sequence Database
For homology search (RNAcmap-style strategy), you can use:
- RNAcentral database
- Custom RNA sequence database

## Setup Script

See `SETUP_DATABASES.md` in the project root for automated setup instructions.
