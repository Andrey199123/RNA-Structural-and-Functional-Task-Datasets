#!/usr/bin/env python3
"""Test RNA generation with just 2 RNAs"""

import subprocess
import sys
from pathlib import Path
from Bio import SeqIO

# Add src to path
sys.path.insert(0, 'RNA-Protein/src')
from rna_msa.alignment_formats import stockholm_to_a3m, fasta_to_a3m, validate_a3m

# Paths
RFAM_CM = Path('RNA-Protein/data/rfam/Rfam.cm')
RFAM_SEED = Path('RNA-Protein/data/rfam/Rfam.seed')
RNACENTRAL_DB = Path('RNA-Protein/data/rnacentral/rnacentral_active')
TRROSETTA_DIR = Path('trRosettaRNA2')

# Setup
output_dir = Path("RNA-Protein/data/test_rna_generation")
msa_dir = output_dir / "msas"
pdb_dir = output_dir / "structures"
msa_dir.mkdir(parents=True, exist_ok=True)
pdb_dir.mkdir(parents=True, exist_ok=True)

# Get 2 short RNAs (under 500 nt for faster testing)
all_rnas = list(SeqIO.parse('RNA-Protein/data/RPI7317/RPI7317_IncRNA_sequence.fasta', 'fasta'))
rnas = [r for r in all_rnas if 200 <= len(r.seq) <= 500][:2]

print("=" * 80)
print(f"Testing RNA Structure Generation with {len(rnas)} RNAs")
print("=" * 80)
print()

for i, rna in enumerate(rnas, 1):
    rna_id = rna.id
    sequence = str(rna.seq)
    
    print(f"[{i}/{len(rnas)}] Processing: {rna_id} ({len(sequence)} nt)")
    
    # Write query
    query_file = msa_dir / f"{rna_id}_query.fa"
    query_file.write_text(f">{rna_id}\n{sequence}\n")
    print(f"  ✓ Query file written")
    
    # Try Infernal
    print(f"  [1] Trying Infernal+Rfam...", end='', flush=True)
    msa_success = False
    method = None
    a3m_file = None
    
    try:
        infernal_out = msa_dir / f"{rna_id}_infernal.tblout"
        cmd = [
            'cmscan',
            '--tblout', str(infernal_out),
            '--cut_ga',
            '--nohmmonly',
            '--cpu', '2',
            str(RFAM_CM),
            str(query_file)
        ]
        
        result = subprocess.run(cmd, capture_output=True, timeout=60)
        print(f" ran", end='', flush=True)
        
        # Parse for families
        families = []
        with open(infernal_out) as f:
            for line in f:
                if not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 16:
                        families.append(parts[1])
        
        if families:
            best_acc = families[0]
            print(f" found {best_acc}", end='', flush=True)
            
            # Extract seed alignment
            sto_file = msa_dir / f"{rna_id}_{best_acc}.sto"
            with open(RFAM_SEED, errors='ignore') as fin:
                with open(sto_file, 'w') as fout:
                    in_family = False
                    for line in fin:
                        if line.startswith(f'#=GF AC   {best_acc}'):
                            in_family = True
                            fout.write('# STOCKHOLM 1.0\n')
                        if in_family:
                            fout.write(line)
                            if line.startswith('//'):
                                break
            
            # Convert to A3M
            a3m_file = msa_dir / f"{rna_id}.a3m"
            stockholm_to_a3m(sto_file, a3m_file, query_id=None)
            
            stats = validate_a3m(a3m_file)
            if stats['valid']:
                msa_success = True
                method = "infernal"
                print(f" ✓ Success ({stats['num_sequences']} seqs)")
            else:
                print(f" ✗ Invalid A3M")
        else:
            print(f" no family found")
    except Exception as e:
        print(f" ✗ Error: {e}")
    
    # Try BLAST if Infernal failed
    if not msa_success:
        print(f"  [2] Trying BLAST+RNAcentral...", end='', flush=True)
        try:
            blast_out = msa_dir / f"{rna_id}_blast.txt"
            cmd = [
                'blastn',
                '-query', str(query_file),
                '-db', str(RNACENTRAL_DB),
                '-evalue', '0.001',
                '-outfmt', '6 sseqid pident',
                '-max_target_seqs', '100',
                '-num_threads', '2'
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            
            # Get homologs
            homolog_ids = set()
            for line in result.stdout.strip().split('\n'):
                if line:
                    parts = line.split('\t')
                    if len(parts) >= 2 and float(parts[1]) >= 70:
                        homolog_ids.add(parts[0])
            
            print(f" {len(homolog_ids)} homologs", end='', flush=True)
            
            if len(homolog_ids) >= 10:
                # Extract sequences
                id_file = msa_dir / f"{rna_id}_ids.txt"
                id_file.write_text('\n'.join(homolog_ids))
                
                homologs_file = msa_dir / f"{rna_id}_homologs.fa"
                cmd = [
                    'blastdbcmd',
                    '-db', str(RNACENTRAL_DB),
                    '-entry_batch', str(id_file),
                    '-line_length', '100000',
                    '-out', str(homologs_file)
                ]
                subprocess.run(cmd, capture_output=True, timeout=60)
                
                # Add query sequence to homologs
                with open(homologs_file, 'r') as f:
                    homolog_content = f.read()
                with open(homologs_file, 'w') as f:
                    f.write(f">{rna_id}\n{sequence}\n")
                    f.write(homolog_content)
                
                # Align with MAFFT
                cmd = ['mafft', '--auto', '--thread', '2', str(homologs_file)]
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
                
                aligned_file = msa_dir / f"{rna_id}_aligned.fa"
                aligned_file.write_text(result.stdout)
                
                # Convert to A3M
                a3m_file = msa_dir / f"{rna_id}.a3m"
                fasta_to_a3m(aligned_file, a3m_file, query_id=rna_id)
                
                stats = validate_a3m(a3m_file)
                if stats['valid']:
                    msa_success = True
                    method = "blast"
                    print(f" ✓ Success ({stats['n_sequences']} seqs)")
                else:
                    print(f" ✗ Invalid A3M")
            else:
                print(f" ✗ Insufficient homologs")
        except Exception as e:
            print(f" ✗ Error: {e}")
    
    # Predict structure if MSA succeeded
    structure_success = False
    if msa_success and a3m_file and a3m_file.exists():
        print(f"  [3] Predicting structure with trRosettaRNA2...", end='', flush=True)
        try:
            # Use absolute paths since we change working directory
            abs_a3m = a3m_file.resolve()
            abs_out = (pdb_dir / rna_id).resolve()
            abs_params = (TRROSETTA_DIR / 'params').resolve()
            
            cmd = [
                'python', '-m', 'trRNA2.predict',
                '-mdir', str(abs_params),
                '-i', str(abs_a3m),
                '-o', str(abs_out),
                '-relax_steps', '0'
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600, cwd=str(TRROSETTA_DIR))
            
            if result.returncode != 0:
                print(f" ✗ trRosetta failed: {result.stderr[:100]}")
            else:
                # trRosettaRNA2 outputs to model_1_unrelaxed.pdb
                predicted_pdb = pdb_dir / rna_id / "model_1_unrelaxed.pdb"
                pdb_file = pdb_dir / f"{rna_id}.pdb"
                
                if predicted_pdb.exists():
                    pdb_file.write_text(predicted_pdb.read_text())
                    structure_success = True
                    print(f" ✓ Success")
                else:
                    print(f" ✗ PDB not found at {predicted_pdb}")
        except subprocess.TimeoutExpired:
            print(f" ✗ Timeout (>10 min)")
        except Exception as e:
            print(f" ✗ Error: {e}")
    
    # Summary for this RNA
    if structure_success:
        print(f"  ✅ COMPLETE: MSA ({method}) + Structure")
    elif msa_success:
        print(f"  ⚠️  PARTIAL: MSA only ({method})")
    else:
        print(f"  ❌ FAILED: No MSA generated")
    
    print()

print("=" * 80)
print("Test complete! Check output in:", output_dir)
print("=" * 80)
