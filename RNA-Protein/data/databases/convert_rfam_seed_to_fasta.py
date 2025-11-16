#!/usr/bin/env python3
"""Convert Rfam.seed (Stockholm) to FASTA format"""

def stockholm_to_fasta(seed_file, output_file):
    """Convert Stockholm format to FASTA"""
    
    sequences = {}
    current_family = None
    
    print(f"Reading {seed_file}...")
    
    # Try different encodings
    for encoding in ['utf-8', 'latin-1', 'iso-8859-1']:
        try:
            with open(seed_file, encoding=encoding, errors='ignore') as f:
                for line in f:
                    line = line.rstrip()
                    
                    # New family
                    if line.startswith('# STOCKHOLM'):
                        continue
                    elif line.startswith('#=GF ID'):
                        parts = line.split()
                        if len(parts) >= 3:
                            current_family = parts[2]
                    elif line.startswith('//'):
                        current_family = None
                    elif line.startswith('#'):
                        continue
                    elif line and current_family and not line.startswith('#'):
                        # Sequence line
                        parts = line.split()
                        if len(parts) >= 2:
                            seq_id = f"{current_family}_{parts[0]}"
                            seq = parts[1].replace('.', '').replace('-', '').upper()
                            # Convert T to U for RNA
                            seq = seq.replace('T', 'U')
                            # Only keep valid RNA nucleotides
                            seq = ''.join(c for c in seq if c in 'ACGU')
                            if seq and len(seq) > 10:  # Skip very short sequences
                                sequences[seq_id] = seq
            break  # Success, exit encoding loop
        except UnicodeDecodeError:
            print(f"Failed with {encoding}, trying next...")
            continue
    
    # Write FASTA
    print(f"Writing {len(sequences)} sequences to {output_file}...")
    with open(output_file, 'w') as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n{seq}\n")
    
    print(f"✓ Converted {len(sequences)} sequences to {output_file}")
    
    # Calculate actual file size
    import os
    size_mb = os.path.getsize(output_file) / 1024 / 1024
    print(f"✓ File size: {size_mb:.1f} MB")

if __name__ == '__main__':
    stockholm_to_fasta('Rfam.seed', 'rfam_sequences.fasta')
