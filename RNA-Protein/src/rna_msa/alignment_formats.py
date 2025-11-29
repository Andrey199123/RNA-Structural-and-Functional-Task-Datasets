#!/usr/bin/env python3
"""
Converts between Stockholm, FASTA, and A3M alignment formats.
"""

from pathlib import Path
from typing import Dict, List, Tuple
import re


def parse_stockholm(sto_file: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    sequences = {}
    metadata = {}
    
    with open(sto_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if not line or line == '//':
                continue
            
            if line.startswith('#'):
                if line.startswith('#=GF'):
                    parts = line.split(None, 2)
                    if len(parts) >= 3:
                        key = parts[1]
                        value = parts[2]
                        metadata[key] = value
                continue
            
            parts = line.split(None, 1)
            if len(parts) == 2:
                seq_id, seq = parts
                if seq_id in sequences:
                    sequences[seq_id] += seq
                else:
                    sequences[seq_id] = seq
    
    return sequences, metadata


def stockholm_to_a3m(sto_file: Path, a3m_file: Path, query_id: str = None):
    sequences, metadata = parse_stockholm(sto_file)
    
    if not sequences:
        raise ValueError(f"No sequences found in {sto_file}")
    
    if query_id is None:
        query_id = list(sequences.keys())[0]
    
    if query_id not in sequences:
        raise ValueError(f"Query ID {query_id} not found in alignment")
    
    query_seq = sequences[query_id]
    
    query_positions = []
    for i, char in enumerate(query_seq):
        if char not in ['-', '.']:
            query_positions.append(i)
    
    a3m_sequences = {}
    
    for seq_id, seq in sequences.items():
        aligned_seq = ''.join([seq[i] for i in query_positions])
        
        if seq_id == query_id:
            a3m_seq = aligned_seq.replace('-', '').replace('.', '').upper()
        else:
            a3m_seq = convert_to_a3m_format(aligned_seq, query_seq, query_positions)
        
        a3m_sequences[seq_id] = a3m_seq
    
    with open(a3m_file, 'w') as f:
        for seq_id, seq in a3m_sequences.items():
            f.write(f">{seq_id}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')


def convert_to_a3m_format(seq: str, query_seq: str, query_positions: List[int]) -> str:
    a3m_seq = []
    for char in seq:
        if char in ['-', '.']:
            continue
        else:
            a3m_seq.append(char.upper())
    return ''.join(a3m_seq)


def fasta_to_a3m(fasta_file: Path, a3m_file: Path, query_id: str = None):
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    if not sequences:
        raise ValueError(f"No sequences found in {fasta_file}")
    
    if query_id is None:
        query_id = list(sequences.keys())[0]
    
    if query_id not in sequences:
        raise ValueError(f"Query ID {query_id} not found in alignment")
    
    query_seq = sequences[query_id]
    query_positions = [i for i, char in enumerate(query_seq) if char not in ['-', '.']]
    
    with open(a3m_file, 'w') as f:
        for seq_id, seq in sequences.items():
            aligned_seq = ''.join([seq[i] for i in query_positions])
            
            if seq_id == query_id:
                a3m_seq = aligned_seq.replace('-', '').replace('.', '').upper()
            else:
                a3m_seq = aligned_seq.replace('-', '').replace('.', '').upper()
            
            f.write(f">{seq_id}\n")
            for i in range(0, len(a3m_seq), 80):
                f.write(a3m_seq[i:i+80] + '\n')


def a3m_to_fasta(a3m_file: Path, fasta_file: Path):
    with open(a3m_file, 'r') as fin, open(fasta_file, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                fout.write(line)
            else:
                fout.write(line.upper())


def validate_a3m(a3m_file: Path) -> Dict[str, any]:
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(a3m_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    if not sequences:
        return {'valid': False, 'error': 'No sequences found'}
    
    query_id = list(sequences.keys())[0]
    query_seq = sequences[query_id]
    query_len = len(query_seq)
    
    if query_seq != query_seq.upper():
        return {'valid': False, 'error': 'Query sequence contains lowercase letters'}
    
    stats = {
        'valid': True,
        'n_sequences': len(sequences),
        'query_id': query_id,
        'query_length': query_len,
        'sequence_lengths': [len(seq) for seq in sequences.values()],
        'min_length': min(len(seq) for seq in sequences.values()),
        'max_length': max(len(seq) for seq in sequences.values()),
        'avg_length': sum(len(seq) for seq in sequences.values()) / len(sequences)
    }
    
    return stats


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 3:
        print("Usage: python alignment_formats.py <input.sto> <output.a3m>")
        sys.exit(1)
    
    sto_file = Path(sys.argv[1])
    a3m_file = Path(sys.argv[2])
    
    print(f"Converting {sto_file} to {a3m_file}...")
    stockholm_to_a3m(sto_file, a3m_file)
    
    print("Validating A3M file...")
    stats = validate_a3m(a3m_file)
    
    if stats['valid']:
        print("✓ Valid A3M file")
        print(f"  Sequences: {stats['n_sequences']}")
        print(f"  Query length: {stats['query_length']}")
        print(f"  Length range: {stats['min_length']}-{stats['max_length']}")
    else:
        print(f"✗ Invalid: {stats['error']}")
