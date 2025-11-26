#!/usr/bin/env python3
"""
Alignment format conversion utilities
Converts between Stockholm, FASTA, and A3M formats
"""

from pathlib import Path
from typing import Dict, List, Tuple
import re


def parse_stockholm(sto_file: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Parse Stockholm format alignment
    
    Args:
        sto_file: Path to Stockholm format file
        
    Returns:
        Tuple of (sequences dict, metadata dict)
    """
    sequences = {}
    metadata = {}
    
    with open(sto_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            
            # Skip empty lines and end marker
            if not line or line == '//':
                continue
            
            # Parse metadata
            if line.startswith('#'):
                if line.startswith('#=GF'):
                    parts = line.split(None, 2)
                    if len(parts) >= 3:
                        key = parts[1]
                        value = parts[2]
                        metadata[key] = value
                continue
            
            # Parse sequence line
            parts = line.split(None, 1)
            if len(parts) == 2:
                seq_id, seq = parts
                if seq_id in sequences:
                    sequences[seq_id] += seq
                else:
                    sequences[seq_id] = seq
    
    return sequences, metadata


def stockholm_to_a3m(sto_file: Path, a3m_file: Path, query_id: str = None):
    """
    Convert Stockholm alignment to A3M format
    
    A3M format rules:
    - First sequence is the query (reference)
    - Uppercase letters = aligned positions (match states)
    - Lowercase letters = insertions relative to query
    - Gaps in query are removed from all sequences
    
    Args:
        sto_file: Input Stockholm file
        a3m_file: Output A3M file
        query_id: ID of query sequence (if None, uses first sequence)
    """
    # Parse Stockholm file
    sequences, metadata = parse_stockholm(sto_file)
    
    if not sequences:
        raise ValueError(f"No sequences found in {sto_file}")
    
    # Identify query sequence
    if query_id is None:
        query_id = list(sequences.keys())[0]
    
    if query_id not in sequences:
        raise ValueError(f"Query ID {query_id} not found in alignment")
    
    query_seq = sequences[query_id]
    
    # Find positions where query has gaps
    query_positions = []
    for i, char in enumerate(query_seq):
        if char not in ['-', '.']:
            query_positions.append(i)
    
    # Convert alignment to A3M format
    a3m_sequences = {}
    
    for seq_id, seq in sequences.items():
        # Extract only positions where query is not gapped
        aligned_seq = ''.join([seq[i] for i in query_positions])
        
        # For query: remove all gaps, keep uppercase
        if seq_id == query_id:
            a3m_seq = aligned_seq.replace('-', '').replace('.', '').upper()
        else:
            # For other sequences: mark insertions as lowercase
            a3m_seq = convert_to_a3m_format(aligned_seq, query_seq, query_positions)
        
        a3m_sequences[seq_id] = a3m_seq
    
    # Write A3M file
    with open(a3m_file, 'w') as f:
        for seq_id, seq in a3m_sequences.items():
            f.write(f">{seq_id}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')


def convert_to_a3m_format(seq: str, query_seq: str, query_positions: List[int]) -> str:
    """
    Convert aligned sequence to A3M format
    
    Rules:
    - Match to query = uppercase
    - Gap in query = lowercase (insertion)
    - Gap in sequence = remove
    """
    a3m_seq = []
    
    for char in seq:
        if char in ['-', '.']:
            # Gap in this sequence - skip
            continue
        else:
            # Has a residue - keep as uppercase (match state)
            a3m_seq.append(char.upper())
    
    return ''.join(a3m_seq)


def fasta_to_a3m(fasta_file: Path, a3m_file: Path, query_id: str = None):
    """
    Convert FASTA alignment to A3M format
    
    Args:
        fasta_file: Input FASTA file (aligned)
        a3m_file: Output A3M file
        query_id: ID of query sequence (if None, uses first sequence)
    """
    # Parse FASTA
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
    
    # Identify query
    if query_id is None:
        query_id = list(sequences.keys())[0]
    
    if query_id not in sequences:
        raise ValueError(f"Query ID {query_id} not found in alignment")
    
    query_seq = sequences[query_id]
    
    # Find non-gap positions in query
    query_positions = [i for i, char in enumerate(query_seq) if char not in ['-', '.']]
    
    # Convert to A3M
    with open(a3m_file, 'w') as f:
        for seq_id, seq in sequences.items():
            # Extract positions where query is not gapped
            aligned_seq = ''.join([seq[i] for i in query_positions])
            
            if seq_id == query_id:
                # Query: remove gaps, uppercase
                a3m_seq = aligned_seq.replace('-', '').replace('.', '').upper()
            else:
                # Others: remove gaps, uppercase
                a3m_seq = aligned_seq.replace('-', '').replace('.', '').upper()
            
            f.write(f">{seq_id}\n")
            for i in range(0, len(a3m_seq), 80):
                f.write(a3m_seq[i:i+80] + '\n')


def a3m_to_fasta(a3m_file: Path, fasta_file: Path):
    """
    Convert A3M format to FASTA (removes case information)
    
    Args:
        a3m_file: Input A3M file
        fasta_file: Output FASTA file
    """
    with open(a3m_file, 'r') as fin, open(fasta_file, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                fout.write(line)
            else:
                # Convert all to uppercase
                fout.write(line.upper())


def validate_a3m(a3m_file: Path) -> Dict[str, any]:
    """
    Validate A3M format and return statistics
    
    Args:
        a3m_file: A3M file to validate
        
    Returns:
        Dictionary with validation results and statistics
    """
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
    
    # Get query (first sequence)
    query_id = list(sequences.keys())[0]
    query_seq = sequences[query_id]
    query_len = len(query_seq)
    
    # Check query has no lowercase (no insertions)
    if query_seq != query_seq.upper():
        return {'valid': False, 'error': 'Query sequence contains lowercase letters'}
    
    # Statistics
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
    # Test conversion
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
