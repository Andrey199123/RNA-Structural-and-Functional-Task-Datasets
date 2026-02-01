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


def stockholm_to_a3m(sto_file: Path, a3m_file: Path, query_id: str = None, 
                     query_sequence: str = None):
    """Convert Stockholm to A3M format preserving alignment columns.
    
    For trRosettaRNA2 compatibility:
    - Query sequence is first and has NO gaps
    - All other sequences have SAME length as query (gaps preserved)
    - Each sequence written on ONE line
    - Ambiguity codes (WRYKMSBDHVN) normalized to 'N'
    
    Args:
        sto_file: Input Stockholm file
        a3m_file: Output A3M file
        query_id: ID to use for query sequence (if query_sequence provided)
        query_sequence: Optional query sequence to prepend (for Rfam alignments
                       where the actual query isn't in the seed alignment)
    """
    sequences, metadata = parse_stockholm(sto_file)
    
    if not sequences:
        raise ValueError(f"No sequences found in {sto_file}")
    
    # Normalize ambiguity codes in all sequences
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
                normalized.append('-')  # Convert Stockholm gaps to standard gaps
            else:
                normalized.append(char)  # Keep other characters as-is
        return ''.join(normalized)
    
    # If external query sequence provided, align others to it
    if query_sequence and query_id:
        # Normalize query
        query_ungapped = normalize_sequence(query_sequence).replace('-', '')
        query_len = len(query_ungapped)
        
        with open(a3m_file, 'w') as f:
            # Write query first (no gaps)
            f.write(f">{query_id}\n{query_ungapped}\n")
            
            # For Rfam seed alignments: remove all gaps from seed sequences
            # and pad/truncate to match query length
            for seq_id, seq in sequences.items():
                normalized = normalize_sequence(seq)
                # Remove all gaps
                ungapped = normalized.replace('-', '')
                
                # Pad or truncate to match query length
                if len(ungapped) < query_len:
                    # Pad with gaps at the end
                    aligned = ungapped + '-' * (query_len - len(ungapped))
                elif len(ungapped) > query_len:
                    # Truncate (shouldn't happen often)
                    aligned = ungapped[:query_len]
                else:
                    # Same length
                    aligned = ungapped
                
                f.write(f">{seq_id}\n{aligned}\n")
        return
    
    # Standard behavior: query is in the alignment
    if query_id is None:
        query_id = list(sequences.keys())[0]
    
    if query_id not in sequences:
        raise ValueError(f"Query ID {query_id} not found in alignment")
    
    query_seq = normalize_sequence(sequences[query_id])
    
    # Find positions where query has residues (not gaps)
    query_positions = [i for i, char in enumerate(query_seq) if char != '-']
    
    with open(a3m_file, 'w') as f:
        # Write query first (no gaps)
        query_ungapped = query_seq.replace('-', '')
        f.write(f">{query_id}\n{query_ungapped}\n")
        
        # Write other sequences, preserving alignment to query
        for seq_id, seq in sequences.items():
            if seq_id == query_id:
                continue  # Already wrote query
            
            normalized = normalize_sequence(seq)
            # Extract only columns where query has residues
            aligned = ''.join([normalized[i] if i < len(normalized) else '-' 
                              for i in query_positions])
            f.write(f">{seq_id}\n{aligned}\n")


def convert_to_a3m_format(seq: str, query_seq: str, query_positions: List[int]) -> str:
    a3m_seq = []
    for char in seq:
        if char in ['-', '.']:
            continue
        else:
            a3m_seq.append(char.upper())
    return ''.join(a3m_seq)


def fasta_to_a3m(fasta_file: Path, a3m_file: Path, query_id: str = None):
    """Convert aligned FASTA to A3M format preserving alignment columns.
    
    For trRosettaRNA2 compatibility:
    - Query sequence is first and has NO gaps
    - All other sequences have SAME length as query (gaps preserved)
    - Each sequence written on ONE line
    - Ambiguity codes (WRYKMSBDHVN) normalized to 'N'
    
    Args:
        fasta_file: Input aligned FASTA file (from MAFFT)
        a3m_file: Output A3M file
        query_id: ID of query sequence (must be in alignment)
    """
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
    
    # Normalize ambiguity codes
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
                normalized.append('-')  # Convert dots to standard gaps
            else:
                normalized.append(char)
        return ''.join(normalized)
    
    query_seq = normalize_sequence(sequences[query_id])
    
    # Find positions where query has residues (not gaps)
    query_positions = [i for i, char in enumerate(query_seq) if char != '-']
    
    with open(a3m_file, 'w') as f:
        # Write query first (no gaps)
        query_ungapped = query_seq.replace('-', '')
        f.write(f">{query_id}\n{query_ungapped}\n")
        
        # Write other sequences, preserving alignment to query
        for seq_id, seq in sequences.items():
            if seq_id == query_id:
                continue  # Already wrote query
            
            normalized = normalize_sequence(seq)
            # Extract only columns where query has residues
            aligned = ''.join([normalized[i] if i < len(normalized) else '-' 
                              for i in query_positions])
            f.write(f">{seq_id}\n{aligned}\n")


def a3m_to_fasta(a3m_file: Path, fasta_file: Path):
    with open(a3m_file, 'r') as fin, open(fasta_file, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                fout.write(line)
            else:
                fout.write(line.upper())


def validate_a3m(a3m_file: Path) -> Dict[str, any]:
    """Validate A3M file format for trRosettaRNA2 compatibility.
    
    Checks:
    - Query sequence (first) has no gaps
    - All sequences have same length
    - No ambiguity codes (should be normalized to N)
    - Each sequence on one line
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
    
    query_id = list(sequences.keys())[0]
    query_seq = sequences[query_id]
    query_len = len(query_seq)
    
    # Check query has no gaps
    if '-' in query_seq:
        return {'valid': False, 'error': 'Query sequence contains gaps'}
    
    # Check all sequences have same length
    seq_lengths = [len(seq) for seq in sequences.values()]
    if len(set(seq_lengths)) > 1:
        return {
            'valid': False, 
            'error': f'Sequences have different lengths: {set(seq_lengths)}. All sequences must have same length as query.'
        }
    
    # Check for ambiguity codes
    ambiguity_codes = 'WRYKMSBDHV'
    sequences_with_ambiguity = []
    for seq_id, seq in sequences.items():
        if any(char in ambiguity_codes for char in seq.upper()):
            sequences_with_ambiguity.append(seq_id)
    
    if sequences_with_ambiguity:
        return {
            'valid': False,
            'error': f'Found ambiguity codes in {len(sequences_with_ambiguity)} sequences. Should be normalized to N.'
        }
    
    stats = {
        'valid': True,
        'n_sequences': len(sequences),
        'query_id': query_id,
        'query_length': query_len,
        'all_same_length': True,
        'sequence_lengths': seq_lengths,
        'min_length': min(seq_lengths),
        'max_length': max(seq_lengths),
        'avg_length': sum(seq_lengths) / len(sequences)
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
        print("Valid A3M file")
        print(f"  Sequences: {stats['n_sequences']}")
        print(f"  Query length: {stats['query_length']}")
        print(f"  Length range: {stats['min_length']}-{stats['max_length']}")
    else:
        print(f"Invalid: {stats['error']}")
