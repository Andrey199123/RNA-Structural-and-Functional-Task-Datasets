#!/usr/bin/env python3
"""
Utility functions for MSA processing
"""

import logging
from pathlib import Path
import subprocess
from typing import List, Tuple
import numpy as np

logger = logging.getLogger(__name__)


def clean_msa(msa_path: Path, output_path: Path, 
              max_identity: float = 0.95,
              min_coverage: float = 0.5,
              max_gap_fraction: float = 0.8) -> Path:
    """
    Clean MSA by removing redundant and low-quality sequences
    
    Args:
        msa_path: Input MSA file
        output_path: Output cleaned MSA file
        max_identity: Maximum sequence identity threshold
        min_coverage: Minimum coverage of query sequence
        max_gap_fraction: Maximum fraction of gaps allowed
        
    Returns:
        Path to cleaned MSA
    """
    logger.info(f"Cleaning MSA: {msa_path}")
    
    # Load MSA
    # TODO: Implement load_msa in feature_extraction module
    # For now, skip MSA loading
    logger.warning("MSA loading not implemented yet, skipping clean_msa")
    return msa_path
    
    if len(msa) == 0:
        raise ValueError("Empty MSA")
    
    # Keep query (first sequence)
    query = msa[0]
    cleaned = [query]
    
    # Filter sequences
    for i in range(1, len(msa)):
        seq = msa[i]
        
        # Check coverage
        non_gap = (seq < 4).sum()
        coverage = non_gap / len(seq)
        if coverage < min_coverage:
            continue
        
        # Check gap fraction
        gap_frac = (seq == 4).sum() / len(seq)
        if gap_frac > max_gap_fraction:
            continue
        
        # Check identity to query
        identity = compute_identity(query, seq)
        if identity > max_identity:
            continue
        
        # Check identity to already selected sequences
        too_similar = False
        for selected in cleaned[1:]:
            if compute_identity(selected, seq) > max_identity:
                too_similar = True
                break
        
        if not too_similar:
            cleaned.append(seq)
    
    logger.info(f"Cleaned MSA: {len(msa)} -> {len(cleaned)} sequences")
    
    # Save cleaned MSA
    save_msa(cleaned, output_path)
    
    return output_path


def compute_identity(seq1: np.ndarray, seq2: np.ndarray) -> float:
    """Compute sequence identity between two sequences"""
    # Only compare non-gap positions
    mask = (seq1 < 4) & (seq2 < 4)
    if mask.sum() == 0:
        return 0.0
    
    matches = (seq1[mask] == seq2[mask]).sum()
    return matches / mask.sum()


def save_msa(msa: List[np.ndarray], output_path: Path, format: str = "fasta"):
    """
    Save MSA to file
    
    Args:
        msa: List of sequence arrays
        output_path: Output file path
        format: Output format (fasta, a3m, stockholm)
    """
    decode = {0: 'A', 1: 'C', 2: 'G', 3: 'U', 4: '-', 5: 'N'}
    
    with open(output_path, 'w') as f:
        if format == "stockholm":
            f.write("# STOCKHOLM 1.0\n")
        
        for i, seq in enumerate(msa):
            seq_str = ''.join(decode.get(nt, 'N') for nt in seq)
            
            if format == "fasta" or format == "a3m":
                f.write(f">seq{i}\n{seq_str}\n")
            elif format == "stockholm":
                f.write(f"seq{i:<30} {seq_str}\n")
        
        if format == "stockholm":
            f.write("//\n")
    
    logger.info(f"Saved MSA to {output_path}")


def remove_gappy_columns(msa: np.ndarray, max_gap_fraction: float = 0.5) -> np.ndarray:
    """
    Remove columns with too many gaps
    
    Args:
        msa: MSA array
        max_gap_fraction: Maximum fraction of gaps allowed per column
        
    Returns:
        Filtered MSA
    """
    n_seq, length = msa.shape
    
    # Compute gap fraction per column
    gap_fractions = (msa == 4).sum(axis=0) / n_seq
    
    # Keep columns with acceptable gap fraction
    keep_cols = gap_fractions <= max_gap_fraction
    
    filtered_msa = msa[:, keep_cols]
    
    logger.info(f"Removed gappy columns: {length} -> {filtered_msa.shape[1]} positions")
    
    return filtered_msa


def realign_msa(msa_path: Path, output_path: Path, method: str = "mafft") -> Path:
    """
    Realign MSA using structure-aware aligner
    
    Args:
        msa_path: Input MSA file
        output_path: Output aligned MSA file
        method: Alignment method (mafft, muscle, etc.)
        
    Returns:
        Path to realigned MSA
    """
    logger.info(f"Realigning MSA with {method}")
    
    if method == "mafft":
        cmd = [
            "mafft",
            "--maxiterate", "1000",
            "--localpair",
            "--thread", "4",
            str(msa_path)
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            with open(output_path, 'w') as f:
                f.write(result.stdout)
            
            logger.info(f"Realigned MSA saved to {output_path}")
            return output_path
            
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            logger.error(f"MAFFT failed: {e}")
            # Copy input as fallback
            import shutil
            shutil.copy(msa_path, output_path)
            return output_path
    
    else:
        raise ValueError(f"Unknown alignment method: {method}")


def compute_msa_statistics(msa_path: Path) -> dict:
    """
    Compute statistics for an MSA
    
    Returns:
        Dictionary with MSA statistics
    """
    # TODO: Implement load_msa in feature_extraction module
    # For now, return basic stats
    logger.warning("MSA loading not implemented yet, returning basic stats")
    msa = None
    
    if msa is None or len(msa) == 0:
        # Count sequences manually from file
        n_seq = 0
        with open(msa_path) as f:
            for line in f:
                if line.startswith('>') or (not line.startswith('#') and not line.startswith('//') and line.strip()):
                    n_seq += 1
        return {'n_sequences': max(1, n_seq // 2), 'length': 0}  # Rough estimate
    
    n_seq, length = msa.shape
    
    # Compute statistics
    stats = {
        'n_sequences': n_seq,
        'length': length,
        'avg_identity': 0.0,
        'avg_coverage': 0.0,
        'avg_gap_fraction': 0.0
    }
    
    if n_seq > 1:
        query = msa[0]
        identities = []
        coverages = []
        gap_fractions = []
        
        for i in range(1, n_seq):
            seq = msa[i]
            
            # Identity to query
            identities.append(compute_identity(query, seq))
            
            # Coverage
            non_gap = (seq < 4).sum()
            coverages.append(non_gap / length)
            
            # Gap fraction
            gap_fractions.append((seq == 4).sum() / length)
        
        stats['avg_identity'] = np.mean(identities)
        stats['avg_coverage'] = np.mean(coverages)
        stats['avg_gap_fraction'] = np.mean(gap_fractions)
    
    return stats


def download_rfam_database(output_dir: Path) -> Path:
    """
    Download Rfam covariance model database
    
    Args:
        output_dir: Directory to save database
        
    Returns:
        Path to Rfam.cm file
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    rfam_cm = output_dir / "Rfam.cm"
    
    if rfam_cm.exists():
        logger.info(f"Rfam database already exists: {rfam_cm}")
        return rfam_cm
    
    logger.info("Downloading Rfam database...")
    
    # Download from Rfam FTP
    import urllib.request
    import gzip
    
    url = "ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"
    gz_file = output_dir / "Rfam.cm.gz"
    
    try:
        urllib.request.urlretrieve(url, gz_file)
        logger.info(f"Downloaded Rfam.cm.gz")
        
        # Decompress
        with gzip.open(gz_file, 'rb') as f_in:
            with open(rfam_cm, 'wb') as f_out:
                f_out.write(f_in.read())
        
        logger.info(f"Decompressed to {rfam_cm}")
        
        # Remove gz file
        gz_file.unlink()
        
        # Press CM database for faster searching
        logger.info("Pressing Rfam database...")
        subprocess.run(["cmpress", str(rfam_cm)], check=True)
        
        return rfam_cm
        
    except Exception as e:
        logger.error(f"Failed to download Rfam database: {e}")
        raise


def convert_stockholm_to_a3m(sto_path: Path, a3m_path: Path, query_id: str = None) -> Path:
    """
    Convert Stockholm alignment to A3M format
    
    A3M format:
    - First sequence is query (no gaps)
    - Other sequences: '-' for deletions, lowercase for insertions
    
    Args:
        sto_path: Input Stockholm file
        a3m_path: Output A3M file
        query_id: ID of query sequence (if None, uses first sequence)
        
    Returns:
        Path to A3M file
    """
    logger.info(f"Converting {sto_path} to A3M format")
    
    # Read Stockholm alignment
    sequences = {}
    seq_order = []
    
    with open(sto_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or line.startswith('//') or not line:
                continue
            
            parts = line.split()
            if len(parts) >= 2:
                seq_id = parts[0]
                seq = parts[1]
                
                if seq_id not in sequences:
                    sequences[seq_id] = []
                    seq_order.append(seq_id)
                sequences[seq_id].append(seq)
    
    # Concatenate sequence fragments
    for seq_id in sequences:
        sequences[seq_id] = ''.join(sequences[seq_id])
    
    if not sequences:
        raise ValueError(f"No sequences found in {sto_path}")
    
    # Determine query
    if query_id is None or query_id not in sequences:
        query_id = seq_order[0]
    
    query_seq = sequences[query_id]
    
    # Write A3M format
    with open(a3m_path, 'w') as f:
        # Write query without gaps
        query_ungapped = query_seq.replace('-', '').replace('.', '')
        f.write(f">{query_id}\n{query_ungapped}\n")
        
        # Write other sequences in A3M format
        for seq_id in seq_order:
            if seq_id == query_id:
                continue
            
            seq = sequences[seq_id]
            
            # Convert to A3M format
            a3m_seq = []
            for q_char, s_char in zip(query_seq, seq):
                if q_char in '-.':  # Query has gap = insertion in subject
                    if s_char not in '-.':
                        a3m_seq.append(s_char.lower())
                else:  # Query has residue
                    if s_char in '-.':
                        a3m_seq.append('-')
                    else:
                        a3m_seq.append(s_char.upper())
            
            f.write(f">{seq_id}\n{''.join(a3m_seq)}\n")
    
    logger.info(f"Converted to A3M: {len(sequences)} sequences")
    return a3m_path


def check_dependencies() -> dict:
    """
    Check if required external tools are installed
    
    Returns:
        Dictionary with tool availability
    """
    tools = {
        'infernal': ['cmscan', '--version'],
        'mafft': ['mafft', '--version'],
        'RNAfold': ['RNAfold', '--version'],
        'python_RNA': None  # Check via import
    }
    
    status = {}
    
    for tool, cmd in tools.items():
        if cmd is None:
            # Check Python module
            try:
                if tool == 'python_RNA':
                    import RNA
                status[tool] = True
            except ImportError:
                status[tool] = False
        else:
            # Check command-line tool
            try:
                result = subprocess.run(cmd, capture_output=True, text=True)
                # Tool exists if it runs (even if exit code is non-zero for --version)
                status[tool] = True
            except FileNotFoundError:
                status[tool] = False
    
    return status
