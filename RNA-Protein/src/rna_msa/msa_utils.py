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
    from .feature_extraction import load_msa
    msa = load_msa(msa_path)
    
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
    from .feature_extraction import load_msa
    msa = load_msa(msa_path)
    
    if len(msa) == 0:
        return {'n_sequences': 0, 'length': 0}
    
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
