#!/usr/bin/env python3
"""
2D Feature Extraction from RNA MSAs
Generates image-like features for vision-based models
"""

import logging
from pathlib import Path
from typing import Dict, Tuple
import numpy as np

logger = logging.getLogger(__name__)


def extract_2d_features(msa_path: Path, output_path: Path = None) -> Dict[str, np.ndarray]:
    """
    Extract 2D feature maps from RNA MSA
    
    Generates three types of features:
    1. Conservation map (per-position)
    2. Coevolution/contact map (pairwise)
    3. Pairing/accessibility map
    
    Args:
        msa_path: Path to MSA file (.a3m or Stockholm format)
        output_path: Optional path to save features
        
    Returns:
        Dictionary with feature arrays
    """
    logger.info(f"Extracting 2D features from {msa_path}")
    
    # Load MSA
    msa = load_msa(msa_path)
    
    if len(msa) == 0:
        raise ValueError(f"Empty MSA: {msa_path}")
    
    # Extract features
    conservation = compute_conservation_map(msa)
    coevolution = compute_coevolution_map(msa)
    accessibility = compute_accessibility_map(msa)
    
    features = {
        'conservation': conservation,
        'coevolution': coevolution,
        'accessibility': accessibility,
        'msa_depth': len(msa)
    }
    
    # Save if requested
    if output_path:
        np.savez_compressed(output_path, **features)
        logger.info(f"Saved features to {output_path}")
    
    return features


def load_msa(msa_path: Path) -> np.ndarray:
    """
    Load MSA from file
    
    Args:
        msa_path: Path to MSA file
        
    Returns:
        MSA as numpy array (n_sequences x length)
    """
    # Ensure msa_path is a Path object
    msa_path = Path(msa_path)
    suffix = msa_path.suffix.lower()
    
    if suffix == '.a3m':
        return load_a3m(msa_path)
    elif suffix in ['.sto', '.stockholm']:
        return load_stockholm(msa_path)
    elif suffix in ['.fasta', '.fa']:
        return load_fasta_msa(msa_path)
    else:
        # Try to auto-detect
        try:
            return load_a3m(msa_path)
        except:
            try:
                return load_stockholm(msa_path)
            except:
                return load_fasta_msa(msa_path)


def load_a3m(a3m_path: Path) -> np.ndarray:
    """Load A3M format MSA"""
    sequences = []
    current_seq = []
    
    with open(a3m_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append(''.join(current_seq))
                current_seq = []
            else:
                # Remove lowercase (insertions) and convert to uppercase
                seq = ''.join(c.upper() if c.isupper() or c == '-' else '' for c in line)
                current_seq.append(seq)
        
        if current_seq:
            sequences.append(''.join(current_seq))
    
    # Convert to numerical array
    return sequences_to_array(sequences)


def load_stockholm(sto_path: Path) -> np.ndarray:
    """Load Stockholm format MSA"""
    sequences = {}
    
    with open(sto_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or line.startswith('//') or not line:
                continue
            
            parts = line.split()
            if len(parts) >= 2:
                seq_id = parts[0]
                seq = parts[1]
                
                if seq_id in sequences:
                    sequences[seq_id] += seq
                else:
                    sequences[seq_id] = seq
    
    return sequences_to_array(list(sequences.values()))


def load_fasta_msa(fasta_path: Path) -> np.ndarray:
    """Load FASTA format MSA"""
    sequences = []
    current_seq = []
    
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append(''.join(current_seq))
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_seq:
            sequences.append(''.join(current_seq))
    
    return sequences_to_array(sequences)


def sequences_to_array(sequences: list) -> np.ndarray:
    """
    Convert list of sequence strings to numerical array
    
    Encoding: A=0, C=1, G=2, U/T=3, -=4, other=5
    """
    if not sequences:
        return np.array([])
    
    # Find maximum length
    max_len = max(len(seq) for seq in sequences)
    
    # Encoding map
    encode = {'A': 0, 'C': 1, 'G': 2, 'U': 3, 'T': 3, '-': 4, '.': 4}
    
    # Convert to array
    n_seq = len(sequences)
    msa_array = np.full((n_seq, max_len), 5, dtype=np.int8)
    
    for i, seq in enumerate(sequences):
        for j, char in enumerate(seq):
            msa_array[i, j] = encode.get(char.upper(), 5)
    
    return msa_array


def compute_conservation_map(msa: np.ndarray) -> np.ndarray:
    """
    Compute per-position conservation scores
    
    Args:
        msa: MSA array (n_sequences x length)
        
    Returns:
        Conservation scores (length,)
    """
    n_seq, length = msa.shape
    conservation = np.zeros(length)
    
    for j in range(length):
        column = msa[:, j]
        
        # Remove gaps
        column_no_gaps = column[column < 4]
        
        if len(column_no_gaps) == 0:
            conservation[j] = 0.0
            continue
        
        # Compute entropy
        counts = np.bincount(column_no_gaps, minlength=4)
        freqs = counts / counts.sum()
        
        # Shannon entropy
        entropy = -np.sum(freqs * np.log2(freqs + 1e-10))
        
        # Convert to conservation (0 = variable, 1 = conserved)
        max_entropy = np.log2(4)  # 4 nucleotides
        conservation[j] = 1.0 - (entropy / max_entropy)
    
    return conservation


def compute_coevolution_map(msa: np.ndarray, method: str = "mi") -> np.ndarray:
    """
    Compute pairwise coevolution/contact map
    
    Args:
        msa: MSA array (n_sequences x length)
        method: Coevolution method (mi, apc, dca)
        
    Returns:
        Coevolution matrix (length x length)
    """
    n_seq, length = msa.shape
    
    if method == "mi":
        coevolution = compute_mutual_information(msa)
    elif method == "apc":
        mi = compute_mutual_information(msa)
        coevolution = apply_apc_correction(mi)
    elif method == "dca":
        coevolution = compute_dca(msa)
    else:
        raise ValueError(f"Unknown coevolution method: {method}")
    
    return coevolution


def compute_mutual_information(msa: np.ndarray) -> np.ndarray:
    """
    Compute mutual information between all position pairs
    
    MI(i,j) measures how much knowing position i tells us about position j
    High MI suggests coevolution (potential base pairing)
    """
    n_seq, length = msa.shape
    mi_matrix = np.zeros((length, length))
    
    for i in range(length):
        for j in range(i + 1, length):
            col_i = msa[:, i]
            col_j = msa[:, j]
            
            # Remove sequences with gaps in either position
            mask = (col_i < 4) & (col_j < 4)
            col_i_clean = col_i[mask]
            col_j_clean = col_j[mask]
            
            if len(col_i_clean) < 5:  # Need minimum sequences
                continue
            
            # Compute MI
            mi = calculate_mi(col_i_clean, col_j_clean)
            mi_matrix[i, j] = mi
            mi_matrix[j, i] = mi
    
    return mi_matrix


def calculate_mi(x: np.ndarray, y: np.ndarray) -> float:
    """Calculate mutual information between two columns"""
    # Joint distribution
    joint_counts = np.zeros((4, 4))
    for xi, yi in zip(x, y):
        if xi < 4 and yi < 4:
            joint_counts[xi, yi] += 1
    
    n = joint_counts.sum()
    if n == 0:
        return 0.0
    
    joint_probs = joint_counts / n
    
    # Marginal distributions
    px = joint_probs.sum(axis=1)
    py = joint_probs.sum(axis=0)
    
    # Mutual information
    mi = 0.0
    for i in range(4):
        for j in range(4):
            if joint_probs[i, j] > 0:
                mi += joint_probs[i, j] * np.log2(
                    joint_probs[i, j] / (px[i] * py[j] + 1e-10)
                )
    
    return mi


def apply_apc_correction(mi_matrix: np.ndarray) -> np.ndarray:
    """
    Apply Average Product Correction to remove background MI
    
    APC corrects for phylogenetic and entropic biases
    """
    length = mi_matrix.shape[0]
    
    # Compute average MI for each position
    avg_mi = mi_matrix.mean(axis=1)
    
    # Compute global average
    global_avg = mi_matrix.mean()
    
    # Apply APC
    apc_matrix = np.zeros_like(mi_matrix)
    for i in range(length):
        for j in range(length):
            if i != j:
                apc_matrix[i, j] = mi_matrix[i, j] - (avg_mi[i] * avg_mi[j] / global_avg)
    
    # Set negative values to 0
    apc_matrix = np.maximum(apc_matrix, 0)
    
    return apc_matrix


def compute_dca(msa: np.ndarray) -> np.ndarray:
    """
    Compute Direct Coupling Analysis (DCA)
    
    DCA removes indirect correlations to identify direct contacts
    This is a simplified implementation
    """
    # For full DCA, would use plmc or similar
    # Here we use APC-corrected MI as approximation
    mi = compute_mutual_information(msa)
    return apply_apc_correction(mi)


def compute_accessibility_map(msa: np.ndarray) -> np.ndarray:
    """
    Compute nucleotide accessibility/pairing map
    
    Args:
        msa: MSA array
        
    Returns:
        Accessibility scores (length,)
    """
    n_seq, length = msa.shape
    
    # Estimate accessibility from gap frequency
    # Positions with many gaps are likely in loops (accessible)
    # Positions with few gaps are likely in stems (paired)
    
    gap_freq = (msa == 4).sum(axis=0) / n_seq
    
    # High gap frequency = high accessibility
    accessibility = gap_freq
    
    # Also consider conservation
    # Highly conserved positions in stems should be less accessible
    conservation = compute_conservation_map(msa)
    
    # Combine: conserved + few gaps = paired (low accessibility)
    #          variable or many gaps = unpaired (high accessibility)
    accessibility = 0.7 * gap_freq + 0.3 * (1 - conservation)
    
    return accessibility


def structure_to_features(structure_dict: Dict) -> Dict[str, np.ndarray]:
    """
    Convert single-sequence structure prediction to 2D features
    
    Args:
        structure_dict: Dictionary from single_sequence.predict_structure_spotRNA
        
    Returns:
        Dictionary with feature arrays
    """
    sequence = structure_dict['sequence']
    structure = structure_dict['structure']
    length = len(sequence)
    
    # Conservation: use GC content as proxy
    conservation = np.zeros(length)
    for i, nt in enumerate(sequence):
        if nt in 'GCgc':
            conservation[i] = 0.8
        else:
            conservation[i] = 0.4
    
    # Coevolution: use pairing probability if available
    if 'pairing_prob' in structure_dict:
        coevolution = structure_dict['pairing_prob']
    else:
        # Convert structure to contact map
        from .single_sequence import structure_to_contact_map
        coevolution = structure_to_contact_map(structure)
    
    # Accessibility
    if 'accessibility' in structure_dict:
        accessibility = structure_dict['accessibility']
    else:
        # Compute from structure
        from .single_sequence import compute_accessibility
        accessibility = compute_accessibility(sequence, structure)
    
    features = {
        'conservation': conservation,
        'coevolution': coevolution,
        'accessibility': accessibility,
        'msa_depth': 0  # Single sequence
    }
    
    return features


def combine_domain_features(domain_features: list, domain_positions: list) -> Dict[str, np.ndarray]:
    """
    Combine features from multiple domains into full-length features
    
    Args:
        domain_features: List of feature dictionaries for each domain
        domain_positions: List of (start, end) tuples for each domain
        
    Returns:
        Combined feature dictionary
    """
    # Determine full length
    max_end = max(end for start, end in domain_positions)
    
    # Initialize combined features
    conservation = np.zeros(max_end)
    accessibility = np.zeros(max_end)
    coevolution = np.zeros((max_end, max_end))
    
    # Combine domain features with overlap handling
    counts = np.zeros(max_end)
    
    for features, (start, end) in zip(domain_features, domain_positions):
        domain_len = end - start
        
        # Add conservation and accessibility
        conservation[start:end] += features['conservation'][:domain_len]
        accessibility[start:end] += features['accessibility'][:domain_len]
        counts[start:end] += 1
        
        # Add coevolution
        coev = features['coevolution']
        coevolution[start:end, start:end] += coev[:domain_len, :domain_len]
    
    # Average overlapping regions
    mask = counts > 0
    conservation[mask] /= counts[mask]
    accessibility[mask] /= counts[mask]
    
    # Normalize coevolution
    for i in range(max_end):
        for j in range(max_end):
            if counts[i] > 0 and counts[j] > 0:
                coevolution[i, j] /= min(counts[i], counts[j])
    
    return {
        'conservation': conservation,
        'coevolution': coevolution,
        'accessibility': accessibility,
        'msa_depth': sum(f['msa_depth'] for f in domain_features) / len(domain_features)
    }
