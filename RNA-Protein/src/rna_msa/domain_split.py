#!/usr/bin/env python3
"""
Domain Split Strategy
For long RNAs (>1000 nt)
"""

import logging
from pathlib import Path
from typing import List
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class RNADomain:
    """RNA domain/segment"""
    id: str
    sequence: str
    length: int
    start: int  # Position in parent RNA
    end: int
    parent_id: str
    description: str = None


def split_into_domains(rna, window_size: int = 300, overlap: int = 50) -> List[RNADomain]:
    """
    Split long RNA into overlapping domains
    
    Args:
        rna: RNASequence object
        window_size: Size of each domain window
        overlap: Overlap between adjacent windows
        
    Returns:
        List of RNADomain objects
    """
    domains = []
    step = window_size - overlap
    
    for i, start in enumerate(range(0, rna.length, step)):
        end = min(start + window_size, rna.length)
        
        if end - start < 100:  # Skip very short terminal segments
            break
        
        domain_seq = rna.sequence[start:end]
        domain = RNADomain(
            id=f"{rna.id}_domain{i}",
            sequence=domain_seq,
            length=len(domain_seq),
            start=start,
            end=end,
            parent_id=rna.id,
            description=f"Domain {i} of {rna.id} ({start}-{end})"
        )
        domains.append(domain)
    
    logger.info(f"Split {rna.id} ({rna.length}nt) into {len(domains)} domains")
    
    return domains


def predict_domain_boundaries(rna, output_dir: Path) -> List[tuple]:
    """
    Predict structural domain boundaries using folding energy
    
    This is more sophisticated than simple windowing - it tries to
    identify natural structural domains based on folding patterns.
    
    Args:
        rna: RNASequence object
        output_dir: Output directory
        
    Returns:
        List of (start, end) tuples for domain boundaries
    """
    # Use sliding window to compute local folding stability
    window_size = 200
    step = 50
    
    energies = []
    positions = []
    
    for start in range(0, rna.length - window_size, step):
        end = start + window_size
        segment = rna.sequence[start:end]
        
        # Compute folding energy (placeholder - would use ViennaRNA)
        energy = compute_folding_energy(segment)
        energies.append(energy)
        positions.append((start, end))
    
    # Find local minima (stable domains)
    domain_boundaries = find_stable_regions(energies, positions)
    
    logger.info(f"Predicted {len(domain_boundaries)} domain boundaries for {rna.id}")
    
    return domain_boundaries


def compute_folding_energy(sequence: str) -> float:
    """
    Compute RNA folding free energy
    
    Args:
        sequence: RNA sequence
        
    Returns:
        Folding energy (kcal/mol)
    """
    try:
        import RNA  # ViennaRNA Python bindings
        structure, mfe = RNA.fold(sequence)
        return mfe
    except ImportError:
        # Fallback: simple GC content heuristic
        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
        return -gc_content * 10  # Rough approximation


def find_stable_regions(energies: List[float], positions: List[tuple]) -> List[tuple]:
    """
    Find stable folding regions from energy profile
    
    Args:
        energies: List of folding energies
        positions: List of (start, end) positions
        
    Returns:
        List of domain boundaries
    """
    if not energies:
        return []
    
    # Simple approach: find local minima
    import numpy as np
    energies = np.array(energies)
    
    # Smooth energy profile
    from scipy.ndimage import gaussian_filter1d
    smoothed = gaussian_filter1d(energies, sigma=2)
    
    # Find local minima
    from scipy.signal import find_peaks
    peaks, _ = find_peaks(-smoothed, distance=3)
    
    # Convert peaks to domain boundaries
    domains = []
    for i in range(len(peaks) - 1):
        start_idx = peaks[i]
        end_idx = peaks[i + 1]
        
        start_pos = positions[start_idx][0]
        end_pos = positions[end_idx][1]
        
        domains.append((start_pos, end_pos))
    
    return domains


def merge_overlapping_domains(domains: List[RNADomain], min_overlap: int = 30) -> List[RNADomain]:
    """
    Merge domains with significant overlap
    
    Args:
        domains: List of RNADomain objects
        min_overlap: Minimum overlap to trigger merge
        
    Returns:
        List of merged domains
    """
    if not domains:
        return []
    
    # Sort by start position
    sorted_domains = sorted(domains, key=lambda d: d.start)
    
    merged = [sorted_domains[0]]
    
    for current in sorted_domains[1:]:
        previous = merged[-1]
        
        overlap = previous.end - current.start
        
        if overlap >= min_overlap:
            # Merge domains
            merged_seq = previous.sequence + current.sequence[overlap:]
            merged_domain = RNADomain(
                id=f"{previous.parent_id}_merged_{len(merged)}",
                sequence=merged_seq,
                length=len(merged_seq),
                start=previous.start,
                end=current.end,
                parent_id=previous.parent_id,
                description=f"Merged domain {previous.start}-{current.end}"
            )
            merged[-1] = merged_domain
        else:
            merged.append(current)
    
    return merged


def process_domain(domain: RNADomain, strategy: str = "rnacmap") -> dict:
    """
    Process a single domain
    
    Args:
        domain: RNADomain object
        strategy: Processing strategy
        
    Returns:
        Dictionary with domain results
    """
    # This would call the appropriate pipeline strategy
    # Placeholder for now
    logger.info(f"Processing domain {domain.id} ({domain.length}nt)")
    
    return {
        'domain_id': domain.id,
        'start': domain.start,
        'end': domain.end,
        'length': domain.length,
        'strategy': strategy
    }
