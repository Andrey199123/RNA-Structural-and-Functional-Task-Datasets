#!/usr/bin/env python3
"""
Splits long RNAs into smaller chunks and processes each one separately.
"""

import logging
from pathlib import Path
from typing import List
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class RNADomain:
    id: str
    sequence: str
    length: int
    start: int  # Position in parent RNA
    end: int
    parent_id: str
    description: str = None


def split_into_domains(rna, window_size: int = 300, overlap: int = 50) -> List[RNADomain]:
    domains = []
    step = window_size - overlap
    
    for i, start in enumerate(range(0, rna.length, step)):
        end = min(start + window_size, rna.length)
        
        if end - start < 100:
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
    window_size = 200
    step = 50
    
    energies = []
    positions = []
    
    for start in range(0, rna.length - window_size, step):
        end = start + window_size
        segment = rna.sequence[start:end]
        energy = compute_folding_energy(segment)
        energies.append(energy)
        positions.append((start, end))
    
    domain_boundaries = find_stable_regions(energies, positions)
    
    logger.info(f"Predicted {len(domain_boundaries)} domain boundaries for {rna.id}")
    
    return domain_boundaries


def compute_folding_energy(sequence: str) -> float:
    try:
        import RNA
        structure, mfe = RNA.fold(sequence)
        return mfe
    except ImportError:
        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
        return -gc_content * 10


def find_stable_regions(energies: List[float], positions: List[tuple]) -> List[tuple]:
    if not energies:
        return []
    
    import numpy as np
    energies = np.array(energies)
    
    from scipy.ndimage import gaussian_filter1d
    smoothed = gaussian_filter1d(energies, sigma=2)
    
    from scipy.signal import find_peaks
    peaks, _ = find_peaks(-smoothed, distance=3)
    
    domains = []
    for i in range(len(peaks) - 1):
        start_idx = peaks[i]
        end_idx = peaks[i + 1]
        
        start_pos = positions[start_idx][0]
        end_pos = positions[end_idx][1]
        
        domains.append((start_pos, end_pos))
    
    return domains


def merge_overlapping_domains(domains: List[RNADomain], min_overlap: int = 30) -> List[RNADomain]:
    if not domains:
        return []
    
    sorted_domains = sorted(domains, key=lambda d: d.start)
    
    merged = [sorted_domains[0]]
    
    for current in sorted_domains[1:]:
        previous = merged[-1]
        
        overlap = previous.end - current.start
        
        if overlap >= min_overlap:
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
    logger.info(f"Processing domain {domain.id} ({domain.length}nt)")
    
    return {
        'domain_id': domain.id,
        'start': domain.start,
        'end': domain.end,
        'length': domain.length,
        'strategy': strategy
    }
