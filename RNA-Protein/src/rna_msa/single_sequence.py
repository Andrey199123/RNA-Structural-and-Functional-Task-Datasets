#!/usr/bin/env python3
"""
Predicts structure for RNAs that have no homologs, so no MSA can be built.
"""

import subprocess
import logging
from pathlib import Path
from typing import Dict, Tuple
import numpy as np

logger = logging.getLogger(__name__)


def structure_to_contact_map(structure: str) -> np.ndarray:
    n = len(structure)
    contact_map = np.zeros((n, n))
    
    stack = []
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                j = stack.pop()
                contact_map[i, j] = 1
                contact_map[j, i] = 1
    
    return contact_map


def predict_structure_spotRNA(rna, output_dir: Path) -> Dict:
    logger.info(f"Predicting structure for {rna.id} using single-sequence method")
    
    try:
        structure = run_spotRNA(rna, output_dir)
        if structure:
            return structure
    except Exception as e:
        logger.warning(f"SPOT-RNA failed for {rna.id}: {e}")
    
    try:
        structure = run_RNAfold(rna, output_dir)
        return structure
    except Exception as e:
        logger.error(f"RNAfold also failed for {rna.id}: {e}")
        raise


def run_spotRNA(rna, output_dir: Path) -> Dict:
    input_file = output_dir / f"{rna.id}_input.fasta"
    with open(input_file, 'w') as f:
        f.write(f">{rna.id}\n{rna.sequence}\n")
    
    output_prefix = output_dir / f"{rna.id}_spotRNA"
    
    logger.warning("SPOT-RNA not available, using RNAfold")
    return run_RNAfold(rna, output_dir)


def run_RNAfold(rna, output_dir: Path) -> Dict:
    try:
        import RNA
        
        structure, mfe = RNA.fold(rna.sequence)
        
        fc = RNA.fold_compound(rna.sequence)
        fc.pf()
        
        n = rna.length
        prob_matrix = np.zeros((n, n))
        
        try:
            for i in range(1, n + 1):
                for j in range(i + 1, n + 1):
                    prob = fc.pr(i, j)
                    if prob > 0.01:
                        prob_matrix[i-1, j-1] = prob
                        prob_matrix[j-1, i-1] = prob
        except (AttributeError, TypeError):
            try:
                bpp = fc.bpp()
                for i in range(1, n + 1):
                    for j in range(i + 1, n + 1):
                        if i < len(bpp) and j < len(bpp[i]):
                            prob = bpp[i][j]
                            if prob > 0.01:
                                prob_matrix[i-1, j-1] = prob
                                prob_matrix[j-1, i-1] = prob
            except (AttributeError, TypeError):
                logger.warning("Could not get base pairing probabilities, using structure only")
                prob_matrix = structure_to_contact_map(structure)
        
        accessibility = compute_accessibility(rna.sequence, structure)
        
        result = {
            'sequence': rna.sequence,
            'structure': structure,
            'mfe': mfe,
            'pairing_prob': prob_matrix,
            'accessibility': accessibility,
            'method': 'RNAfold'
        }
        
        output_file = output_dir / f"{rna.id}_structure.npz"
        np.savez_compressed(
            output_file,
            structure=structure,
            mfe=mfe,
            pairing_prob=prob_matrix,
            accessibility=accessibility
        )
        
        logger.info(f"Predicted structure for {rna.id}: MFE={mfe:.2f}")
        
        return result
        
    except ImportError:
        logger.error("ViennaRNA Python bindings not available")
        return run_RNAfold_cli(rna, output_dir)


def run_RNAfold_cli(rna, output_dir: Path) -> Dict:
    input_file = output_dir / f"{rna.id}_input.fasta"
    with open(input_file, 'w') as f:
        f.write(f">{rna.id}\n{rna.sequence}\n")
    
    try:
        cmd = ["RNAfold", "--noPS", "-i", str(input_file)]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, cwd=output_dir)
        
        lines = result.stdout.strip().split('\n')
        structure = None
        mfe = None
        
        for line in lines:
            if line and not line.startswith('>'):
                parts = line.split()
                if len(parts) >= 2:
                    structure = parts[0]
                    energy_str = parts[1].strip('()')
                    try:
                        mfe = float(energy_str)
                    except ValueError:
                        pass
                break
        
        if not structure:
            raise ValueError("Could not parse RNAfold output")
        
        accessibility = compute_accessibility(rna.sequence, structure)
        
        result = {
            'sequence': rna.sequence,
            'structure': structure,
            'mfe': mfe or 0.0,
            'accessibility': accessibility,
            'method': 'RNAfold_cli'
        }
        
        logger.info(f"Predicted structure for {rna.id} using RNAfold CLI")
        
        return result
        
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.error(f"RNAfold CLI failed: {e}")
        return {
            'sequence': rna.sequence,
            'structure': '.' * rna.length,
            'mfe': 0.0,
            'accessibility': np.ones(rna.length),
            'method': 'fallback'
        }


def compute_accessibility(sequence: str, structure: str) -> np.ndarray:
    n = len(sequence)
    accessibility = np.zeros(n)
    
    for i, char in enumerate(structure):
        if char == '.':
            accessibility[i] = 1.0
        elif char in '()':
            accessibility[i] = 0.3
        else:
            accessibility[i] = 0.5
    
    return accessibility


def compute_conservation_single(sequence: str) -> np.ndarray:
    n = len(sequence)
    conservation = np.zeros(n)
    
    window = 5
    
    for i in range(n):
        start = max(0, i - window)
        end = min(n, i + window + 1)
        local_seq = sequence[start:end]
        
        gc_count = local_seq.count('G') + local_seq.count('C')
        gc_content = gc_count / len(local_seq)
        
        conservation[i] = gc_content
    
    return conservation


def predict_binding_sites(rna, structure_dict: Dict) -> np.ndarray:
    structure = structure_dict['structure']
    accessibility = structure_dict.get('accessibility', np.ones(rna.length))
    
    n = rna.length
    binding_score = np.zeros(n)
    
    for i, char in enumerate(structure):
        if char == '.':
            binding_score[i] = 0.8
        else:
            binding_score[i] = 0.2
    
    binding_score *= accessibility
    
    motifs = identify_structural_motifs(structure)
    for motif_start, motif_end in motifs:
        binding_score[motif_start:motif_end] *= 1.5
    
    binding_score = np.clip(binding_score, 0, 1)
    
    return binding_score


def identify_structural_motifs(structure: str) -> list:
    motifs = []
    
    i = 0
    while i < len(structure):
        if structure[i] == '(':
            depth = 1
            j = i + 1
            while j < len(structure) and depth > 0:
                if structure[j] == '(':
                    depth += 1
                elif structure[j] == ')':
                    depth -= 1
                j += 1
            
            loop_size = j - i
            if 3 <= loop_size <= 10:
                motifs.append((i, j))
            
            i = j
        else:
            i += 1
    
    return motifs
