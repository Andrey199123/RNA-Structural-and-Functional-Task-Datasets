#!/usr/bin/env python3
"""
RNAcmap-Style Pipeline
For novel RNAs without Rfam matches
"""

import subprocess
import logging
from pathlib import Path
from typing import Tuple
import tempfile

logger = logging.getLogger(__name__)


def predict_secondary_structure(rna, output_dir: Path, method: str = "spotRNA") -> Tuple[str, Path]:
    """
    Predict RNA secondary structure
    
    Args:
        rna: RNASequence object
        output_dir: Output directory
        method: Prediction method (spotRNA, RNAfold, etc.)
        
    Returns:
        Tuple of (dot-bracket structure, structure file path)
    """
    if method == "spotRNA":
        return predict_spotRNA(rna, output_dir)
    elif method == "RNAfold":
        return predict_RNAfold(rna, output_dir)
    else:
        raise ValueError(f"Unknown structure prediction method: {method}")


def predict_spotRNA(rna, output_dir: Path) -> Tuple[str, Path]:
    """Predict structure using SPOT-RNA"""
    # Write sequence
    seq_file = output_dir / f"{rna.id}_seq.fa"
    with open(seq_file, 'w') as f:
        f.write(f">{rna.id}\n{rna.sequence}\n")
    
    # Run SPOT-RNA (placeholder - actual implementation depends on installation)
    # In practice: python SPOT-RNA.py --input seq_file --output output_dir
    
    logger.info(f"Predicting structure for {rna.id} using SPOT-RNA")
    
    # For now, use RNAfold as fallback
    return predict_RNAfold(rna, output_dir)


def predict_RNAfold(rna, output_dir: Path) -> Tuple[str, Path]:
    """Predict structure using ViennaRNA RNAfold"""
    seq_file = output_dir / f"{rna.id}_seq.fa"
    with open(seq_file, 'w') as f:
        f.write(f">{rna.id}\n{rna.sequence}\n")
    
    output_file = output_dir / f"{rna.id}_structure.txt"
    
    try:
        cmd = ["RNAfold", "--noPS", "-i", str(seq_file)]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Parse output
        lines = result.stdout.strip().split('\n')
        structure = None
        for line in lines:
            if line and not line.startswith('>'):
                # Structure line contains dot-bracket and energy
                structure = line.split()[0]
                break
        
        if structure:
            with open(output_file, 'w') as f:
                f.write(f">{rna.id}\n{rna.sequence}\n{structure}\n")
            
            logger.info(f"Predicted structure for {rna.id}")
            return structure, output_file
        
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.warning(f"RNAfold failed for {rna.id}: {e}")
    
    # Fallback: unpaired structure
    structure = '.' * rna.length
    with open(output_file, 'w') as f:
        f.write(f">{rna.id}\n{rna.sequence}\n{structure}\n")
    
    return structure, output_file


def build_custom_cm(rna, output_dir: Path, n_cpu: int = 4) -> Path:
    """
    Build custom covariance model from predicted structure
    
    Args:
        rna: RNASequence object
        output_dir: Output directory
        n_cpu: Number of CPUs
        
    Returns:
        Path to covariance model file
    """
    # Predict structure
    structure, struct_file = predict_secondary_structure(rna, output_dir)
    
    # Clean structure - remove any non-structure characters
    # Valid structure chars: . ( ) < > [ ] { } (and sometimes letters for pseudoknots)
    # For basic structure, keep only . ( )
    cleaned_structure = ''.join(c if c in '.()' else '.' for c in structure)
    
    # Ensure structure is same length as sequence
    if len(cleaned_structure) != len(rna.sequence):
        logger.warning(f"Structure length mismatch for {rna.id}, using all unpaired")
        cleaned_structure = '.' * len(rna.sequence)
    
    # Create Stockholm format alignment (single sequence)
    sto_file = output_dir / f"{rna.id}_single.sto"
    
    with open(sto_file, 'w') as f:
        f.write("# STOCKHOLM 1.0\n")
        f.write(f"#=GF ID {rna.id}\n")
        # Ensure proper spacing - sequence ID should be exactly as is
        f.write(f"{rna.id} {rna.sequence}\n")
        f.write(f"#=GC SS_cons {cleaned_structure}\n")
        f.write("//\n")
    
    # Build CM using cmbuild
    cm_file = output_dir / f"{rna.id}_custom.cm"
    
    try:
        cmd = ["cmbuild", "-F", str(cm_file), str(sto_file)]
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"Built custom CM for {rna.id}")
        
        # Calibrate CM
        cmd = ["cmcalibrate", "--cpu", str(n_cpu), str(cm_file)]
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"Calibrated CM for {rna.id}")
        
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.error(f"CM build failed for {rna.id}: {e}")
        raise
    
    return cm_file


def search_homologs(rna, cm_file: Path, sequence_db: Path, output_dir: Path,
                   n_cpu: int = 4, evalue_threshold: float = 0.01) -> Path:
    """
    Search for homologs using custom CM
    
    Args:
        rna: RNASequence object
        cm_file: Covariance model file
        sequence_db: Sequence database to search
        output_dir: Output directory
        n_cpu: Number of CPUs
        evalue_threshold: E-value threshold
        
    Returns:
        Path to homolog sequences file
    """
    hits_file = output_dir / f"{rna.id}_homologs.fa"
    tblout_file = output_dir / f"{rna.id}_homologs.tblout"
    
    try:
        # Run cmsearch
        cmd = [
            "cmsearch",
            "--cpu", str(n_cpu),
            "-E", str(evalue_threshold),
            "--tblout", str(tblout_file),
            "-A", str(hits_file),  # Output alignment
            str(cm_file),
            str(sequence_db)
        ]
        
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"Found homologs for {rna.id}")
        
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.error(f"Homolog search failed for {rna.id}: {e}")
        raise
    
    return hits_file


def align_sequences(sequences_file: Path, output_dir: Path, rna_id: str) -> Path:
    """
    Align sequences using structure-aware aligner
    
    Args:
        sequences_file: FASTA file with sequences to align
        output_dir: Output directory
        rna_id: RNA identifier
        
    Returns:
        Path to alignment file (Stockholm format)
    """
    aligned_file = output_dir / f"{rna_id}_aligned.sto"
    
    try:
        # Use MAFFT with RNA-specific options
        cmd = [
            "mafft",
            "--maxiterate", "1000",
            "--localpair",
            "--thread", "4",
            str(sequences_file)
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Convert FASTA output to Stockholm
        with open(aligned_file, 'w') as f:
            f.write("# STOCKHOLM 1.0\n")
            f.write(result.stdout)
            f.write("//\n")
        
        logger.info(f"Aligned sequences for {rna_id}")
        
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.warning(f"MAFFT alignment failed for {rna_id}: {e}")
        # Copy input as-is
        import shutil
        shutil.copy(sequences_file, aligned_file)
    
    return aligned_file


def refine_structure_from_msa(msa_file: Path, output_dir: Path, rna_id: str) -> str:
    """
    Refine secondary structure based on covariation in MSA
    
    Args:
        msa_file: MSA file (Stockholm format)
        output_dir: Output directory
        rna_id: RNA identifier
        
    Returns:
        Refined dot-bracket structure
    """
    # Use R-scape or similar to detect covarying pairs
    # For now, return placeholder
    logger.info(f"Refining structure from MSA for {rna_id}")
    
    # TODO: Implement R-scape integration
    # cmd = ["R-scape", "--outdir", str(output_dir), str(msa_file)]
    
    return None  # Will use original structure


def iterative_msa_refinement(rna, cm_file: Path, sequence_db: Path, output_dir: Path,
                             n_iterations: int = 3, n_cpu: int = 4) -> Tuple[Path, int]:
    """
    Iteratively refine MSA and structure
    
    Args:
        rna: RNASequence object
        cm_file: Initial covariance model
        sequence_db: Sequence database
        output_dir: Output directory
        n_iterations: Number of refinement iterations
        n_cpu: Number of CPUs
        
    Returns:
        Tuple of (final MSA path, MSA depth)
    """
    current_cm = cm_file
    
    for iteration in range(n_iterations):
        logger.info(f"MSA refinement iteration {iteration + 1}/{n_iterations} for {rna.id}")
        
        # Search for homologs
        try:
            homologs_file = search_homologs(
                rna, current_cm, sequence_db, output_dir, n_cpu=n_cpu
            )
        except Exception as e:
            logger.warning(f"Homolog search failed in iteration {iteration + 1}: {e}")
            break
        
        # Align sequences
        msa_file = align_sequences(homologs_file, output_dir, f"{rna.id}_iter{iteration}")
        
        # Count sequences in MSA
        msa_depth = count_sequences_in_msa(msa_file)
        logger.info(f"MSA depth after iteration {iteration + 1}: {msa_depth}")
        
        if msa_depth < 5:
            logger.warning(f"Low MSA depth ({msa_depth}), stopping iterations")
            break
        
        # Refine structure
        refined_structure = refine_structure_from_msa(msa_file, output_dir, rna.id)
        
        # Rebuild CM if structure was refined
        if refined_structure and iteration < n_iterations - 1:
            # Update CM for next iteration
            pass
    
    # Return final MSA
    final_msa = output_dir / f"{rna.id}_final.sto"
    if msa_file.exists():
        import shutil
        shutil.copy(msa_file, final_msa)
    
    return final_msa, msa_depth


def count_sequences_in_msa(msa_file: Path) -> int:
    """Count number of sequences in MSA"""
    count = 0
    with open(msa_file) as f:
        for line in f:
            if line.startswith('>'):
                count += 1
            elif not line.startswith('#') and not line.startswith('//'):
                if line.strip():
                    count += 1
    return max(1, count)
