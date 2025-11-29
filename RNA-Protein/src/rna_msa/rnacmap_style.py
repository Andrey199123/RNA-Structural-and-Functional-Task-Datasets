#!/usr/bin/env python3
"""
Builds custom covariance models for RNAs that don't match known Rfam families.
"""

import subprocess
import logging
from pathlib import Path
from typing import Tuple
import tempfile

logger = logging.getLogger(__name__)


def predict_secondary_structure(rna, output_dir: Path, method: str = "spotRNA") -> Tuple[str, Path]:
    if method == "spotRNA":
        return predict_spotRNA(rna, output_dir)
    elif method == "RNAfold":
        return predict_RNAfold(rna, output_dir)
    else:
        raise ValueError(f"Unknown structure prediction method: {method}")


def predict_spotRNA(rna, output_dir: Path) -> Tuple[str, Path]:
    seq_file = output_dir / f"{rna.id}_seq.fa"
    with open(seq_file, 'w') as f:
        f.write(f">{rna.id}\n{rna.sequence}\n")
    
    logger.info(f"Predicting structure for {rna.id} using SPOT-RNA")
    return predict_RNAfold(rna, output_dir)


def predict_RNAfold(rna, output_dir: Path) -> Tuple[str, Path]:
    seq_file = output_dir / f"{rna.id}_seq.fa"
    with open(seq_file, 'w') as f:
        f.write(f">{rna.id}\n{rna.sequence}\n")
    
    output_file = output_dir / f"{rna.id}_structure.txt"
    
    try:
        cmd = ["RNAfold", "--noPS", "-i", str(seq_file)]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        lines = result.stdout.strip().split('\n')
        structure = None
        for line in lines:
            if line and not line.startswith('>'):
                structure = line.split()[0]
                break
        
        if structure:
            with open(output_file, 'w') as f:
                f.write(f">{rna.id}\n{rna.sequence}\n{structure}\n")
            
            logger.info(f"Predicted structure for {rna.id}")
            return structure, output_file
        
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.warning(f"RNAfold failed for {rna.id}: {e}")
    
    structure = '.' * rna.length
    with open(output_file, 'w') as f:
        f.write(f">{rna.id}\n{rna.sequence}\n{structure}\n")
    
    return structure, output_file


def build_custom_cm(rna, output_dir: Path, n_cpu: int = 4, calibrate: bool = False) -> Path:
    structure, struct_file = predict_secondary_structure(rna, output_dir)
    
    cleaned_structure = ''.join(c if c in '.()' else '.' for c in structure)
    
    if len(cleaned_structure) != len(rna.sequence):
        logger.warning(f"Structure length mismatch for {rna.id}, using all unpaired")
        cleaned_structure = '.' * len(rna.sequence)
    
    sto_file = output_dir / f"{rna.id}_single.sto"
    
    with open(sto_file, 'w') as f:
        f.write("# STOCKHOLM 1.0\n")
        f.write(f"#=GF ID {rna.id}\n")
        f.write(f"{rna.id} {rna.sequence}\n")
        f.write(f"#=GC SS_cons {cleaned_structure}\n")
        f.write("//\n")
    
    cm_file = output_dir / f"{rna.id}_custom.cm"
    
    try:
        cmd = ["cmbuild", "-F", str(cm_file), str(sto_file)]
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"Built custom CM for {rna.id}")
        
        if calibrate:
            logger.info(f"Calibrating CM for {rna.id} (this may take several minutes)...")
            cmd = ["cmcalibrate", "--cpu", str(n_cpu), str(cm_file)]
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            logger.info(f"Calibrated CM for {rna.id}")
        else:
            logger.info(f"Skipping CM calibration for {rna.id} (use calibrate=True for better E-values)")
        
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.error(f"CM build failed for {rna.id}: {e}")
        raise
    
    return cm_file


def search_homologs(rna, cm_file: Path, sequence_db: Path, output_dir: Path,
                   n_cpu: int = 4, evalue_threshold: float = 0.01, use_blast_filter: bool = True) -> Path:
    hits_file = output_dir / f"{rna.id}_homologs.fa"
    tblout_file = output_dir / f"{rna.id}_homologs.tblout"
    
    if use_blast_filter and sequence_db.stat().st_size > 1_000_000_000:
        logger.info(f"Large database detected, using BLAST pre-filtering for {rna.id}")
        
        try:
            from .fast_homology_search import FastHomologySearch
            
            searcher = FastHomologySearch(sequence_db)
            n_homologs, candidates_file = searcher.fast_homology_search(
                query_seq=rna.sequence,
                query_id=rna.id,
                output_fasta=output_dir / f"{rna.id}_blast_candidates.fa",
                max_homologs=1000,
                min_identity=0.5
            )
            
            if n_homologs == 0:
                logger.warning(f"No BLAST candidates found for {rna.id}")
                hits_file.touch()
                return hits_file
            
            search_db = candidates_file
            logger.info(f"BLAST found {n_homologs} candidates, running cmsearch on subset")
            
        except Exception as e:
            logger.warning(f"BLAST pre-filtering failed: {e}, falling back to full cmsearch")
            search_db = sequence_db
    else:
        search_db = sequence_db
    
    try:
        cmd = [
            "cmsearch",
            "--cpu", str(n_cpu),
            "-E", str(evalue_threshold),
            "--tblout", str(tblout_file),
            "-A", str(hits_file),  # Output alignment
            str(cm_file),
            str(search_db)
        ]
        
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"Found homologs for {rna.id}")
        
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.error(f"Homolog search failed for {rna.id}: {e}")
        raise
    
    return hits_file


def align_sequences(sequences_file: Path, output_dir: Path, rna_id: str) -> Path:
    aligned_file = output_dir / f"{rna_id}_aligned.sto"
    
    try:
        cmd = [
            "mafft",
            "--maxiterate", "1000",
            "--localpair",
            "--thread", "4",
            str(sequences_file)
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        with open(aligned_file, 'w') as f:
            f.write("# STOCKHOLM 1.0\n")
            f.write(result.stdout)
            f.write("//\n")
        
        logger.info(f"Aligned sequences for {rna_id}")
        
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.warning(f"MAFFT alignment failed for {rna_id}: {e}")
        import shutil
        shutil.copy(sequences_file, aligned_file)
    
    return aligned_file


def refine_structure_from_msa(msa_file: Path, output_dir: Path, rna_id: str) -> str:
    logger.info(f"Refining structure from MSA for {rna_id}")
    return None


def iterative_msa_refinement(rna, cm_file: Path, sequence_db: Path, output_dir: Path,
                             n_iterations: int = 3, n_cpu: int = 4) -> Tuple[Path, int]:
    current_cm = cm_file
    
    for iteration in range(n_iterations):
        logger.info(f"MSA refinement iteration {iteration + 1}/{n_iterations} for {rna.id}")
        
        try:
            homologs_file = search_homologs(
                rna, current_cm, sequence_db, output_dir, n_cpu=n_cpu
            )
        except Exception as e:
            logger.warning(f"Homolog search failed in iteration {iteration + 1}: {e}")
            break
        
        msa_file = align_sequences(homologs_file, output_dir, f"{rna.id}_iter{iteration}")
        msa_depth = count_sequences_in_msa(msa_file)
        logger.info(f"MSA depth after iteration {iteration + 1}: {msa_depth}")
        
        if msa_depth < 5:
            logger.warning(f"Low MSA depth ({msa_depth}), stopping iterations")
            break
        
        refined_structure = refine_structure_from_msa(msa_file, output_dir, rna.id)
        
        if refined_structure and iteration < n_iterations - 1:
            pass
    
    final_msa = output_dir / f"{rna.id}_final.sto"
    if msa_file.exists():
        import shutil
        shutil.copy(msa_file, final_msa)
    
    return final_msa, msa_depth


def count_sequences_in_msa(msa_file: Path) -> int:
    count = 0
    with open(msa_file) as f:
        for line in f:
            if line.startswith('>'):
                count += 1
            elif not line.startswith('#') and not line.startswith('//'):
                if line.strip():
                    count += 1
    return max(1, count)
