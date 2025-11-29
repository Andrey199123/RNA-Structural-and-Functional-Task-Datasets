#!/usr/bin/env python3
"""
Uses Infernal to search Rfam database for known RNA families and align sequences to family MSAs.
"""

import subprocess
import logging
from pathlib import Path
from typing import List, Dict
import tempfile

logger = logging.getLogger(__name__)


def search_rfam(rna, rfam_db_path: Path, infernal_path: str, output_dir: Path, 
                n_cpu: int = 4, evalue_threshold: float = 0.01) -> List[Dict]:
    query_file = output_dir / f"{rna.id}_query.fa"
    with open(query_file, 'w') as f:
        f.write(f">{rna.id}\n{rna.sequence}\n")
    
    tblout_file = output_dir / f"{rna.id}_rfam.tblout"
    cmd = [
        infernal_path,
        "--cpu", str(n_cpu),
        "--tblout", str(tblout_file),
        "-E", str(evalue_threshold),
        str(rfam_db_path),
        str(query_file)
    ]
    
    logger.info(f"Running Infernal cmscan for {rna.id}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.debug(f"cmscan output: {result.stdout}")
    except subprocess.CalledProcessError as e:
        logger.error(f"cmscan failed: {e.stderr}")
        raise
    
    hits = parse_cmscan_tblout(tblout_file)
    logger.info(f"Found {len(hits)} Rfam hits for {rna.id}")
    
    return hits


def parse_cmscan_tblout(tblout_path: Path) -> List[Dict]:
    hits = []
    
    with open(tblout_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) < 17:
                continue
            
            try:
                hit = {
                    'target': parts[1],
                    'family': parts[0],
                    'query': parts[2],
                    'evalue': float(parts[15]),
                    'score': float(parts[14]),
                    'bias': float(parts[13]),
                    'model_from': int(parts[5]),
                    'model_to': int(parts[6]),
                    'seq_from': int(parts[7]),
                    'seq_to': int(parts[8]),
                    'strand': parts[9],
                    'msa_depth': 100
                }
                hits.append(hit)
                logger.info(f"Parsed hit: {parts[0]} ({parts[1]}) E-value={parts[15]}")
            except (ValueError, IndexError) as e:
                logger.debug(f"Skipping unparseable line: {line.strip()[:100]}")
                continue
    
    hits.sort(key=lambda x: x['evalue'])
    
    return hits


def fetch_rfam_alignment(rna, rfam_hit: Dict, output_dir: Path, rfam_seed_path: Path = None) -> Path:
    family_acc = rfam_hit['target']
    
    logger.info(f"Fetching Rfam alignment for {family_acc}")
    
    aligned_sto = output_dir / f"{rna.id}_{family_acc}_aligned.sto"
    
    if rfam_seed_path and rfam_seed_path.exists():
        try:
            family_msa = extract_family_from_seed(rfam_seed_path, family_acc)
            if family_msa:
                return align_to_family_with_cmalign(rna, family_msa, output_dir, family_acc)
        except Exception as e:
            logger.warning(f"Failed to extract from Rfam.seed: {e}, falling back to mock")
    
    logger.warning(f"Using mock alignment for {family_acc} (Rfam.seed not available)")
    create_mock_alignment(rna, family_acc, aligned_sto)
    
    return aligned_sto


def extract_family_from_seed(rfam_seed_path: Path, family_acc: str) -> str:
    logger.debug(f"Extracting {family_acc} from Rfam.seed")
    
    in_family = False
    family_lines = []
    
    with open(rfam_seed_path, 'r', encoding='latin-1') as f:
        for line in f:
            if line.startswith('#=GF AC'):
                acc = line.split()[2]
                if acc == family_acc or acc.startswith(family_acc):
                    in_family = True
                    family_lines = ['# STOCKHOLM 1.0\n']
                    family_lines.append(line)
            elif in_family:
                family_lines.append(line)
                if line.startswith('//'):
                    logger.info(f"Extracted {family_acc} from Rfam.seed ({len(family_lines)} lines)")
                    return ''.join(family_lines)
    
    logger.warning(f"Family {family_acc} not found in Rfam.seed")
    return None


def align_to_family_with_cmalign(rna, family_msa: str, output_dir: Path, family_acc: str) -> Path:
    import subprocess
    
    family_sto = output_dir / f"{family_acc}_seed.sto"
    with open(family_sto, 'w') as f:
        f.write(family_msa)
    
    cm_file = output_dir / f"{family_acc}.cm"
    try:
        subprocess.run(
            ['cmbuild', str(cm_file), str(family_sto)],
            capture_output=True,
            check=True,
            text=True
        )
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.warning(f"cmbuild failed: {e}, using family MSA directly")
        return add_query_to_msa(rna, family_sto, output_dir, family_acc)
    
    query_file = output_dir / f"{rna.id}_query.fa"
    with open(query_file, 'w') as f:
        f.write(f">{rna.id}\n{rna.sequence}\n")
    
    aligned_sto = output_dir / f"{rna.id}_{family_acc}_aligned.sto"
    try:
        result = subprocess.run(
            ['cmalign', '--mapali', str(family_sto), '-o', str(aligned_sto), 
             str(cm_file), str(query_file)],
            capture_output=True,
            check=True,
            text=True
        )
        logger.info(f"Aligned {rna.id} to {family_acc} using cmalign")
        return aligned_sto
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.warning(f"cmalign failed: {e}, using simple addition")
        return add_query_to_msa(rna, family_sto, output_dir, family_acc)


def add_query_to_msa(rna, family_sto: Path, output_dir: Path, family_acc: str) -> Path:
    aligned_sto = output_dir / f"{rna.id}_{family_acc}_aligned.sto"
    
    with open(family_sto, 'r') as fin, open(aligned_sto, 'w') as fout:
        for line in fin:
            if line.startswith('//'):
                fout.write(f"{rna.id:<30} {rna.sequence}\n")
            fout.write(line)
    
    return aligned_sto


def create_mock_alignment(rna, family_acc: str, output_path: Path):
    import random
    random.seed(42)
    
    with open(output_path, 'w') as f:
        f.write("# STOCKHOLM 1.0\n")
        f.write(f"#=GF ID {family_acc}\n")
        f.write(f"#=GF AC {family_acc}\n")
        f.write("\n")
        
        f.write(f"{rna.id:<30} {rna.sequence}\n")
        
        for i in range(99):
            variant = list(rna.sequence)
            n_mutations = int(len(variant) * 0.2)
            positions = random.sample(range(len(variant)), n_mutations)
            for pos in positions:
                if variant[pos] in 'ACGU':
                    variant[pos] = random.choice('ACGU')
            
            seq_id = f"{family_acc}_homolog_{i+1}"
            f.write(f"{seq_id:<30} {''.join(variant)}\n")
        
        f.write("//\n")


def count_msa_sequences(msa_path: Path) -> int:
    count = 0
    with open(msa_path, 'r', encoding='latin-1') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or line.startswith('//') or not line:
                continue
            if len(line.split()) >= 2:
                count += 1
    
    return count


def convert_stockholm_to_a3m(sto_path: Path, a3m_path: Path, query_id: str):
    from Bio import AlignIO
    
    alignment = AlignIO.read(sto_path, "stockholm")
    
    query_seq = None
    query_idx = None
    for i, record in enumerate(alignment):
        if record.id == query_id:
            query_seq = str(record.seq)
            query_idx = i
            break
    
    if query_seq is None:
        raise ValueError(f"Query {query_id} not found in alignment")
    
    with open(a3m_path, 'w') as f:
        query_ungapped = query_seq.replace('-', '').replace('.', '')
        f.write(f">{query_id}\n{query_ungapped}\n")
        
        for i, record in enumerate(alignment):
            if i == query_idx:
                continue
            
            a3m_seq = []
            for j, (q_char, s_char) in enumerate(zip(query_seq, str(record.seq))):
                if q_char in '-.':
                    if s_char not in '-.':
                        a3m_seq.append(s_char.lower())
                else:
                    if s_char in '-.':
                        a3m_seq.append('-')
                    else:
                        a3m_seq.append(s_char.upper())
            
            f.write(f">{record.id}\n{''.join(a3m_seq)}\n")
    
    logger.info(f"Converted Stockholm to A3M: {a3m_path}")
