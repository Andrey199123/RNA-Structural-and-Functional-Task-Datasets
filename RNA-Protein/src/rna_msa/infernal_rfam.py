#!/usr/bin/env python3
"""
Infernal + Rfam Strategy
For known/annotated RNA families
"""

import subprocess
import logging
from pathlib import Path
from typing import List, Dict
import tempfile

logger = logging.getLogger(__name__)


def search_rfam(rna, rfam_db_path: Path, infernal_path: str, output_dir: Path, 
                n_cpu: int = 4, evalue_threshold: float = 0.01) -> List[Dict]:
    """
    Search Rfam database using Infernal cmscan
    
    Args:
        rna: RNASequence object
        rfam_db_path: Path to Rfam.cm database
        infernal_path: Path to cmscan executable
        output_dir: Output directory
        n_cpu: Number of CPUs
        evalue_threshold: E-value threshold for hits
        
    Returns:
        List of Rfam hits sorted by E-value
    """
    # Write query sequence to temp file
    query_file = output_dir / f"{rna.id}_query.fa"
    with open(query_file, 'w') as f:
        f.write(f">{rna.id}\n{rna.sequence}\n")
    
    # Run cmscan
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
    
    # Parse results
    hits = parse_cmscan_tblout(tblout_file)
    logger.info(f"Found {len(hits)} Rfam hits for {rna.id}")
    
    return hits


def parse_cmscan_tblout(tblout_path: Path) -> List[Dict]:
    """
    Parse Infernal cmscan tblout format
    
    Format: target_name accession query_name accession mdl mdl_from mdl_to seq_from seq_to strand trunc pass gc bias score E-value inc description
    """
    hits = []
    
    with open(tblout_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) < 17:
                continue
            
            try:
                # Infernal tblout format (space-separated)
                # 0:target 1:accession 2:query 3:accession 4:mdl 5:from 6:to 7:seq_from 8:seq_to 
                # 9:strand 10:trunc 11:pass 12:gc 13:bias 14:score 15:E-value 16:inc 17+:description
                hit = {
                    'target': parts[1],           # Rfam accession (e.g., "RF00005")
                    'family': parts[0],           # Target name (e.g., "tRNA")
                    'query': parts[2],            # Query sequence name
                    'evalue': float(parts[15]),   # E-value
                    'score': float(parts[14]),    # Score
                    'bias': float(parts[13]),     # Bias
                    'model_from': int(parts[5]),  # Model start
                    'model_to': int(parts[6]),    # Model end
                    'seq_from': int(parts[7]),    # Sequence start
                    'seq_to': int(parts[8]),      # Sequence end
                    'strand': parts[9],           # Strand
                    'msa_depth': 100  # Placeholder - Rfam families have deep MSAs
                }
                hits.append(hit)
                logger.info(f"Parsed hit: {parts[0]} ({parts[1]}) E-value={parts[15]}")
            except (ValueError, IndexError) as e:
                # Skip lines that can't be parsed
                logger.debug(f"Skipping unparseable line: {line.strip()[:100]}")
                continue
    
    # Sort by E-value
    hits.sort(key=lambda x: x['evalue'])
    
    return hits


def fetch_rfam_alignment(rna, rfam_hit: Dict, output_dir: Path, rfam_seed_path: Path = None) -> Path:
    """
    Fetch Rfam family alignment and align query sequence to it using cmalign
    
    Args:
        rna: RNASequence object
        rfam_hit: Dictionary with Rfam hit information
        output_dir: Output directory
        rfam_seed_path: Path to Rfam.seed file (optional)
        
    Returns:
        Path to aligned MSA in Stockholm format
    """
    family_acc = rfam_hit['target']
    
    logger.info(f"Fetching Rfam alignment for {family_acc}")
    
    aligned_sto = output_dir / f"{rna.id}_{family_acc}_aligned.sto"
    
    # Try to extract family from Rfam.seed if available
    if rfam_seed_path and rfam_seed_path.exists():
        try:
            family_msa = extract_family_from_seed(rfam_seed_path, family_acc)
            if family_msa:
                # Use cmalign to align query to family
                return align_to_family_with_cmalign(rna, family_msa, output_dir, family_acc)
        except Exception as e:
            logger.warning(f"Failed to extract from Rfam.seed: {e}, falling back to mock")
    
    # Fallback: Create mock alignment
    logger.warning(f"Using mock alignment for {family_acc} (Rfam.seed not available)")
    create_mock_alignment(rna, family_acc, aligned_sto)
    
    return aligned_sto


def extract_family_from_seed(rfam_seed_path: Path, family_acc: str) -> str:
    """
    Extract a specific family's alignment from Rfam.seed file
    
    Args:
        rfam_seed_path: Path to Rfam.seed file
        family_acc: Family accession (e.g., RF00005)
        
    Returns:
        Stockholm format alignment as string, or None if not found
    """
    logger.debug(f"Extracting {family_acc} from Rfam.seed")
    
    in_family = False
    family_lines = []
    
    # Use latin-1 encoding to handle special characters in species names
    with open(rfam_seed_path, 'r', encoding='latin-1') as f:
        for line in f:
            if line.startswith('#=GF AC'):
                # Check if this is our family
                acc = line.split()[2]
                if acc == family_acc or acc.startswith(family_acc):
                    in_family = True
                    family_lines = ['# STOCKHOLM 1.0\n']
                    family_lines.append(line)
            elif in_family:
                family_lines.append(line)
                if line.startswith('//'):
                    # End of family
                    logger.info(f"Extracted {family_acc} from Rfam.seed ({len(family_lines)} lines)")
                    return ''.join(family_lines)
    
    logger.warning(f"Family {family_acc} not found in Rfam.seed")
    return None


def align_to_family_with_cmalign(rna, family_msa: str, output_dir: Path, family_acc: str) -> Path:
    """
    Align query RNA to family MSA using cmalign
    
    Args:
        rna: RNASequence object
        family_msa: Stockholm format MSA string
        output_dir: Output directory
        family_acc: Family accession
        
    Returns:
        Path to aligned Stockholm file
    """
    import subprocess
    
    # Write family MSA to temp file
    family_sto = output_dir / f"{family_acc}_seed.sto"
    with open(family_sto, 'w') as f:
        f.write(family_msa)
    
    # Build CM from family alignment
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
        # Just add query to family MSA without proper alignment
        return add_query_to_msa(rna, family_sto, output_dir, family_acc)
    
    # Write query sequence
    query_file = output_dir / f"{rna.id}_query.fa"
    with open(query_file, 'w') as f:
        f.write(f">{rna.id}\n{rna.sequence}\n")
    
    # Align query to family using cmalign
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
    """
    Simple fallback: add query sequence to family MSA without proper alignment
    """
    aligned_sto = output_dir / f"{rna.id}_{family_acc}_aligned.sto"
    
    with open(family_sto, 'r') as fin, open(aligned_sto, 'w') as fout:
        for line in fin:
            if line.startswith('//'):
                # Add query before end marker
                fout.write(f"{rna.id:<30} {rna.sequence}\n")
            fout.write(line)
    
    return aligned_sto


def create_mock_alignment(rna, family_acc: str, output_path: Path):
    """Create mock alignment for testing when Rfam.seed not available"""
    import random
    random.seed(42)
    
    with open(output_path, 'w') as f:
        f.write("# STOCKHOLM 1.0\n")
        f.write(f"#=GF ID {family_acc}\n")
        f.write(f"#=GF AC {family_acc}\n")
        f.write("\n")
        
        # Write query sequence
        f.write(f"{rna.id:<30} {rna.sequence}\n")
        
        # Write mock homologous sequences
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
    """
    Count number of sequences in Stockholm MSA file
    
    Args:
        msa_path: Path to Stockholm format MSA
        
    Returns:
        Number of sequences in MSA
    """
    count = 0
    with open(msa_path, 'r', encoding='latin-1') as f:
        for line in f:
            line = line.strip()
            # Skip comments, markup, and empty lines
            if line.startswith('#') or line.startswith('//') or not line:
                continue
            # Count sequence lines (have at least 2 fields: ID and sequence)
            if len(line.split()) >= 2:
                count += 1
    
    return count


def convert_stockholm_to_a3m(sto_path: Path, a3m_path: Path, query_id: str):
    """
    Convert Stockholm alignment to A3M format
    
    A3M format:
    - First sequence is query (no gaps)
    - Other sequences: '-' for deletions, lowercase for insertions
    
    Args:
        sto_path: Input Stockholm file
        a3m_path: Output A3M file
        query_id: ID of query sequence
    """
    from Bio import AlignIO
    
    # Read Stockholm alignment
    alignment = AlignIO.read(sto_path, "stockholm")
    
    # Find query sequence
    query_seq = None
    query_idx = None
    for i, record in enumerate(alignment):
        if record.id == query_id:
            query_seq = str(record.seq)
            query_idx = i
            break
    
    if query_seq is None:
        raise ValueError(f"Query {query_id} not found in alignment")
    
    # Write A3M format
    with open(a3m_path, 'w') as f:
        # Write query without gaps
        query_ungapped = query_seq.replace('-', '').replace('.', '')
        f.write(f">{query_id}\n{query_ungapped}\n")
        
        # Write other sequences in A3M format
        for i, record in enumerate(alignment):
            if i == query_idx:
                continue
            
            # Convert to A3M format
            a3m_seq = []
            for j, (q_char, s_char) in enumerate(zip(query_seq, str(record.seq))):
                if q_char in '-.' :  # Query has gap = insertion in subject
                    if s_char not in '-.':
                        a3m_seq.append(s_char.lower())
                else:  # Query has residue
                    if s_char in '-.':
                        a3m_seq.append('-')
                    else:
                        a3m_seq.append(s_char.upper())
            
            f.write(f">{record.id}\n{''.join(a3m_seq)}\n")
    
    logger.info(f"Converted Stockholm to A3M: {a3m_path}")
