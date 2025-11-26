#!/usr/bin/env python3
"""
Main RNA MSA Pipeline
Orchestrates the entire workflow from sequence to 2D features
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from enum import Enum

logger = logging.getLogger(__name__)


class RNAStrategy(Enum):
    """Strategy for generating RNA MSA"""
    INFERNAL_RFAM = "infernal_rfam"      # Known RNA families
    RNACMAP_STYLE = "rnacmap_style"      # Novel RNAs
    DOMAIN_SPLIT = "domain_split"        # Long RNAs (>1000nt)
    SINGLE_SEQUENCE = "single_sequence"  # Orphan/low-homology RNAs


@dataclass
class RNASequence:
    """RNA sequence with metadata"""
    id: str
    sequence: str
    length: int
    description: Optional[str] = None
    
    @classmethod
    def from_fasta(cls, fasta_path: Path):
        """Load RNA sequences from FASTA file"""
        sequences = []
        with open(fasta_path) as f:
            current_id = None
            current_seq = []
            current_desc = None
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        seq = ''.join(current_seq)
                        sequences.append(cls(
                            id=current_id,
                            sequence=seq,
                            length=len(seq),
                            description=current_desc
                        ))
                    parts = line[1:].split(None, 1)
                    current_id = parts[0]
                    current_desc = parts[1] if len(parts) > 1 else None
                    current_seq = []
                else:
                    current_seq.append(line)
            
            if current_id:
                seq = ''.join(current_seq)
                sequences.append(cls(
                    id=current_id,
                    sequence=seq,
                    length=len(seq),
                    description=current_desc
                ))
        
        return sequences


class RNAMSAPipeline:
    """
    Main pipeline for generating RNA MSAs and extracting 2D features
    """
    
    def __init__(
        self,
        rfam_db_path: Optional[Path] = None,
        rfam_seed_path: Optional[Path] = None,
        infernal_path: str = "cmscan",
        sequence_db_path: Optional[Path] = None,
        output_dir: Path = Path("data/rna_msa_output"),
        n_cpu: int = 4
    ):
        """
        Initialize RNA MSA Pipeline
        
        Args:
            rfam_db_path: Path to Rfam covariance model database
            rfam_seed_path: Path to Rfam.seed file with family alignments
            infernal_path: Path to Infernal cmscan executable
            sequence_db_path: Path to sequence database for homology search
            output_dir: Directory for output files
            n_cpu: Number of CPUs for parallel processing
        """
        self.rfam_db_path = rfam_db_path
        self.rfam_seed_path = rfam_seed_path
        self.infernal_path = infernal_path
        self.sequence_db_path = sequence_db_path
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.n_cpu = n_cpu
        
        logger.info(f"Initialized RNA MSA Pipeline with output dir: {self.output_dir}")
    
    def determine_strategy(self, rna: RNASequence) -> RNAStrategy:
        """
        Determine the best MSA generation strategy for an RNA sequence
        
        Args:
            rna: RNA sequence object
            
        Returns:
            Appropriate RNAStrategy
        """
        # Long RNAs need domain splitting
        if rna.length > 1000:
            logger.info(f"RNA {rna.id} is long ({rna.length}nt), using DOMAIN_SPLIT strategy")
            return RNAStrategy.DOMAIN_SPLIT
        
        # Try Rfam first for shorter RNAs
        if self.rfam_db_path and self.rfam_db_path.exists():
            logger.info(f"RNA {rna.id} will try INFERNAL_RFAM strategy first")
            return RNAStrategy.INFERNAL_RFAM
        
        # Default to RNAcmap-style for novel RNAs
        logger.info(f"RNA {rna.id} using RNACMAP_STYLE strategy")
        return RNAStrategy.RNACMAP_STYLE
    
    def process_rna(
        self,
        rna: RNASequence,
        strategy: Optional[RNAStrategy] = None,
        force_strategy: bool = True
    ) -> Dict:
        """
        Process a single RNA sequence through the pipeline
        
        Args:
            rna: RNA sequence to process
            strategy: Override automatic strategy selection
            force_strategy: If True, don't fall back to other strategies
            
        Returns:
            Dictionary with MSA path, feature maps, and metadata
        """
        if strategy is None:
            strategy = self.determine_strategy(rna)
        
        logger.info(f"Processing RNA {rna.id} with strategy: {strategy.value}")
        
        result = {
            'rna_id': rna.id,
            'length': rna.length,
            'strategy': strategy.value,
            'success': False
        }
        
        try:
            if strategy == RNAStrategy.INFERNAL_RFAM:
                result.update(self._process_infernal_rfam(rna))
            elif strategy == RNAStrategy.RNACMAP_STYLE:
                result.update(self._process_rnacmap_style(rna))
            elif strategy == RNAStrategy.DOMAIN_SPLIT:
                result.update(self._process_domain_split(rna))
            elif strategy == RNAStrategy.SINGLE_SEQUENCE:
                result.update(self._process_single_sequence(rna))
            
            result['success'] = True
            
        except Exception as e:
            logger.error(f"Failed to process {rna.id} with {strategy.value}: {e}")
            result['error'] = str(e)
            
            # Fall back to single-sequence if not forcing strategy
            if not force_strategy and strategy != RNAStrategy.SINGLE_SEQUENCE:
                logger.info(f"Falling back to SINGLE_SEQUENCE for {rna.id}")
                try:
                    result.update(self._process_single_sequence(rna))
                    result['strategy'] = f"{strategy.value}_fallback_single"
                    result['success'] = True
                except Exception as e2:
                    logger.error(f"Fallback also failed for {rna.id}: {e2}")
                    result['error'] = str(e2)
        
        return result
    
    def _process_infernal_rfam(self, rna: RNASequence) -> Dict:
        """Process RNA using Infernal + Rfam"""
        from .infernal_rfam import search_rfam, fetch_rfam_alignment, count_msa_sequences
        
        # Search Rfam for matches
        rfam_hits = search_rfam(
            rna,
            self.rfam_db_path,
            self.infernal_path,
            self.output_dir,
            n_cpu=self.n_cpu
        )
        
        if not rfam_hits:
            raise ValueError(f"No Rfam hits found for {rna.id}")
        
        # Use best hit
        best_hit = rfam_hits[0]
        logger.info(f"Best Rfam hit for {rna.id}: {best_hit['family']} (E={best_hit['evalue']})")
        
        # Fetch and align to family MSA
        msa_path = fetch_rfam_alignment(
            rna,
            best_hit,
            self.output_dir,
            rfam_seed_path=self.rfam_seed_path
        )
        
        # Count actual MSA depth
        msa_depth = count_msa_sequences(msa_path)
        
        return {
            'msa_path': str(msa_path),
            'rfam_family': best_hit['family'],
            'rfam_evalue': best_hit['evalue'],
            'msa_depth': msa_depth
        }
    
    def _process_rnacmap_style(self, rna: RNASequence) -> Dict:
        """Process novel RNA using RNAcmap-style iterative approach"""
        from .rnacmap_style import build_custom_cm, iterative_msa_refinement
        
        # Build custom covariance model (skip calibration for speed)
        cm_path = build_custom_cm(
            rna,
            self.output_dir,
            n_cpu=self.n_cpu,
            calibrate=False  # Skip slow calibration for testing
        )
        
        # Iterative MSA refinement
        msa_path, msa_depth = iterative_msa_refinement(
            rna,
            cm_path,
            self.sequence_db_path,
            self.output_dir,
            n_iterations=3,
            n_cpu=self.n_cpu
        )
        
        # Fail if MSA depth is too low (< 5 sequences means no real homologs)
        if msa_depth < 5:
            raise ValueError(f"Insufficient MSA depth ({msa_depth} sequences). Need at least 5 homologs for meaningful alignment.")
        
        return {
            'msa_path': str(msa_path),
            'cm_path': str(cm_path),
            'msa_depth': msa_depth
        }
    
    def _process_domain_split(self, rna: RNASequence) -> Dict:
        """Process long RNA by splitting into domains"""
        from .domain_split import split_into_domains, process_domain
        
        # Split RNA into structural domains
        domains = split_into_domains(
            rna,
            window_size=300,
            overlap=50
        )
        
        logger.info(f"Split {rna.id} into {len(domains)} domains")
        
        # Process each domain
        domain_results = []
        for i, domain in enumerate(domains):
            logger.info(f"Processing domain {i+1}/{len(domains)} for {rna.id}")
            
            # Use single-sequence for domains (RNAcmap is too slow without sequence DB)
            try:
                # Skip RNAcmap if no sequence database
                if self.sequence_db_path and self.sequence_db_path.exists():
                    domain_result = self._process_rnacmap_style(domain)
                else:
                    domain_result = self._process_single_sequence(domain)
            except Exception as e:
                logger.warning(f"Domain {i} failed: {e}, using single-sequence")
                domain_result = self._process_single_sequence(domain)
            
            domain_result['domain_index'] = i
            domain_result['domain_start'] = domain.start
            domain_result['domain_end'] = domain.end
            domain_results.append(domain_result)
        
        return {
            'domains': domain_results,
            'n_domains': len(domains)
        }
    
    def _process_single_sequence(self, rna: RNASequence) -> Dict:
        """Process orphan RNA with single-sequence structure prediction"""
        from .single_sequence import predict_structure_spotRNA
        
        # Predict structure using SPOT-RNA or similar
        structure = predict_structure_spotRNA(
            rna,
            self.output_dir
        )
        
        return {
            'structure': structure,
            'method': 'single_sequence',
            'msa_depth': 0,
            'note': 'No MSA available - using single-sequence prediction'
        }
    
    def process_dataset(
        self,
        fasta_path: Path,
        output_features_path: Optional[Path] = None
    ) -> Dict[str, Dict]:
        """
        Process entire dataset of RNA sequences
        
        Args:
            fasta_path: Path to FASTA file with RNA sequences
            output_features_path: Optional path to save extracted features
            
        Returns:
            Dictionary mapping RNA IDs to their results
        """
        logger.info(f"Processing dataset from {fasta_path}")
        
        # Load sequences
        sequences = RNASequence.from_fasta(fasta_path)
        logger.info(f"Loaded {len(sequences)} RNA sequences")
        
        # Process each sequence
        results = {}
        for i, rna in enumerate(sequences, 1):
            logger.info(f"Processing {i}/{len(sequences)}: {rna.id}")
            result = self.process_rna(rna)
            results[rna.id] = result
        
        return results

