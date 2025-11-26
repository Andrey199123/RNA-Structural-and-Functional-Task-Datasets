#!/usr/bin/env python3
"""
Optimized homology search using BLAST pre-filtering
Two-stage approach: BLAST (fast) → cmsearch (accurate)
"""

import subprocess
from pathlib import Path
from typing import List, Optional, Tuple
import tempfile
import logging

logger = logging.getLogger(__name__)


class FastHomologySearch:
    """
    Fast homology search using BLAST pre-filtering
    """
    
    def __init__(self, sequence_db: Path, blast_db: Optional[Path] = None):
        """
        Initialize fast homology search
        
        Args:
            sequence_db: Path to sequence database (FASTA)
            blast_db: Path to BLAST database (if None, will create)
        """
        self.sequence_db = Path(sequence_db)
        self.blast_db = blast_db or self.sequence_db.with_suffix('')
        
        # Check if BLAST database exists
        if not self._blast_db_exists():
            logger.info(f"BLAST database not found, will create on first search")
    
    def _blast_db_exists(self) -> bool:
        """Check if BLAST database files exist"""
        extensions = ['.nhr', '.nin', '.nsq']
        return all((Path(str(self.blast_db) + ext).exists() for ext in extensions))
    
    def create_blast_db(self, force: bool = False):
        """
        Create BLAST database from sequence database
        
        Args:
            force: Force recreation even if exists
        """
        if self._blast_db_exists() and not force:
            logger.info(f"BLAST database already exists: {self.blast_db}")
            return
        
        logger.info(f"Creating BLAST database from {self.sequence_db}...")
        logger.info("This may take 10-20 minutes for large databases...")
        
        cmd = [
            'makeblastdb',
            '-in', str(self.sequence_db),
            '-dbtype', 'nucl',
            '-out', str(self.blast_db),
            '-parse_seqids'
        ]
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            logger.info("✓ BLAST database created successfully")
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to create BLAST database: {e.stderr}")
            raise
        except FileNotFoundError:
            logger.error("makeblastdb not found. Install BLAST+: conda install -c bioconda blast")
            raise
    
    def blast_search(
        self,
        query_seq: str,
        query_id: str,
        max_targets: int = 1000,
        evalue: float = 0.001,
        min_identity: float = 0.5
    ) -> List[str]:
        """
        Fast BLAST search to find candidate homologs
        
        Args:
            query_seq: Query RNA sequence
            query_id: Query identifier
            max_targets: Maximum number of hits to return
            evalue: E-value threshold
            min_identity: Minimum sequence identity (0-1)
            
        Returns:
            List of sequence IDs of candidate homologs
        """
        # Ensure BLAST database exists
        if not self._blast_db_exists():
            self.create_blast_db()
        
        # Create temporary query file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write(f">{query_id}\n{query_seq}\n")
            query_file = f.name
        
        try:
            # Run BLAST
            cmd = [
                'blastn',
                '-query', query_file,
                '-db', str(self.blast_db),
                '-evalue', str(evalue),
                '-outfmt', '6 sseqid pident',  # Tab format: subject_id, percent_identity
                '-max_target_seqs', str(max_targets),
                '-num_threads', '4'
            ]
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            # Parse results
            candidates = []
            for line in result.stdout.strip().split('\n'):
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= 2:
                    seq_id = parts[0]
                    identity = float(parts[1]) / 100.0  # Convert percentage to fraction
                    
                    if identity >= min_identity:
                        candidates.append(seq_id)
            
            logger.info(f"BLAST found {len(candidates)} candidates (identity >= {min_identity*100}%)")
            return candidates
            
        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST search failed: {e.stderr}")
            return []
        except FileNotFoundError:
            logger.error("blastn not found. Install BLAST+: conda install -c bioconda blast")
            return []
        finally:
            # Clean up temp file
            Path(query_file).unlink(missing_ok=True)
    
    def extract_sequences(self, seq_ids: List[str], output_file: Path) -> int:
        """
        Extract sequences from database by ID
        
        Args:
            seq_ids: List of sequence IDs to extract
            output_file: Output FASTA file
            
        Returns:
            Number of sequences extracted
        """
        if not seq_ids:
            return 0
        
        # Create ID list file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            for seq_id in seq_ids:
                f.write(f"{seq_id}\n")
            id_file = f.name
        
        try:
            # Use blastdbcmd to extract sequences
            cmd = [
                'blastdbcmd',
                '-db', str(self.blast_db),
                '-entry_batch', id_file,
                '-out', str(output_file)
            ]
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            # Count extracted sequences
            n_seqs = 0
            with open(output_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        n_seqs += 1
            
            logger.info(f"Extracted {n_seqs} sequences to {output_file}")
            return n_seqs
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Sequence extraction failed: {e.stderr}")
            return 0
        except FileNotFoundError:
            logger.error("blastdbcmd not found. Install BLAST+: conda install -c bioconda blast")
            return 0
        finally:
            Path(id_file).unlink(missing_ok=True)
    
    def fast_homology_search(
        self,
        query_seq: str,
        query_id: str,
        output_fasta: Path,
        max_homologs: int = 1000,
        min_identity: float = 0.5
    ) -> Tuple[int, Path]:
        """
        Complete fast homology search pipeline
        
        Args:
            query_seq: Query RNA sequence
            query_id: Query identifier
            output_fasta: Output file for homolog sequences
            max_homologs: Maximum number of homologs
            min_identity: Minimum sequence identity
            
        Returns:
            Tuple of (number of homologs found, output file path)
        """
        logger.info(f"Fast homology search for {query_id}...")
        
        # Stage 1: BLAST pre-filter
        candidates = self.blast_search(
            query_seq,
            query_id,
            max_targets=max_homologs * 2,  # Get extra candidates
            min_identity=min_identity
        )
        
        if not candidates:
            logger.warning(f"No homologs found for {query_id}")
            return 0, output_fasta
        
        # Stage 2: Extract sequences
        n_extracted = self.extract_sequences(candidates[:max_homologs], output_fasta)
        
        return n_extracted, output_fasta


def check_blast_available() -> bool:
    """Check if BLAST tools are available"""
    tools = ['makeblastdb', 'blastn', 'blastdbcmd']
    
    for tool in tools:
        try:
            subprocess.run(
                [tool, '-version'],
                capture_output=True,
                check=True
            )
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False
    
    return True


if __name__ == "__main__":
    # Test BLAST availability
    import sys
    
    if check_blast_available():
        print("✓ BLAST tools are available")
    else:
        print("✗ BLAST tools not found")
        print("Install with: conda install -c bioconda blast")
        sys.exit(1)
    
    # Test database creation
    if len(sys.argv) > 1:
        db_path = Path(sys.argv[1])
        searcher = FastHomologySearch(db_path)
        searcher.create_blast_db()
