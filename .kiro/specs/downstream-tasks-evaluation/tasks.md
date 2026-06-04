# Implementation Plan: Downstream Tasks Evaluation System

## Overview

This implementation plan breaks down the development of a comprehensive evaluation system for a pretrained RNA structure prediction model across seven downstream tasks. The system follows a modular pipeline architecture with independent task pipelines that can execute in parallel. Each task downloads datasets from published papers, preprocesses them according to paper-specific methodologies, and prepares them for zero-shot evaluation or fine-tuning.

The implementation follows a bottom-up approach: first building core utilities and components, then individual task pipelines, and finally the orchestration layer that coordinates all tasks.

## Tasks

- [x] 1. Set up project structure and core utilities
  - Create directory structure for source code, data, logs, reports, and configuration
  - Implement file utilities for path handling, file I/O, and archive extraction
  - Implement sequence utilities for RNA sequence validation and manipulation
  - Set up logging configuration with file and console handlers
  - Create configuration file schemas for tasks, preprocessing, and validation
  - _Requirements: 14.5, 15.4, 15.5_

- [x] 2. Implement Dataset Downloader component
  - [x] 2.1 Create DatasetDownloader class with download and verification methods
    - Implement HTTP/FTP download with progress tracking
    - Add retry logic with exponential backoff (3 attempts)
    - Implement checksum verification for file integrity
    - Add support for resuming partial downloads
    - _Requirements: 1.1-1.7, 15.1_
  
  - [x] 2.2 Implement archive extraction functionality
    - Support .tar.gz, .zip, .tar.bz2 formats
    - Extract to specified directories with validation
    - Handle nested archives if present
    - _Requirements: 1.9_
  
  - [x] 2.3 Add download metadata tracking
    - Store source URL, download timestamp, file checksums
    - Save metadata to JSON files alongside downloaded data
    - _Requirements: 1.10, 14.2_
  
  - [x]* 2.4 Write unit tests for DatasetDownloader
    - Test download with mock URLs and simulated failures
    - Test retry logic and exponential backoff
    - Test checksum verification with valid and invalid checksums
    - Test archive extraction for all supported formats
    - _Requirements: 1.8, 15.1_

- [x] 3. Implement Preprocessor Engine component
  - [x] 3.1 Create PreprocessorEngine class with filtering methods
    - Implement sequence length filtering
    - Implement quality score filtering
    - Implement redundancy removal using CD-HIT-EST
    - Handle missing data with configurable strategies
    - _Requirements: 3.1, 3.2_
  
  - [x] 3.2 Implement normalization methods
    - Per-nucleotide normalization for RSA values
    - Distance normalization for contact maps
    - Feature scaling (z-score, min-max)
    - _Requirements: 3.3_
  
  - [x] 3.3 Implement data splitting functionality
    - Stratified splitting for classification tasks
    - Random splitting for regression tasks
    - Support for paper-specified splits (TR89/CN48, TR221/TS80)
    - _Requirements: 3.4_
  
  - [x] 3.4 Implement sequence and label encoding
    - One-hot encoding for sequences
    - Token ID encoding for model input
    - Binary, multi-class, and multi-label encoding
    - Structure encoding (dot-bracket, contact matrix)
    - _Requirements: 3.5, 4.6, 4.7_
  
  - [x]* 3.5 Write unit tests for PreprocessorEngine
    - Test each filtering operation with sample data
    - Test normalization methods with known inputs/outputs
    - Test splitting strategies for correctness and reproducibility
    - Test encoding schemes for all supported formats
    - _Requirements: 3.8, 3.9_

- [x] 4. Implement Data Validator component
  - [x] 4.1 Create DataValidator class with validation methods
    - Implement sample count validation with tolerance
    - Implement sequence length distribution checks
    - Implement class distribution validation
    - Implement data leakage detection between splits
    - Implement format and schema validation
    - _Requirements: 13.1-13.3, 13.5, 13.6_
  
  - [x] 4.2 Implement validation report generation
    - Generate detailed validation reports with pass/fail status
    - Document all checks performed and their results
    - Highlight discrepancies with paper specifications
    - _Requirements: 13.4_
  
  - [x]* 4.3 Write unit tests for DataValidator
    - Test sample count validation with various tolerances
    - Test distribution checks with matching and mismatched data
    - Test data leakage detection with overlapping and non-overlapping splits
    - Test format validation with valid and invalid schemas
    - _Requirements: 13.1-13.6_

- [x] 5. Checkpoint - Ensure core components are functional
  - Ensure all tests pass, ask the user if questions arise.

- [x] 6. Implement Task Pipeline base class
  - [x] 6.1 Create TaskPipeline base class with stage execution
    - Implement monolithic task execution with error handling
    - Persist `status.json` and support resume/force semantics (skip when prior status is SUCCESS unless forced)
    - _Requirements: 15.2, 15.3_
  
  - [x] 6.2 Implement progress tracking and logging
    - Track task start/end timestamps and status state
    - Log all operations with timestamps and context
    - Generate task-specific logs
    - _Requirements: 5.4, 14.1_
  
  - [x]* 6.3 Write unit tests for TaskPipeline base class
    - Test stage execution with mock stages
    - Test checkpoint saving and resumption
    - Test error handling and recovery
    - Test progress tracking accuracy
    - _Requirements: 15.2, 15.3_

- [x] 7. Implement Task 1: RSA Prediction Pipeline
  - [x] 7.1 Create RSATaskPipeline class extending TaskPipeline
    - **Primary benchmark (matches Azam's project notes)**: RNAsnap/RNAsol RSA benchmark (TR89 / TS44 / CN48)
    - Train: TR89 (89 protein-bound RNA chains)
    - Test: TS44 (44 protein-bound RNA chains)
    - Test: CN48 (48 protein-free RNA chains)
    - Parse TINKER ASA benchmark records and normalize ASA per nucleotide (A/G=400Å², U/C=350Å²; clamp to [0,1])
    - Output RSA data in model-compatible JSON format (normalized ASA + label mask)
    - **Fallbacks (kept for convenience only)**: RNAsnap2_data.zip (TR95/VL24/TS45/TS31) and `benchmarking/` (TS45/TS31-only)
    - **Status**: Ready for preparation
    - _Requirements: 6.1-6.4, 1.1, 3.1-3.7_
  
  - [x]* 7.2 Write integration tests for RSA pipeline
    - Test complete pipeline with sample RNAsnap2 data
    - Verify output format matches specification
    - Validate SASA values are correctly normalized
    - Check splits (TR89, TS44, CN48) are correctly extracted (val split is empty for the official benchmark)
    - _Requirements: 6.1-6.4, 13.1-13.6_

- [x] 8. Implement Task 2: Subcellular Localization Pipeline (DM3Loc/LOCAS; trunc8000)
  - [x] 8.1 Investigation + benchmark selection
    - Identify the paper-aligned benchmark dataset used by DM3Loc/LOCAS (RNALocate-derived; sequences + labels included in DM3Loc FASTA)
    - Analyze sequence length distribution (many mRNAs exceed common context limits; some extremely long)
    - Decide to proceed using a pragmatic truncation heuristic to keep the task runnable
    - _Requirements: 7.1-7.6, 1.2_

  - [x] 8.2 Data preprocessing (DM3Loc-style)
    - Use DM3Loc nonredundant FASTA with 7-bit header labels
    - Predict 6 compartments (exclude cytoplasm as an auxiliary class used in DM3Loc training)
    - Apply DM3Loc-style end-truncation to max 8000 nt (first 4000 + last 4000)
    - Drop “cytoplasm-only” records by default (expected n=17,277); allow keeping all (n=17,306) via config/flag
    - Generate deterministic 5-fold CV splits (seeded)
    - Output gzipped JSONL + folds + summary + validation report
    - _Requirements: 7.1-7.6, 3.1-3.7, 13.1-13.6_

- [x] 9. Implement Task 3: Contact Map Prediction Pipeline
  - [x] 9.1 Create ContactMapTaskPipeline class
    - **Source**: Download the public RNAcontact benchmark archives (train/test) from yanglab
    - Train set: 211 sequences; Test set: 80 sequences (as provided by the public benchmark archives)
    - Split train set into train/val deterministically (seeded) for convenience
    - Calculate minimum atom-pair distances
    - Label contacts within 8Å as positive
    - Filter for long-range contacts (|i-j| ≥ 24)
    - **Status**: Ready for preparation
    - _Requirements: 8.1-8.4, 1.3, 3.1-3.7_
  
  - [x]* 9.2 Write integration tests for contact map pipeline
    - Test complete pipeline with sample RNAcontact sequences
    - Validate contact map generation accuracy
    - Check long-range contact filtering (|i-j| ≥ 24)
    - _Requirements: 8.1-8.4, 13.1-13.6_

- [~] 10. Checkpoint - Ensure first three task pipelines are complete
  - Ensure all tests pass, ask the user if questions arise.

- [x] 11. Implement Task 4: RNA-Protein Binding Pipeline
  - [x] 11.1 Create RBPBindingTaskPipeline class
    - **Source**: Download PrismNet dataset from github.com/kuixu/PrismNet
    - Filter for 17 RBPs in HeLa cell environment
    - Extract 101 nt fixed-length sequences (1,827 to 15,002 sequences per RBP)
    - Parse HDF5 (.h5) file format
    - **CRITICAL**: Use exact same train/val/test split as ERNIE-RNA paper for fair comparison
    - Optionally extract icSHAPE features for comparison
    - Output RBP data in model-compatible format with binding site labels
    - **Status**: Ready for preparation
    - _Requirements: 9.1-9.4, 1.4, 3.1-3.7_
  
  - [x]* 11.2 Write integration tests for RBP binding pipeline
    - Test complete pipeline with sample PrismNet data
    - Verify fixed sequence length enforcement (101 nt)
    - Validate exact match with ERNIE-RNA paper train/val/test splits
    - Check icSHAPE feature extraction if enabled
    - Verify per-RBP data organization (17 RBPs)
    - _Requirements: 9.1-9.4, 13.1-13.6_

- [x] 12. Implement Task 5: RNA-RNA Interaction Pipeline
  - [x] 12.1 Create RRITaskPipeline class
    - **Source**: DeepMirTar_SdA repo snapshot (engineered feature pickles)
    - Note: this implementation prepares the available feature matrices (not the original 27,720 raw sequence pairs)
    - Output: `deepmirtar_features.npz` with stratified train/val/test indices
    - **Status**: Ready for preparation
    - _Requirements: 10.1-10.4, 1.5, 3.1-3.7_
  
  - [x]* 12.2 Write integration tests for RRI pipeline
    - Test complete pipeline with sample DeepMirTar data
    - Verify feature matrix + labels are saved and split indices are non-overlapping
    - _Requirements: 10.1-10.4, 13.1-13.6_

- [ ] 13. Task 6: Secondary Structure Prediction (ON HOLD - BLOCKED)
  - [~] 13.1 BLOCKED - Awaiting approval for SecondaryStructureTaskPipeline
    - **Status**: ON HOLD - May not be meaningful (model trained on secondary structure)
    - **DO NOT PROCEED** until explicit approval received
    - Potential datasets if approved:
      - bpRNA-1m with TR0/TV0/TS0 splits (Danaee et al., NAR 2018)
      - RNAStralign (Tan et al., NAR 2017)
      - ArchiveII (Sloma and Mathews, RNA 2016)
      - Szikszai 2022 Cross-Family Split (Szikszai et al., Bioinformatics 2022)
    - **NOTE**: Task inclusion pending decision on evaluation meaningfulness
    - _Requirements: 11.1-11.4, 1.6_
  
  - [ ]* 13.2 BLOCKED - Integration tests pending task approval
    - Cannot proceed until Task 6 is approved for implementation
    - _Requirements: 11.1-11.4, 13.1-13.6_

- [~] 14. Implement Task 7: Long-Range Base Pair Prediction Pipeline
  - [~] 14.1 Create LongRangeBasePairTaskPipeline class
    - **Dataset A (PDB Benchmark)**:
      - Download 532 non-redundant RNA sequences from PDB
      - Remove protein and DNA binding partners
      - Apply redundancy filtering using CD-HIT-EST at 80% identity
      - Annotate secondary structures using DSSR tool
      - Identify 158 sequences containing pseudoknots
      - **NOTE**: Not yet implemented paper-exact; requires the official PDB/chain manifest to reproduce 532/158 exactly.
    - **Dataset B (bpRNA-1m TR0)**:
      - Download TR0 split containing 10,814 sequences
      - Identify 999 sequences containing pseudoknots
      - Prepare exclusively for zero-shot evaluation (no fine-tuning)
    - Filter base pairs: keep only |i-j| ≥ 101 for long-range evaluation
    - Create separate evaluation set for long-range pairs
    - Output long-range base pair data in model-compatible format
    - **Status**: Ready for preparation
    - _Requirements: 12.1-12.10, 1.7, 3.1-3.7_
  
  - [x]* 14.2 Write integration tests for long-range base pair pipeline
    - Test TR0 JSONL ingestion + long-range pair export on a tiny synthetic dataset
    - Check long-range pair filtering accuracy (|i-j| ≥ min_sep)
    - _Requirements: 12.1-12.10, 13.1-13.6_

- [~] 15. Checkpoint - Ensure all seven task pipelines are complete
  - Ensure all tests pass, ask the user if questions arise.

- [x] 16. Implement Task Pipeline Orchestrator
  - [x] 16.1 Create TaskPipelineOrchestrator class
    - Implement task registry mapping task IDs to pipeline classes
    - Add execute_all_tasks method for parallel execution
    - Add execute_single_task method for individual task execution
    - Implement task status tracking with persistent storage
    - Add get_task_status and get_all_statuses methods
    - _Requirements: 5.1-5.6_
  
  - [x] 16.2 Implement consolidated reporting
    - Generate consolidated report summarizing all task statuses
    - Include completion status, elapsed time, and output locations
    - Document any failures or warnings
    - _Requirements: 5.5, 14.3_
  
  - [x]* 16.3 Write integration tests for orchestrator
    - Test parallel execution of multiple tasks
    - Test individual task execution
    - Test status tracking and persistence
    - Test consolidated report generation
    - _Requirements: 5.1-5.6_

- [x] 17. Implement configuration management
  - [x] 17.1 Create configuration file templates
    - Create task_configs.yaml with all seven task configurations
    - Create preprocessing_configs.yaml with task-specific preprocessing parameters
    - Create validation_specs.yaml with validation criteria for each task
    - _Requirements: 14.5_
  
  - [x] 17.2 Implement configuration loading and validation
    - Load YAML configuration files
    - Validate configuration schemas
    - Provide default values for optional parameters
    - Support custom user configurations
    - _Requirements: 14.5_

- [x] 18. Implement external tool integrations
  - [x] 18.1 Integrate DSSR for structure annotation
    - Check DSSR installation and version
    - Implement Python wrapper for DSSR command-line tool
    - Parse DSSR output and extract base pairs
    - **NOTE**: Required for Task 7 (Long-Range Base Pair Prediction)
    - _Requirements: 12.4_
  
  - [x] 18.2 Integrate CD-HIT-EST for redundancy removal
    - Check CD-HIT-EST installation and version
    - Implement Python wrapper for CD-HIT-EST command-line tool
    - Parse CD-HIT-EST output and extract non-redundant sequences
    - **NOTE**: Required for Task 7 (PDB Benchmark redundancy filtering at 80% identity)
    - _Requirements: 12.3_

- [x] 19. Implement reproducibility and documentation features
  - [x] 19.1 Create preprocessing log generation
    - Log every transformation with timestamps
    - Include input/output statistics for each step
    - Store logs in task-specific log files
    - _Requirements: 14.1_
  
  - [x] 19.2 Generate README files for preprocessed datasets
    - Document data format and schema
    - Describe label encoding schemes
    - Provide usage instructions and examples
    - _Requirements: 14.3_
  
  - [x] 19.3 Create reproducibility reports
    - Document all parameters and random seeds used
    - Record software versions and dependencies
    - Include hardware specifications and evaluation date
    - _Requirements: 14.4, 14.6_

- [x] 20. Implement command-line interface
  - [x] 20.1 Create CLI for orchestrator
    - Add command to execute all tasks: `python -m downstream_tasks.cli run-all`
    - Add command to execute single task: `python -m downstream_tasks.cli run-task <task_id>`
    - Add command to check task status: `python -m downstream_tasks.cli status [task_id]`
    - Add command to generate reports: `python -m downstream_tasks.cli report`
    - _Requirements: 5.1, 5.2, 5.3, 5.5_
  
  - [x] 20.2 Add CLI options for configuration
    - Support custom configuration file paths
    - Support output directory specification
    - Support parallel vs sequential execution mode
    - Support resume from checkpoint flag
    - _Requirements: 14.5, 15.3_

- [x] 21. Final integration and end-to-end testing
  - [x] 21.1 Create end-to-end test suite
    - Add smoke tests for key pipelines (Task3/4/5/7) + resume/consolidation paths
    - _Requirements: 5.1-5.6, 15.1-15.3_
  
  - [x] 21.2 Perform validation against paper specifications
    - Generate per-task validation reports (counts, leakage checks, basic schema checks)
    - _Requirements: 13.1-13.4_
  
  - [x] 21.3 Generate final documentation
    - Create comprehensive README with installation instructions
    - Document all configuration options
    - Provide usage examples for common scenarios
    - Include troubleshooting guide
    - _Requirements: 14.3, 14.5_

- [x] 22. Final checkpoint - System ready for deployment
  - Ensure all tests pass, ask the user if questions arise.

## Notes

- Tasks marked with `*` are optional test tasks and can be skipped for faster MVP delivery
- Each task references specific requirements for traceability
- Checkpoints ensure incremental validation at key milestones
- The implementation follows a bottom-up approach: utilities → components → task pipelines → orchestration

**Key Updates Based on Azam's Notes**:
- **Task 1 (RSA)**: Uses the RNAsnap/RNAsol RSA benchmark (TR89 / TS44 / CN48); labels provided; no FreeSASA
- **Task 2 (Subcellular Localization)**: Implemented using DM3Loc/LOCAS-style mRNA benchmark with DM3Loc end-truncation to 8000 nt (first 4000 + last 4000); requires Azam sign-off as an acceptable comparison heuristic
- **Task 3 (Contact Map)**: Uses the public RNAcontact benchmark archives (train=211, test=80) and a deterministic split for train/val
- **Task 4 (RBP Binding)**: CRITICAL - Must use exact same train/val/test split as ERNIE-RNA paper
- **Task 5 (RRI)**: Current implementation prepares available engineered feature pickles (not raw sequence pairs)
- **Task 6 (Secondary Structure)**: ON HOLD - DO NOT PROCEED until explicit approval (model trained on secondary structure)
- **Task 7 (Long-Range Base Pairs)**: TR0 (10,814 seqs) implemented; PDB Benchmark remains pending the official manifest to reproduce 532/158 exactly

**Task Priorities**:
- Ready for preparation: Tasks 1, 2, 3, 4, 5, 7
- Blocked pending approval: Task 6

**Critical Requirements**:
- Tasks 3, 4, 5 must match exact paper splits for fair comparison
- Task 6 requires explicit approval before implementation
- Send completed tasks to Azam incrementally (not all at once)

- All task pipelines (except Task 6) are independent and can be developed in parallel after core components are complete
- External tool integrations (DSSR, CD-HIT-EST) should be tested with mock data before full pipeline integration
- Configuration files provide flexibility while maintaining paper methodology defaults
- Comprehensive logging and documentation enable reproducibility and debugging

## Task Dependency Graph

```json
{
  "waves": [
    { "id": 0, "tasks": ["1"] },
    { "id": 1, "tasks": ["2.1", "2.2", "2.3", "3.1", "3.2", "3.3", "3.4", "4.1", "4.2"] },
    { "id": 2, "tasks": ["2.4", "3.5", "4.3", "6.1", "6.2"] },
    { "id": 3, "tasks": ["6.3", "7.1", "8.1", "9.1"] },
    { "id": 4, "tasks": ["7.2", "8.2", "9.2", "11.1", "12.1", "14.1"] },
    { "id": 5, "tasks": ["11.2", "12.2", "14.2", "16.1", "16.2", "17.1", "17.2", "18.1", "18.2"] },
    { "id": 6, "tasks": ["16.3", "19.1", "19.2", "19.3", "20.1", "20.2"] },
    { "id": 7, "tasks": ["21.1", "21.2", "21.3"] }
  ]
}
```

**Note**: Task 13 (Secondary Structure - Task 6) is excluded from the dependency graph as it is ON HOLD pending approval. If approved, it would be added to wave 4.
