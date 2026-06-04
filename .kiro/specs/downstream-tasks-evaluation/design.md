# Design Document

## Overview

This design document describes the architecture and implementation approach for a comprehensive evaluation system that prepares seven downstream task datasets for evaluating a pretrained RNA structure prediction model. The system follows a modular pipeline architecture where each task operates independently, allowing parallel execution and incremental completion.

## Architecture

The system follows a modular pipeline architecture with independent task pipelines coordinated by a central orchestrator. Each task pipeline executes a standardized workflow (Download → Extract → Filter → Normalize → Split → Encode → Validate) while applying task-specific preprocessing methodologies from published papers.

**Key Architectural Principles**:
- **Modularity**: Each task pipeline is independent and self-contained
- **Parallelism**: Tasks can execute concurrently without dependencies
- **Checkpointing**: Each pipeline stage saves state for recovery
- **Paper Fidelity**: Preprocessing exactly matches published methodologies
- **Incremental Delivery**: Tasks complete and deliver independently

## Components and Interfaces

The system consists of five core components:

1. **Task Pipeline Orchestrator**: Coordinates execution of all task pipelines
   - Interface: `execute_all_tasks()`, `execute_single_task(task_id)`, `get_task_status(task_id)`
   
2. **Task Pipeline**: Executes complete workflow for a single task
   - Interface: `execute()`, `get_progress()`, `resume_from_checkpoint(checkpoint)`
   
3. **Dataset Downloader**: Retrieves datasets with retry logic
   - Interface: `download(url, output_path)`, `verify_integrity(file_path, checksum)`, `extract_archive(archive_path)`
   
4. **Preprocessor Engine**: Applies paper-specific transformations
   - Interface: `preprocess(input_data)`, `apply_filters(data)`, `normalize(data)`, `create_splits(data)`, `encode(data)`
   
5. **Data Validator**: Ensures quality and paper compliance
   - Interface: `validate(data, expected)`, `check_sample_counts(data, expected)`, `check_data_leakage(train, test)`

## Data Models

**TaskConfig**:
```python
{
  "task_id": str,           # e.g., "task1_rsa"
  "task_name": str,         # e.g., "RNA Solvent Accessibility"
  "paper_source": str,      # e.g., "RNAsnap2 (Sun et al., 2019)"
  "dataset_url": str,
  "preprocessing_config": PreprocessingConfig,
  "output_dir": Path,
  "requires_training_set": bool
}
```

**PreprocessingConfig**:
```python
{
  "sequence_length_min": Optional[int],
  "sequence_length_max": Optional[int],
  "similarity_threshold": Optional[float],
  "normalization_method": Optional[str],
  "split_ratios": Optional[Tuple[float, float, float]],
  "encoding_scheme": str,
  "distance_threshold": Optional[float]
}
```

**ValidationResult**:
```python
{
  "passed": bool,
  "checks_performed": List[str],
  "passed_checks": List[str],
  "failed_checks": List[str],
  "warnings": List[str],
  "statistics": Dict[str, Any],
  "discrepancies": List[str]
}
```

## Correctness Properties

The system must satisfy the following correctness properties to ensure reliable and reproducible dataset preparation:

### Property 1: Paper Fidelity
- For all tasks, preprocessed datasets SHALL match paper-reported statistics within 5% tolerance
- Train/test splits SHALL exactly match paper specifications when explicitly defined (Tasks 3, 4, 5)
- Preprocessing steps SHALL be applied in the same order as described in source papers

**Validates: Requirements 2.1-2.7, 3.1-3.4, 13.1-13.3**

### Property 2: Data Integrity
- No data leakage SHALL exist between train and test sets (verified by sequence overlap detection)
- All sequences SHALL contain only valid RNA nucleotides (A, U, G, C, and ambiguity codes)
- Label counts SHALL match sequence counts for all samples

**Validates: Requirements 13.5, 13.6**

### Property 3: Reproducibility
- Given the same input data and configuration, the system SHALL produce identical output
- All random operations SHALL use configurable seeds
- All preprocessing steps SHALL be logged with sufficient detail for reproduction

**Validates: Requirements 14.1-14.6**

### Property 4: Completeness
- All required fields SHALL be present in output data
- No samples SHALL be silently dropped without logging
- All paper-specified evaluation metrics SHALL be computable from output data

**Validates: Requirements 4.3-4.5, 13.6**

### Property 5: Checkpoint Consistency
- Pipeline state SHALL be recoverable from any checkpoint
- Resumed pipelines SHALL produce identical results to uninterrupted execution
- Checkpoint files SHALL include all information needed for resumption

**Validates: Requirements 15.2, 15.3**

## System Architecture

### High-Level Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                   Evaluation System                          │
│                                                              │
│  ┌────────────────────────────────────────────────────┐    │
│  │         Task Pipeline Orchestrator                  │    │
│  │  - Manages 7 independent task pipelines            │    │
│  │  - Tracks completion status                        │    │
│  │  - Enables parallel execution                      │    │
│  └────────────────────────────────────────────────────┘    │
│                          │                                   │
│         ┌────────────────┼────────────────┐                │
│         ▼                ▼                ▼                 │
│  ┌──────────┐    ┌──────────┐    ┌──────────┐            │
│  │  Task 1  │    │  Task 2  │ ...│  Task 7  │            │
│  │   RSA    │    │ SubCell  │    │Long-Range│            │
│  └──────────┘    └──────────┘    └──────────┘            │
│       │               │                 │                   │
│       └───────────────┴─────────────────┘                  │
│                       │                                     │
│         ┌─────────────┴─────────────┐                     │
│         ▼                           ▼                      │
│  ┌─────────────┐            ┌─────────────┐              │
│  │  Dataset    │            │ Preprocessor │              │
│  │ Downloader  │            │   Engine     │              │
│  └─────────────┘            └─────────────┘              │
│         │                           │                      │
│         └───────────┬───────────────┘                     │
│                     ▼                                      │
│            ┌─────────────────┐                            │
│            │  Data Validator │                            │
│            └─────────────────┘                            │
└─────────────────────────────────────────────────────────────┘
```

### Component Design

#### 1. Task Pipeline Orchestrator

**Purpose**: Coordinates execution of all seven task pipelines and manages their lifecycle.

**Responsibilities**:
- Initialize and launch task pipelines
- Track completion status of each task
- Enable parallel execution without blocking
- Generate consolidated progress reports
- Handle task-level failures gracefully

**Interface**:
```python
class TaskPipelineOrchestrator:
    def execute_all_tasks(self) -> Dict[str, TaskStatus]
    def execute_single_task(self, task_id: str) -> TaskStatus
    def get_task_status(self, task_id: str) -> TaskStatus
    def get_all_statuses(self) -> Dict[str, TaskStatus]
```

**Implementation Notes**:
- Uses a task registry to map task IDs to pipeline configurations
- Supports both sequential and parallel execution modes
- Maintains a status file for each task to enable recovery
- Logs all operations to a central orchestrator log

#### 2. Task Pipeline

**Purpose**: Executes the complete workflow for a single downstream task from download through evaluation preparation.

**Responsibilities**:
- Download dataset using Dataset Downloader
- Extract paper methodology configuration
- Preprocess data using Preprocessor Engine
- Validate output using Data Validator
- Generate task completion report

**Interface**:
```python
class TaskPipeline:
    def __init__(self, task_config: TaskConfig)
    def execute(self) -> TaskResult
    def get_progress(self) -> ProgressInfo
    def resume_from_checkpoint(self, checkpoint: str) -> TaskResult
```

**Task Configuration Schema**:
```python
@dataclass
class TaskConfig:
    task_id: str  # e.g., "task1_rsa"
    task_name: str  # e.g., "RNA Solvent Accessibility"
    paper_source: str  # e.g., "ERNIE-RNA"
    dataset_url: str
    preprocessing_config: PreprocessingConfig
    output_dir: Path
    requires_training_set: bool
```

**Pipeline Stages**:
1. **Download Stage**: Retrieve dataset files
2. **Extraction Stage**: Parse and extract relevant data
3. **Preprocessing Stage**: Apply paper-specific transformations
4. **Validation Stage**: Verify data quality and format
5. **Finalization Stage**: Generate metadata and reports

**Checkpointing**:
- Each stage saves a checkpoint file upon completion
- Failed pipelines can resume from the last successful checkpoint
- Checkpoint files include stage name, timestamp, and output paths

#### 3. Dataset Downloader

**Purpose**: Retrieves datasets from paper repositories, supplementary materials, or public databases.

**Responsibilities**:
- Download files from URLs with retry logic
- Verify file integrity using checksums
- Handle different data sources (HTTP, FTP, APIs)
- Cache downloaded files to avoid redundant downloads
- Extract compressed archives

**Interface**:
```python
class DatasetDownloader:
    def download(self, url: str, output_path: Path) -> DownloadResult
    def verify_integrity(self, file_path: Path, expected_checksum: str) -> bool
    def extract_archive(self, archive_path: Path, extract_dir: Path) -> List[Path]
```

**Download Strategy**:
- Implements exponential backoff retry (3 attempts)
- Supports resume for partial downloads
- Validates file size and format before marking complete
- Stores download metadata (URL, timestamp, checksum)

**Supported Formats**:
- Compressed archives: `.tar.gz`, `.zip`, `.tar.bz2`
- Data formats: `.fasta`, `.pdb`, `.csv`, `.json`, `.txt`
- Structured formats: `.h5`, `.pkl`, `.npz`

#### 4. Preprocessor Engine

**Purpose**: Transforms raw datasets into model-compatible format following paper-specific methodologies.

**Responsibilities**:
- Apply filtering criteria (sequence length, quality)
- Normalize data (SASA values, contact distances)
- Create train/validation/test splits
- Encode sequences and labels
- Apply data augmentation if specified
- Generate preprocessed output files

**Interface**:
```python
class PreprocessorEngine:
    def __init__(self, config: PreprocessingConfig)
    def preprocess(self, input_data: Path) -> PreprocessedData
    def apply_filters(self, data: RawData) -> FilteredData
    def normalize(self, data: FilteredData) -> NormalizedData
    def create_splits(self, data: NormalizedData) -> DataSplits
    def encode(self, data: DataSplits) -> EncodedData
```

**Preprocessing Configuration Schema**:
```python
@dataclass
class PreprocessingConfig:
    sequence_length_min: Optional[int]
    sequence_length_max: Optional[int]
    similarity_threshold: Optional[float]  # For redundancy removal
    normalization_method: Optional[str]
    split_ratios: Optional[Tuple[float, float, float]]  # train, val, test
    encoding_scheme: str  # e.g., "one-hot", "token-ids"
    augmentation_methods: List[str]
    distance_threshold: Optional[float]  # For contact maps
```

**Preprocessing Operations**:

1. **Filtering**:
   - Length-based filtering
   - Quality score filtering
   - Redundancy removal (CD-HIT, BLAST)
   - Missing data handling

2. **Normalization**:
   - Per-nucleotide normalization (RSA)
   - Distance normalization (contact maps)
   - Feature scaling (z-score, min-max)

3. **Splitting**:
   - Stratified splitting for classification tasks
   - Random splitting for regression tasks
   - Preserve paper-specified splits when available

4. **Encoding**:
   - Sequence encoding (one-hot, k-mer, token IDs)
   - Label encoding (binary, multi-class, multi-label)
   - Structure encoding (dot-bracket, contact matrix)

#### 5. Data Validator

**Purpose**: Ensures preprocessed datasets meet quality standards and match paper specifications.

**Responsibilities**:
- Validate sample counts against paper reports
- Check sequence length distributions
- Verify class balance and distributions
- Detect data leakage between splits
- Validate data format and schema
- Generate validation reports

**Interface**:
```python
class DataValidator:
    def validate(self, data: PreprocessedData, expected: ValidationSpec) -> ValidationResult
    def check_sample_counts(self, data: PreprocessedData, expected: int, tolerance: float) -> bool
    def check_length_distribution(self, data: PreprocessedData, expected: LengthDist) -> bool
    def check_class_distribution(self, data: PreprocessedData, expected: ClassDist) -> bool
    def check_data_leakage(self, train: Dataset, test: Dataset) -> bool
    def check_format(self, data: PreprocessedData, schema: DataSchema) -> bool
```

**Validation Checks**:

1. **Quantitative Checks**:
   - Sample count validation (±5% tolerance)
   - Sequence length statistics (min, max, mean, std)
   - Class distribution comparison
   - Missing value detection

2. **Qualitative Checks**:
   - Data format compliance
   - Required field presence
   - Data type correctness
   - Value range validation

3. **Integrity Checks**:
   - No duplicate sequences in test set
   - No overlap between train and test
   - Consistent label encoding
   - Valid sequence characters

**Validation Report**:
```python
@dataclass
class ValidationResult:
    passed: bool
    checks_performed: List[str]
    passed_checks: List[str]
    failed_checks: List[str]
    warnings: List[str]
    statistics: Dict[str, Any]
    discrepancies: List[str]
```

## Task-Specific Designs

### Task 1: RNA Solvent Accessibility (RSA) Prediction

**Dataset**: TR89 (training: 89 RNA chains), CN48 (testing: 48 protein-free RNAs), TS44 (testing: 44 protein-bound RNAs)

**Data Source**: RNAsnap2 benchmark (Yang et al., RNA 2017; Sun et al., Bioinformatics 2019)
**Repository**: github.com/jaswindersingh2/RNAsnap2

**Preprocessing Pipeline**:
1. Download RNAsnap2 benchmark dataset from GitHub repository
2. Extract TR89 training split (89 RNA chains)
3. Extract CN48 test split (48 protein-free RNA chains)
4. Extract TS44 test split (44 protein-bound RNA chains)
5. Parse FASTA sequences paired with position-wise float labels
6. Normalize SASA values per nucleotide type (A, U, G, C)

**Output Format**:
```python
{
    "sequence": "AUGC...",
    "sasa_values": [0.45, 0.23, ...],  # Normalized per position
    "split": "TR89" | "CN48" | "TS44",
    "rna_type": "protein-free" | "protein-bound",
    "source_id": "..."
}
```

**Special Considerations**:
- Three test sets: CN48 (protein-free) and TS44 (protein-bound) for comprehensive evaluation
- Input format: FASTA sequences with position-wise float labels
- Evaluation metric: Pearson correlation
- No FreeSASA computation needed (labels already provided)

### Task 2: RNA Subcellular Localization (DM3Loc/LOCAS; trunc8000)

**Dataset**: DM3Loc/LOCAS-style mRNA multi-label localization benchmark derived from RNALocate

**Data Source**: DM3Loc benchmark FASTA (labels embedded in header; nonredundant file has 17,306 records)

**Status**: Implemented with DM3Loc-style end-truncation (max 8000 nt)

**Preprocessing Pipeline**:
1. Load DM3Loc nonredundant FASTA where each header begins with a 7-bit label string
2. Map 7-bit label vector to 6 prediction targets (exclude cytoplasm as an auxiliary class)
3. Drop records with zero positives among the 6 targets by default (removes “cytoplasm-only” entries; expected n=17,277)
4. Apply DM3Loc-style truncation:
   - If length > 8000: keep first 4000 nt + last 4000 nt (concatenate)
   - Else: keep full sequence (store unpadded; pad in dataloader)
5. Write gzipped JSONL records + deterministic 5-fold CV split manifest (seeded)
6. Emit summary and validation report (counts + max length after truncation)

**Output Format**:
```python
{
    "id": "...",
    "sequence": "AUGC...",          # truncated to <= 8000 nt
    "labels_6": [0, 1, 0, 0, 0, 1], # multi-label targets
    "label_names_6": ["Nucleus", "Exosome", "Cytosol", "Ribosome", "Membrane", "ER"],
    "length_original": 12345,
    "length_used": 8000,
    "truncation": "end_concat_4000_4000" | "none"
}
```

**Special Considerations**:
- Many mRNAs are far longer than typical model context limits; truncation is a pragmatic heuristic to keep the task runnable.
- This may not exactly match every paper’s evaluation protocol; Azam should explicitly approve trunc8000 as an acceptable comparison heuristic.
- Multi-label classification task (functional; not geometric).

### Task 3: RNA Contact Map Prediction

**Dataset**: RNAcontact (TR221 training, TR168 and VL43 validation, TS80 testing)

**Data Source**: RNAcontact dataset (Sun et al., Bioinformatics 2021)
**Initial Dataset**: 1,786 RNA sequences from PDB
**Final Dataset**: 301 filtered sequences

**Preprocessing Pipeline**:
1. Download RNAcontact dataset
2. Extract 301 filtered sequences from initial 1,786 PDB sequences
3. Create TR221 training split
4. Create TR168 and VL43 validation splits
5. Create TS80 test split
6. Calculate minimum atom-pair distances
7. Label contacts within 8Å as positive
8. Filter for long-range contacts (|i-j| ≥ 24)

**Output Format**:
```python
{
    "sequence": "AUGC...",
    "contact_map": [[0, 0, 1, ...], ...],  # Binary matrix
    "distance_map": [[0.0, 5.2, 7.8, ...], ...],  # Distance matrix
    "pdb_id": "2cjr",
    "length": 156,
    "split": "TR221" | "TR168" | "VL43" | "TS80"
}
```

**Special Considerations**:
- Contact definition: Minimum atom-pair distance below 8Å
- Only evaluate long-range contacts (|i-j| ≥ 24)
- Evaluation: Long-Range Top-L precision (L/10, L/5, L/2, L/1)
- Requires ResNet decoder architecture
- **Status**: Ready for preparation

### Task 4: RNA-Protein Binding Prediction

**Dataset**: PrismNet dataset (17 RBPs in HeLa cells)

**Data Source**: PrismNet benchmark (Xu et al., Nucleic Acids Research 2023)
**Repository**: github.com/kuixu/PrismNet

**Preprocessing Pipeline**:
1. Download PrismNet benchmark dataset from GitHub
2. Filter for 17 RBPs in HeLa environment
3. Extract 101 nt fixed-length sequences
4. Parse HDF5 (.h5) file format
5. Split each RBP: Use exact same train/val/test split as ERNIE-RNA paper
6. Optionally extract icSHAPE features for comparison

**Output Format**:
```python
{
    "sequence": "AUGC...",  # 101 nt fixed length
    "binding_labels": [0, 1, 1, 0, ...],  # Per-nucleotide binding
    "rbp_name": "ELAVL1",
    "cell_type": "HeLa",
    "icshape_features": [0.23, 0.45, ...],  # Optional
    "split": "train" | "val" | "test"
}
```

**Special Considerations**:
- Fixed sequence length (101 nt)
- Per-RBP evaluation (17 separate models)
- Sequences per RBP: 1,827 to 15,002
- File format: HDF5 (.h5)
- **CRITICAL**: Must use exact same train/val/test split as ERNIE-RNA paper for fair comparison
- Evaluation: Mean AUPRC across all 17 proteins
- Two model variants: with/without icSHAPE
- **Status**: Ready for preparation

### Task 5: RNA-RNA Interaction

**Dataset**: DeepMirTar (27,720 total pairs: 13,860 positive + 13,860 negative)

**Data Source**: DeepMirTar dataset (Wen et al., Bioinformatics 2018)

**Preprocessing Pipeline**:
1. Download DeepMirTar dataset
2. Extract 27,720 miRNA-mRNA pairs (13,860 positive, 13,860 negative)
3. Filter: miRNA <26 nt, mRNA <53 nt (3' UTR)
4. Parse tab-separated interaction records with labels
5. Extract interaction labels (positive/negative)
6. Use exact same train/val/test split as RNAErnie paper
7. Identify canonical vs non-canonical seed types

**Output Format**:
```python
{
    "mirna_sequence": "AUGC...",  # <26 nt
    "mrna_sequence": "AUGC...",  # <53 nt (3' UTR)
    "interaction": 1,  # 1 = interact, 0 = no interaction
    "seed_type": "canonical" | "non-canonical",
    "mirna_id": "hsa-miR-21",
    "mrna_id": "ENST00000123456",
    "split": "train" | "val" | "test"
}
```

**Special Considerations**:
- Total interaction pairs: 27,720 (balanced: 13,860 positive, 13,860 negative)
- Paired sequence input
- Input format: Tab-separated interaction records
- **CRITICAL**: Must use exact same train/val/test split as RNAErnie paper for fair comparison
- Requires hybrid CNN + Bi-LSTM architecture
- Evaluation: Accuracy, Recall, Precision, F1, AUC
- **Status**: Ready for preparation

### Task 6: RNA Secondary Structure Prediction (ON HOLD)

**Status**: ON HOLD - Task may not be meaningful since pretrained model was already trained on secondary structure data

**Dataset**: bpRNA-1m (TR0/TV0/TS0), RNAStralign, ArchiveII, Szikszai 2022 Cross-Family Split

**Data Sources**:
- Danaee et al., Nucleic Acids Research 2018 (bpRNA-1m with TR0, TV0, TS0 splits)
- Tan et al., Nucleic Acids Research 2017 (RNAStralign)
- Sloma and Mathews, RNA 2016 (ArchiveII)
- Szikszai et al., Bioinformatics 2022 (Cross-Family Split)

**Preprocessing Pipeline** (if task proceeds):
1. Download bpRNA-1m with existing TR0/TV0/TS0 splits
2. Download RNAStralign for additional training
3. Download ArchiveII and create ArchiveII600 subset (<600 bp)
4. Download Szikszai 2022 Cross-Family Split
5. Extract secondary structure using DSSR tool
6. Encode structures as contact matrices

**Output Format**:
```python
{
    "sequence": "AUGC...",
    "structure_dotbracket": "(((...)))...",
    "base_pairs": [(0, 8), (1, 7), (2, 6)],  # 0-indexed pairs
    "contact_matrix": [[0, 0, 1, ...], ...],
    "family": "tRNA",
    "length": 76
}
```

**Special Considerations**:
- **CRITICAL**: Task on hold - may not be meaningful for evaluation
- Model was already trained on secondary structure data
- Awaiting decision on task inclusion
- Multiple datasets for training/testing
- Requires DSSR tool for structure extraction
- Post-processing: E2Efold method
- Evaluation: Precision, Recall, F1 (binary and macro-avg)
- **Status**: DO NOT PROCEED until explicit approval

### Task 7: Long-Range Base Pair Prediction

**Dataset A**: PDB Benchmark (532 non-redundant RNA sequences)
**Dataset B**: bpRNA-1m TR0 (10,814 sequences)

**Data Sources**: RNA-km paper datasets

**Preprocessing Pipeline**:
1. **Dataset A (PDB Benchmark)**:
   - Download 532 non-redundant RNA sequences from PDB
   - Remove protein and DNA binding partners
   - Apply redundancy filtering using CD-HIT-EST at 80% identity
   - Annotate secondary structures using DSSR tool
   - Identify 158 sequences containing pseudoknots

2. **Dataset B (bpRNA-1m TR0)**:
   - Download TR0 split containing 10,814 sequences
   - Identify 999 sequences containing pseudoknots
   - Prepare exclusively for zero-shot evaluation (no fine-tuning)

3. Filter base pairs: keep only |i-j| ≥ 101 for long-range evaluation
4. Create separate evaluation set for long-range pairs

**Output Format**:
```python
{
    "sequence": "AUGC...",
    "long_range_pairs": [(10, 150), (25, 200)],  # |i-j| ≥ 101
    "all_pairs": [(0, 8), (10, 150), ...],
    "contains_pseudoknots": True,
    "dataset": "PDB" | "TR0",
    "length": 250
}
```

**Special Considerations**:
- Two separate datasets: PDB Benchmark (532 seqs) and TR0 (10,814 seqs)
- PDB: 158 sequences with pseudoknots
- TR0: 999 sequences with pseudoknots
- Redundancy filtering: CD-HIT-EST at 80% identity (PDB only)
- DSSR tool required for structure annotation
- Zero-shot evaluation only (no fine-tuning on TR0)
- Evaluation: Accuracy for |i-j| ≥ 101 pairs
- **Status**: Ready for preparation

## Data Flow

### Overall Data Flow

```
Raw Dataset → Download → Extract → Filter → Normalize → Split → Encode → Validate → Output
```

### Detailed Flow for Each Task

1. **Download Phase**:
   - Input: Dataset URL, task configuration
   - Output: Raw data files (PDB, FASTA, CSV, etc.)
   - Artifacts: Download metadata, checksums

2. **Extraction Phase**:
   - Input: Raw data files
   - Output: Structured data (sequences, labels, metadata)
   - Artifacts: Extraction logs, parsing statistics

3. **Filtering Phase**:
   - Input: Structured data
   - Output: Filtered data (length, quality, redundancy)
   - Artifacts: Filter statistics, removed samples log

4. **Normalization Phase**:
   - Input: Filtered data
   - Output: Normalized data (scaled features, standardized labels)
   - Artifacts: Normalization parameters, statistics

5. **Splitting Phase**:
   - Input: Normalized data
   - Output: Train/validation/test splits
   - Artifacts: Split indices, stratification info

6. **Encoding Phase**:
   - Input: Data splits
   - Output: Model-compatible encoded data
   - Artifacts: Encoding vocabulary, label mappings

7. **Validation Phase**:
   - Input: Encoded data, validation specifications
   - Output: Validation report, quality metrics
   - Artifacts: Validation logs, discrepancy reports

## File Organization

### Directory Structure

```
downstream_tasks/
├── config/
│   ├── task_configs.yaml          # Task-specific configurations
│   ├── preprocessing_configs.yaml # Preprocessing parameters
│   └── validation_specs.yaml      # Validation specifications
├── data/
│   ├── raw/                        # Downloaded raw datasets
│   │   ├── task1_rsa/
│   │   ├── task2_subcellular/
│   │   ├── task3_contact_map/
│   │   ├── task4_rbp/
│   │   ├── task5_rri/
│   │   ├── task6_secondary_structure/
│   │   └── task7_long_range/
│   ├── preprocessed/               # Preprocessed datasets
│   │   ├── task1_rsa/
│   │   │   ├── train.json
│   │   │   ├── test.json
│   │   │   └── metadata.json
│   │   └── ...
│   └── checkpoints/                # Pipeline checkpoints
│       ├── task1_checkpoint.json
│       └── ...
├── logs/
│   ├── orchestrator.log            # Main orchestrator log
│   ├── task1_rsa.log
│   └── ...
├── reports/
│   ├── task1_rsa_report.md
│   ├── validation_reports/
│   └── consolidated_report.md
├── src/
│   ├── orchestrator.py
│   ├── task_pipeline.py
│   ├── dataset_downloader.py
│   ├── preprocessor_engine.py
│   ├── data_validator.py
│   └── utils/
│       ├── file_utils.py
│       ├── sequence_utils.py
│       └── structure_utils.py
└── tools/
    ├── freesasa/                   # External tool for RSA
    ├── dssr/                       # External tool for structure
    └── cd-hit/                     # External tool for redundancy
```

### Configuration Files

**task_configs.yaml**:
```yaml
tasks:
  task1_rsa:
    name: "RNA Solvent Accessibility"
    paper_source: "Original downstream tasks"
    dataset_url: "http://rna3dhub.org/..."
    requires_training_set: true
    output_format: "json"
    evaluation_metric: "pearson_correlation"
    
  task2_subcellular:
    name: "RNA Subcellular Localization"
    paper_source: "RNALocate"
    dataset_url: "http://www.rna-society.org/rnalocate/..."
    requires_training_set: true
    output_format: "json"
    evaluation_metric: "accuracy"
    
  # ... (remaining tasks)
```

**preprocessing_configs.yaml**:
```yaml
task1_rsa:
  sequence_length_min: null
  sequence_length_max: null
  normalization_method: "per_nucleotide"
  distance_threshold: 8.0  # Angstroms
  split_ratios: null  # Use paper-specified splits
  encoding_scheme: "token_ids"
  
task3_contact_map:
  sequence_length_min: 32
  sequence_length_max: 1000
  similarity_threshold: 0.30  # CD-HIT-EST
  distance_threshold: 8.0
  min_contacts: 5
  split_ratios: null  # Use TR221/TS80
  encoding_scheme: "token_ids"
  
# ... (remaining tasks)
```

## Error Handling

### Error Categories

1. **Download Errors**:
   - Network failures → Retry with exponential backoff
   - Invalid URLs → Log error, skip task
   - Corrupted files → Re-download, verify checksum
   - Insufficient disk space → Abort with clear message

2. **Preprocessing Errors**:
   - Invalid data format → Log error, skip sample
   - Missing required fields → Log warning, use defaults
   - Computation failures → Save checkpoint, allow resume
   - Tool execution failures → Check dependencies, retry

3. **Validation Errors**:
   - Sample count mismatch → Generate warning, continue
   - Format violations → Abort, fix preprocessing
   - Data leakage detected → Abort, fix splitting
   - Schema violations → Abort, fix encoding

### Recovery Mechanisms

1. **Checkpointing**:
   - Save state after each pipeline stage
   - Enable resume from last successful checkpoint
   - Store checkpoint metadata (timestamp, stage, outputs)

2. **Retry Logic**:
   - Network operations: 3 retries with exponential backoff
   - Tool executions: 2 retries with error logging
   - File operations: 1 retry after brief delay

3. **Graceful Degradation**:
   - Continue other tasks if one fails
   - Generate partial reports for incomplete tasks
   - Preserve all successfully completed work

## Dependencies

### External Tools

1. **FreeSASA**: Solvent accessibility calculation (Task 1)
   - Installation: `conda install -c salilab freesasa`
   - Version: ≥2.0

2. **DSSR**: RNA structure annotation (Task 6)
   - Installation: Download from x3dna.org
   - Version: ≥2.0

3. **CD-HIT-EST**: Sequence redundancy removal (Task 3)
   - Installation: `conda install -c bioconda cd-hit`
   - Version: ≥4.8

4. **BLASTclust**: Alternative redundancy removal (Task 3)
   - Installation: `conda install -c bioconda blast`
   - Version: ≥2.10

### Python Libraries

```
biopython>=1.79      # Sequence parsing and manipulation
numpy>=1.21.0        # Numerical operations
pandas>=1.3.0        # Data manipulation
requests>=2.26.0     # HTTP downloads
tqdm>=4.62.0         # Progress bars
pyyaml>=5.4.0        # Configuration parsing
jsonschema>=3.2.0    # Data validation
h5py>=3.3.0          # HDF5 file handling
scikit-learn>=0.24.0 # Data splitting and metrics
```

## Testing Strategy

### Unit Tests

1. **Dataset Downloader Tests**:
   - Test download with mock URLs
   - Test retry logic with simulated failures
   - Test checksum verification
   - Test archive extraction

2. **Preprocessor Engine Tests**:
   - Test each filtering operation
   - Test normalization methods
   - Test splitting strategies
   - Test encoding schemes

3. **Data Validator Tests**:
   - Test sample count validation
   - Test distribution checks
   - Test data leakage detection
   - Test format validation

### Integration Tests

1. **End-to-End Pipeline Tests**:
   - Test complete pipeline with small dataset
   - Verify output format correctness
   - Check checkpoint and resume functionality
   - Validate error handling

2. **Task-Specific Tests**:
   - Test each task pipeline with sample data
   - Verify paper methodology compliance
   - Check output against expected format

### Validation Tests

1. **Data Quality Tests**:
   - Compare preprocessed data statistics with paper
   - Verify no data leakage between splits
   - Check label distribution correctness

2. **Reproducibility Tests**:
   - Run pipeline multiple times, verify identical output
   - Test with different random seeds
   - Verify deterministic behavior

## Performance Considerations

### Optimization Strategies

1. **Parallel Processing**:
   - Execute independent tasks in parallel
   - Use multiprocessing for data-intensive operations
   - Batch process sequences for encoding

2. **Caching**:
   - Cache downloaded datasets
   - Cache intermediate preprocessing results
   - Cache computed features (SASA, distances)

3. **Memory Management**:
   - Stream large files instead of loading entirely
   - Process data in chunks for large datasets
   - Clear intermediate results after checkpointing

4. **Disk I/O**:
   - Use efficient file formats (HDF5, NPZ)
   - Compress large datasets
   - Minimize redundant file operations

### Expected Resource Requirements

**Per Task Estimates**:

| Task | Download Size | Preprocessed Size | Processing Time | Memory |
|------|---------------|-------------------|-----------------|--------|
| Task 1 (RSA) | ~500 MB | ~100 MB | ~2 hours | ~4 GB |
| Task 2 (SubCell) | ~50 MB | ~20 MB | ~30 min | ~2 GB |
| Task 3 (Contact) | ~2 GB | ~500 MB | ~4 hours | ~8 GB |
| Task 4 (RBP) | ~1 GB | ~300 MB | ~3 hours | ~6 GB |
| Task 5 (RRI) | ~100 MB | ~50 MB | ~1 hour | ~3 GB |
| Task 6 (SecStruct) | ~3 GB | ~1 GB | ~6 hours | ~10 GB |
| Task 7 (LongRange) | Same as Task 6 | ~200 MB | ~1 hour | ~4 GB |

**Total Requirements**:
- Disk space: ~15 GB (raw + preprocessed)
- Peak memory: ~10 GB (Task 6)
- Total processing time: ~17 hours (sequential), ~6 hours (parallel)

## Monitoring and Logging

### Logging Strategy

1. **Orchestrator Logs**:
   - Task start/completion events
   - Overall progress tracking
   - Error summaries

2. **Task-Specific Logs**:
   - Pipeline stage transitions
   - Preprocessing statistics
   - Validation results
   - Error details with stack traces

3. **Performance Logs**:
   - Execution time per stage
   - Memory usage tracking
   - Disk I/O statistics

### Progress Reporting

1. **Real-Time Progress**:
   - Progress bars for long operations
   - Stage completion notifications
   - Estimated time remaining

2. **Status Dashboard**:
   - Task completion matrix
   - Current stage for each task
   - Error count and warnings

3. **Completion Notifications**:
   - Email notification on task completion
   - Summary report generation
   - Output file locations

## Security Considerations

1. **Data Integrity**:
   - Verify checksums for all downloads
   - Validate data format before processing
   - Detect and prevent data corruption

2. **Input Validation**:
   - Sanitize file paths
   - Validate configuration parameters
   - Check for malicious file content

3. **Access Control**:
   - Restrict write access to output directories
   - Log all file operations
   - Prevent unauthorized data access

## Future Enhancements

1. **Extensibility**:
   - Plugin architecture for new tasks
   - Custom preprocessing modules
   - Configurable validation rules

2. **Scalability**:
   - Distributed processing support
   - Cloud storage integration
   - GPU acceleration for preprocessing

3. **Usability**:
   - Web-based dashboard
   - Interactive configuration wizard
   - Automated paper methodology extraction

4. **Integration**:
   - Direct model evaluation integration
   - Results visualization tools
   - Benchmark comparison framework
