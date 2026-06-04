# Requirements Document

## Introduction

This document specifies requirements for implementing a comprehensive evaluation system for a pretrained RNA structure prediction model across seven downstream tasks. The system must download datasets from published papers (ERNIE-RNA, RNAErnie, RNA-km), preprocess them according to the original methodologies, and prepare them for either zero-shot evaluation or fine-tuning. Each task must be evaluated using the exact same datasets and preprocessing steps as reported in the corresponding papers to ensure fair comparison with baseline models.

## Glossary

- **Evaluation_System**: The complete software system that orchestrates dataset acquisition, preprocessing, and evaluation preparation for all seven downstream tasks
- **Dataset_Downloader**: Component responsible for retrieving datasets from paper repositories or supplementary materials
- **Preprocessor**: Component that transforms raw datasets into the format required by the pretrained model
- **Task_Pipeline**: The complete workflow for a single downstream task from dataset acquisition through evaluation preparation
- **Zero_Shot_Evaluation**: Model evaluation without fine-tuning on task-specific training data
- **Fine_Tuning**: Adapting the pretrained model using task-specific training data before evaluation
- **Test_Set**: The held-out dataset used for final model evaluation
- **Train_Set**: The dataset used for fine-tuning the model (when required)
- **Paper_Methodology**: The exact preprocessing steps, data splits, and evaluation protocols described in the original research paper
- **RSA**: RNA Solvent Accessibility - prediction of nucleotide exposure to solvent
- **Contact_Map**: Binary matrix indicating spatial proximity between RNA nucleotide pairs
- **RBP**: RNA-Binding Protein interaction prediction
- **RRI**: RNA-RNA Interaction prediction
- **Secondary_Structure**: RNA base-pairing patterns (stems, loops, bulges)
- **Long_Range_Interactions**: RNA structural contacts between distant nucleotides in sequence
- **Subcellular_Localization**: Prediction of RNA location within the cell (nucleus, cytoplasm, etc.)

## Requirements

### Requirement 1: Dataset Acquisition

**User Story:** As a researcher, I want to automatically download all required datasets for the seven downstream tasks, so that I can evaluate my pretrained RNA model without manual data collection.

#### Acceptance Criteria

1. WHEN the Dataset_Downloader is invoked for Task 1 (RSA prediction), THE Dataset_Downloader SHALL retrieve the RSA dataset with test set and training set if available
2. WHEN the Dataset_Downloader is invoked for Task 2 (subcellular localization), THE Dataset_Downloader SHALL retrieve the subcellular localization dataset with test set and training set if available
3. WHEN the Dataset_Downloader is invoked for Task 3 (contact map prediction), THE Dataset_Downloader SHALL retrieve the contact map dataset from the ERNIE-RNA paper repository with test set and training set if available
4. WHEN the Dataset_Downloader is invoked for Task 4 (RNA-protein binding), THE Dataset_Downloader SHALL retrieve the RBP dataset from the ERNIE-RNA paper repository with test set and training set if available
5. WHEN the Dataset_Downloader is invoked for Task 5 (RNA-RNA interaction), THE Dataset_Downloader SHALL retrieve the RRI dataset from the RNAErnie paper repository with test set and training set if available
6. WHEN the Dataset_Downloader is invoked for Task 6 (secondary structure), THE Dataset_Downloader SHALL retrieve the secondary structure dataset from the RNAErnie paper repository with test set and training set if available
7. WHEN the Dataset_Downloader is invoked for Task 7 (long-range interactions), THE Dataset_Downloader SHALL retrieve the long-range interactions dataset from the RNA-km paper section 3.4 with test set and training set if available
8. WHEN a dataset download fails, THE Dataset_Downloader SHALL return an error message indicating the task number, dataset name, and failure reason
9. WHEN a dataset is successfully downloaded, THE Dataset_Downloader SHALL verify the downloaded files are not corrupted and contain expected data fields
10. FOR ALL downloaded datasets, THE Dataset_Downloader SHALL store metadata including source paper, download date, dataset version, and file checksums

### Requirement 2: Paper Methodology Extraction

**User Story:** As a researcher, I want to extract the exact preprocessing steps and evaluation protocols from each paper, so that I can replicate the baseline results fairly.

#### Acceptance Criteria

1. WHEN the Evaluation_System processes Task 3 (contact map), THE Evaluation_System SHALL extract preprocessing methodology from the ERNIE-RNA paper section describing contact map prediction
2. WHEN the Evaluation_System processes Task 4 (RBP), THE Evaluation_System SHALL extract preprocessing methodology from the ERNIE-RNA paper section describing RNA-protein binding prediction
3. WHEN the Evaluation_System processes Task 5 (RRI), THE Evaluation_System SHALL extract preprocessing methodology from the RNAErnie paper section describing RNA-RNA interaction
4. WHEN the Evaluation_System processes Task 6 (secondary structure), THE Evaluation_System SHALL extract preprocessing methodology from the RNAErnie paper section describing secondary structure prediction
5. WHEN the Evaluation_System processes Task 7 (long-range interactions), THE Evaluation_System SHALL extract preprocessing methodology from the RNA-km paper section 3.4
6. FOR ALL tasks, THE Evaluation_System SHALL document the data split ratios, sequence length limits, filtering criteria, and normalization procedures used in the Paper_Methodology
7. FOR ALL tasks, THE Evaluation_System SHALL document the evaluation metrics reported in the original paper

### Requirement 3: Dataset Preprocessing

**User Story:** As a researcher, I want to preprocess each dataset exactly as described in the corresponding paper, so that my evaluation results are directly comparable to published baselines.

#### Acceptance Criteria

1. WHEN the Preprocessor receives a raw dataset and Paper_Methodology, THE Preprocessor SHALL apply all filtering criteria specified in the Paper_Methodology
2. WHEN the Preprocessor applies sequence length limits, THE Preprocessor SHALL exclude sequences outside the specified length range
3. WHEN the Preprocessor applies normalization, THE Preprocessor SHALL use the same normalization parameters as the Paper_Methodology
4. WHEN the Preprocessor creates data splits, THE Preprocessor SHALL use the same train-validation-test split ratios as the Paper_Methodology
5. IF the Paper_Methodology specifies sequence encoding, THEN THE Preprocessor SHALL encode sequences using the specified encoding scheme
6. IF the Paper_Methodology specifies data augmentation, THEN THE Preprocessor SHALL apply the same augmentation techniques
7. WHEN preprocessing completes, THE Preprocessor SHALL output preprocessed data in a format compatible with the pretrained model input requirements
8. FOR ALL preprocessing operations, THE Preprocessor SHALL log each transformation step with input and output statistics
9. WHEN preprocessing fails, THE Preprocessor SHALL return an error message indicating the failed step and the reason

### Requirement 4: Evaluation Data Preparation

**User Story:** As a researcher, I want to prepare preprocessed datasets for both zero-shot evaluation and fine-tuning, so that I can assess my model's performance under different evaluation scenarios.

#### Acceptance Criteria

1. WHEN the Task_Pipeline prepares data for zero-shot evaluation, THE Task_Pipeline SHALL output only the Test_Set in model-compatible format
2. WHEN the Task_Pipeline prepares data for fine-tuning evaluation, THE Task_Pipeline SHALL output both Train_Set and Test_Set in model-compatible format
3. FOR ALL prepared datasets, THE Task_Pipeline SHALL include metadata files documenting the preprocessing steps applied
4. FOR ALL prepared datasets, THE Task_Pipeline SHALL include a validation script that verifies data format correctness
5. WHEN data preparation completes, THE Task_Pipeline SHALL generate a summary report containing dataset statistics including sample counts, sequence length distributions, and class distributions
6. IF the pretrained model requires specific input tokenization, THEN THE Task_Pipeline SHALL apply the model-specific tokenization to all sequences
7. IF the pretrained model requires specific label encoding, THEN THE Task_Pipeline SHALL encode labels according to the model requirements

### Requirement 5: Task Pipeline Orchestration

**User Story:** As a researcher, I want to execute all seven task pipelines independently and in parallel, so that I can efficiently prepare all evaluation datasets without sequential dependencies.

#### Acceptance Criteria

1. WHEN the Evaluation_System is invoked, THE Evaluation_System SHALL execute all seven Task_Pipeline instances independently
2. WHEN a Task_Pipeline completes successfully, THE Evaluation_System SHALL mark that task as completed and notify the designated recipient
3. WHEN a Task_Pipeline fails, THE Evaluation_System SHALL log the failure without blocking other Task_Pipeline executions
4. FOR ALL Task_Pipeline executions, THE Evaluation_System SHALL track progress including current step, elapsed time, and completion percentage
5. WHEN all Task_Pipeline instances complete, THE Evaluation_System SHALL generate a consolidated report summarizing the status of all seven tasks
6. THE Evaluation_System SHALL allow individual Task_Pipeline instances to be re-executed without affecting completed tasks

### Requirement 6: Task 1 - RSA Prediction Pipeline

**User Story:** As a researcher, I want to prepare the RSA prediction dataset from RNAsnap2 benchmark, so that I can evaluate my model's ability to predict RNA solvent accessibility.

#### Acceptance Criteria

1. WHEN the Task_Pipeline for RSA prediction executes, THE Task_Pipeline SHALL download the RNAsnap2 benchmark dataset from github.com/jaswindersingh2/RNAsnap2
2. WHEN the RSA dataset is downloaded, THE Task_Pipeline SHALL extract TR89 training split containing 89 RNA chains
3. WHEN the RSA dataset is downloaded, THE Task_Pipeline SHALL extract CN48 test split containing 48 protein-free RNA chains
4. WHEN the RSA dataset is downloaded, THE Task_Pipeline SHALL extract TS44 test split containing 44 protein-bound RNA chains
5. WHEN the Preprocessor processes RSA data, THE Preprocessor SHALL parse FASTA sequences paired with position-wise float labels for normalized SASA values
6. WHEN preprocessing completes, THE Task_Pipeline SHALL output RSA data in model-compatible format with per-nucleotide normalized SASA labels
7. WHEN the RSA Task_Pipeline completes, THE Evaluation_System SHALL send the completed task notification to the designated recipient

### Requirement 7: Task 2 - Subcellular Localization Pipeline

**User Story:** As a researcher, I want to prepare the DM3Loc/LOCAS-style mRNA subcellular localization benchmark (RNALocate-derived), so that I can evaluate my model on a runnable multi-label localization task even when raw mRNAs exceed typical context limits.

#### Acceptance Criteria

1. WHEN the Task_Pipeline for subcellular localization executes, THE Task_Pipeline SHALL use the DM3Loc-provided benchmark FASTA with 7-bit label prefixes in each header
2. WHEN preprocessing runs, THE Task_Pipeline SHALL map the 7-bit labels to 6 prediction targets (excluding cytoplasm as an auxiliary class used in DM3Loc training)
3. WHEN preprocessing runs, THE Task_Pipeline SHALL drop records with no positives among the 6 prediction targets by default
4. WHEN an mRNA sequence length exceeds 8000 nt, THE Task_Pipeline SHALL apply DM3Loc-style end-truncation (first 4000 nt + last 4000 nt) so max input length is 8000
5. WHEN preprocessing completes, THE Task_Pipeline SHALL output gzipped JSONL records plus deterministic 5-fold CV split manifest (seeded)
6. WHEN the Task_Pipeline completes, THE Task_Pipeline SHALL write a summary + validation report (counts + max post-truncation length) to support reproducibility

### Requirement 8: Task 3 - Contact Map Prediction Pipeline

**User Story:** As a researcher, I want to prepare the RNAcontact dataset from the ERNIE-RNA paper, so that I can evaluate my model's ability to predict RNA spatial contacts.

#### Acceptance Criteria

1. WHEN the Task_Pipeline for contact map prediction executes, THE Task_Pipeline SHALL download the RNAcontact dataset from Sun et al., Bioinformatics 2021
2. WHEN the contact map dataset is downloaded, THE Task_Pipeline SHALL extract the 301 filtered sequences from the initial 1,786 RNA sequences from PDB
3. WHEN the contact map dataset is processed, THE Task_Pipeline SHALL create TR221 training split, TR168 and VL43 validation splits, and TS80 test split
4. WHEN the Preprocessor processes contact map data, THE Preprocessor SHALL apply contact definition using minimum atom-pair distance below 8Å
5. WHEN preprocessing completes, THE Task_Pipeline SHALL output contact map data in model-compatible format with binary contact matrices
6. WHEN the contact map Task_Pipeline completes, THE Evaluation_System SHALL send the completed task notification to the designated recipient

### Requirement 9: Task 4 - RNA-Protein Binding Pipeline

**User Story:** As a researcher, I want to prepare the PrismNet dataset from the ERNIE-RNA paper, so that I can evaluate my model's ability to predict RBP binding sites.

#### Acceptance Criteria

1. WHEN the Task_Pipeline for RBP prediction executes, THE Task_Pipeline SHALL download the PrismNet dataset from github.com/kuixu/PrismNet
2. WHEN the RBP dataset is downloaded, THE Task_Pipeline SHALL extract data for 17 RBPs in HeLa cell environment
3. WHEN the RBP dataset is processed, THE Task_Pipeline SHALL extract exactly 101 nt sequences with 1,827 to 15,002 sequences per RBP
4. WHEN the RBP dataset is processed, THE Task_Pipeline SHALL parse HDF5 (.h5) file format
5. WHEN the Preprocessor creates data splits, THE Preprocessor SHALL use the exact same train/val/test split as the ERNIE-RNA paper for fair comparison
6. WHEN preprocessing completes, THE Task_Pipeline SHALL output RBP data in model-compatible format with binding site labels
7. WHEN the RBP Task_Pipeline completes, THE Evaluation_System SHALL send the completed task notification to the designated recipient

### Requirement 10: Task 5 - RNA-RNA Interaction Pipeline

**User Story:** As a researcher, I want to prepare the DeepMirTar dataset from the RNAErnie paper, so that I can evaluate my model's ability to predict RNA-RNA interactions.

#### Acceptance Criteria

1. WHEN the Task_Pipeline for RRI prediction executes, THE Task_Pipeline SHALL download the DeepMirTar dataset from Wen et al., Bioinformatics 2018
2. WHEN the RRI dataset is downloaded, THE Task_Pipeline SHALL extract 27,720 total interaction pairs with 13,860 positive and 13,860 negative samples
3. WHEN the RRI dataset is processed, THE Task_Pipeline SHALL filter for miRNA sequences less than 26 nt and mRNA site sequences less than 53 nt
4. WHEN the RRI dataset is processed, THE Task_Pipeline SHALL parse tab-separated interaction records with labels
5. WHEN the Preprocessor creates data splits, THE Preprocessor SHALL use the exact same train/val/test split as the RNAErnie paper for fair comparison
6. WHEN preprocessing completes, THE Task_Pipeline SHALL output RRI data in model-compatible format with interaction labels
7. WHEN the RRI Task_Pipeline completes, THE Evaluation_System SHALL send the completed task notification to the designated recipient

### Requirement 11: Task 6 - Secondary Structure Prediction Pipeline (ON HOLD)

**User Story:** As a researcher, I want to evaluate whether the secondary structure prediction task is meaningful for my pretrained model, given that it was already trained on secondary structure data.

#### Acceptance Criteria

1. THIS TASK IS ON HOLD pending decision on whether to include it in the downstream task list
2. IF the task is approved, THEN THE Task_Pipeline SHALL download datasets from Danaee et al., Nucleic Acids Research 2018 with TR0, TV0, and TS0 splits
3. IF the task is approved, THEN THE Task_Pipeline SHALL download RNAStralign from Tan et al., Nucleic Acids Research 2017
4. IF the task is approved, THEN THE Task_Pipeline SHALL download ArchiveII from Sloma and Mathews, RNA 2016
5. IF the task is approved, THEN THE Task_Pipeline SHALL download Szikszai 2022 Cross-Family Split from Szikszai et al., Bioinformatics 2022
6. THE Evaluation_System SHALL NOT proceed with this task until explicit approval is received

### Requirement 12: Task 7 - Long-Range Base Pair Prediction Pipeline

**User Story:** As a researcher, I want to prepare the long-range base pair prediction datasets, so that I can evaluate my model's ability to capture distant structural contacts.

#### Acceptance Criteria

1. WHEN the Task_Pipeline for long-range base pair prediction executes, THE Task_Pipeline SHALL download Dataset A: PDB Benchmark with 532 non-redundant RNA sequences from PDB
2. WHEN Dataset A is processed, THE Task_Pipeline SHALL ensure protein and DNA binding partners are removed
3. WHEN Dataset A is processed, THE Task_Pipeline SHALL apply redundancy filtering using CD-HIT-EST at 80 percent identity
4. WHEN Dataset A is processed, THE Task_Pipeline SHALL annotate secondary structures using DSSR tool
5. WHEN Dataset A is processed, THE Task_Pipeline SHALL identify the 158 sequences that contain pseudoknots
6. WHEN the Task_Pipeline executes, THE Task_Pipeline SHALL download Dataset B: bpRNA-1m TR0 split containing 10,814 sequences
7. WHEN Dataset B is processed, THE Task_Pipeline SHALL identify the 999 sequences that contain pseudoknots
8. WHEN Dataset B is processed, THE Task_Pipeline SHALL prepare data exclusively for zero-shot evaluation without fine-tuning
9. WHEN preprocessing completes, THE Task_Pipeline SHALL output long-range base pair data in model-compatible format with long-range contact labels
10. WHEN the long-range base pair Task_Pipeline completes, THE Evaluation_System SHALL send the completed task notification to the designated recipient

### Requirement 13: Data Validation and Quality Assurance

**User Story:** As a researcher, I want to validate that all preprocessed datasets match the paper specifications, so that I can trust my evaluation results are fair comparisons.

#### Acceptance Criteria

1. WHEN a Task_Pipeline completes preprocessing, THE Evaluation_System SHALL validate that the Test_Set sample count matches the paper-reported test set size within 5 percent tolerance
2. WHEN a Task_Pipeline completes preprocessing, THE Evaluation_System SHALL validate that sequence length distributions match the Paper_Methodology specifications
3. WHEN a Task_Pipeline completes preprocessing, THE Evaluation_System SHALL validate that class distributions match the Paper_Methodology specifications within 5 percent tolerance
4. IF validation fails for any dataset, THEN THE Evaluation_System SHALL generate a validation report highlighting discrepancies
5. FOR ALL preprocessed datasets, THE Evaluation_System SHALL verify that no data leakage exists between Train_Set and Test_Set
6. FOR ALL preprocessed datasets, THE Evaluation_System SHALL verify that all required data fields are present and correctly formatted

### Requirement 14: Reproducibility and Documentation

**User Story:** As a researcher, I want comprehensive documentation of all preprocessing steps and data transformations, so that I can reproduce the evaluation and debug any issues.

#### Acceptance Criteria

1. FOR ALL Task_Pipeline executions, THE Evaluation_System SHALL generate a preprocessing log documenting every transformation applied with timestamps
2. FOR ALL downloaded datasets, THE Evaluation_System SHALL store the original paper citation, dataset URL, and download timestamp
3. FOR ALL preprocessed datasets, THE Evaluation_System SHALL generate a README file describing the data format, label encoding, and usage instructions
4. WHEN a Task_Pipeline completes, THE Evaluation_System SHALL generate a reproducibility report containing all parameters, random seeds, and software versions used
5. THE Evaluation_System SHALL provide a configuration file format that allows users to specify custom preprocessing parameters while maintaining Paper_Methodology defaults
6. FOR ALL evaluation runs, THE Evaluation_System SHALL log the pretrained model version, evaluation date, and hardware specifications

### Requirement 15: Error Handling and Recovery

**User Story:** As a researcher, I want robust error handling and recovery mechanisms, so that transient failures do not require restarting the entire evaluation pipeline.

#### Acceptance Criteria

1. WHEN a dataset download fails due to network issues, THE Dataset_Downloader SHALL retry the download up to 3 times with exponential backoff
2. WHEN a preprocessing step fails, THE Preprocessor SHALL save intermediate results and allow resumption from the failed step
3. WHEN a Task_Pipeline fails, THE Evaluation_System SHALL preserve all successfully completed steps and allow re-execution of only the failed steps
4. IF disk space is insufficient, THEN THE Evaluation_System SHALL return an error message indicating required space and available space
5. IF a required dependency is missing, THEN THE Evaluation_System SHALL return an error message listing all missing dependencies
6. WHEN an error occurs, THE Evaluation_System SHALL log the full error traceback and context information for debugging
