#!/usr/bin/env python3
"""
Example: Integrating 2D RNA structural features into the model

This script shows how to combine sequence embeddings with 2D structural features
for improved RNA-protein interaction prediction.
"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from pathlib import Path


class RNA2DFeatureEncoder(nn.Module):
    """
    Encoder for 2D RNA structural features
    
    Takes conservation, coevolution, and accessibility maps as input
    and produces a fixed-size embedding.
    """
    
    def __init__(self, hidden_dim=128, use_cnn=True):
        super().__init__()
        self.hidden_dim = hidden_dim
        self.use_cnn = use_cnn
        
        if use_cnn:
            # CNN-based encoder for image-like features
            self.conv1 = nn.Conv1d(3, 64, kernel_size=7, padding=3)
            self.conv2 = nn.Conv1d(64, 128, kernel_size=5, padding=2)
            self.conv3 = nn.Conv1d(128, hidden_dim, kernel_size=3, padding=1)
            
            self.pool = nn.AdaptiveAvgPool1d(1)
            
        else:
            # Simple MLP encoder
            self.fc1 = nn.Linear(3, 64)
            self.fc2 = nn.Linear(64, hidden_dim)
    
    def forward(self, conservation, coevolution, accessibility):
        """
        Args:
            conservation: (batch, length) or (length,)
            coevolution: (batch, length, length) or (length, length)
            accessibility: (batch, length) or (length,)
        
        Returns:
            Structural embedding: (batch, hidden_dim) or (hidden_dim,)
        """
        # Handle single sequence
        if conservation.dim() == 1:
            conservation = conservation.unsqueeze(0)
            coevolution = coevolution.unsqueeze(0)
            accessibility = accessibility.unsqueeze(0)
            squeeze_output = True
        else:
            squeeze_output = False
        
        batch_size, length = conservation.shape
        
        if self.use_cnn:
            # Aggregate coevolution matrix to 1D
            coev_1d = coevolution.mean(dim=2)  # (batch, length)
            
            # Stack features as channels
            x = torch.stack([conservation, coev_1d, accessibility], dim=1)  # (batch, 3, length)
            
            # Apply convolutions
            x = F.relu(self.conv1(x))
            x = F.relu(self.conv2(x))
            x = F.relu(self.conv3(x))
            
            # Global pooling
            x = self.pool(x).squeeze(-1)  # (batch, hidden_dim)
            
        else:
            # Simple aggregation
            coev_1d = coevolution.mean(dim=2)
            x = torch.stack([conservation, coev_1d, accessibility], dim=-1)  # (batch, length, 3)
            
            # Average over length
            x = x.mean(dim=1)  # (batch, 3)
            
            # MLP
            x = F.relu(self.fc1(x))
            x = self.fc2(x)
        
        if squeeze_output:
            x = x.squeeze(0)
        
        return x


class CombinedRNAEncoder(nn.Module):
    """
    Combined encoder that fuses sequence and structural features
    """
    
    def __init__(self, seq_dim=640, struct_dim=128, hidden_dim=256, out_dim=128):
        super().__init__()
        
        # Sequence encoder (uses pre-computed embeddings)
        self.seq_proj = nn.Linear(seq_dim, hidden_dim)
        
        # Structural encoder
        self.struct_encoder = RNA2DFeatureEncoder(hidden_dim=struct_dim)
        self.struct_proj = nn.Linear(struct_dim, hidden_dim)
        
        # Fusion
        self.fusion = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(hidden_dim, out_dim)
        )
    
    def forward(self, seq_emb, conservation, coevolution, accessibility):
        """
        Args:
            seq_emb: Pre-computed sequence embedding (seq_dim,)
            conservation: Conservation map (length,)
            coevolution: Coevolution matrix (length, length)
            accessibility: Accessibility map (length,)
        
        Returns:
            Combined embedding (out_dim,)
        """
        # Process sequence features
        seq_feat = F.relu(self.seq_proj(seq_emb))
        
        # Process structural features
        struct_emb = self.struct_encoder(conservation, coevolution, accessibility)
        struct_feat = F.relu(self.struct_proj(struct_emb))
        
        # Fuse features
        combined = torch.cat([seq_feat, struct_feat], dim=-1)
        output = self.fusion(combined)
        
        return output


def load_features(rna_id, seq_emb_file, struct_feat_file):
    """
    Load sequence and structural features for an RNA
    
    Args:
        rna_id: RNA identifier
        seq_emb_file: Path to sequence embeddings NPZ
        struct_feat_file: Path to 2D structural features NPZ
    
    Returns:
        Dictionary with all features as torch tensors
    """
    # Load sequence embedding
    seq_data = np.load(seq_emb_file)
    if rna_id not in seq_data.files:
        raise ValueError(f"RNA {rna_id} not found in sequence embeddings")
    
    seq_emb = seq_data[rna_id]
    
    # Average over sequence length if needed (for RNA-FM embeddings)
    if seq_emb.ndim == 2:
        seq_emb = seq_emb.mean(axis=0)
    
    # Load structural features
    struct_data = np.load(struct_feat_file, allow_pickle=True)
    if rna_id not in struct_data.files:
        raise ValueError(f"RNA {rna_id} not found in structural features")
    
    struct_dict = struct_data[rna_id].item()
    
    # Convert to tensors
    features = {
        'seq_emb': torch.from_numpy(seq_emb).float(),
        'conservation': torch.from_numpy(struct_dict['conservation']).float(),
        'coevolution': torch.from_numpy(struct_dict['coevolution']).float(),
        'accessibility': torch.from_numpy(struct_dict['accessibility']).float(),
        'msa_depth': struct_dict['msa_depth']
    }
    
    return features


def example_usage():
    """Example of using the combined encoder"""
    
    # File paths
    seq_emb_file = Path("data/embeddings/rna_fm_rpi7317.npz")
    struct_feat_file = Path("data/embeddings/rna_2d_features_rpi7317.npz")
    
    # Check if files exist
    if not seq_emb_file.exists():
        print(f"Sequence embeddings not found: {seq_emb_file}")
        print("Please generate sequence embeddings first")
        return
    
    if not struct_feat_file.exists():
        print(f"Structural features not found: {struct_feat_file}")
        print("Please run RNA MSA pipeline first:")
        print("  python run_rna_msa_pipeline.py --rna_fasta data/raw/RPI7317_IncRNA_sequence.fasta")
        return
    
    # Load data
    struct_data = np.load(struct_feat_file, allow_pickle=True)
    rna_ids = struct_data.files[:5]  # Use first 5 RNAs as example
    
    print(f"Found {len(struct_data.files)} RNAs with structural features")
    print(f"Using {len(rna_ids)} RNAs for example\n")
    
    # Initialize model
    model = CombinedRNAEncoder(
        seq_dim=640,  # RNA-FM embedding dimension
        struct_dim=128,
        hidden_dim=256,
        out_dim=128
    )
    
    print("Model architecture:")
    print(model)
    print()
    
    # Process each RNA
    for rna_id in rna_ids:
        print(f"Processing {rna_id}...")
        
        try:
            # Load features
            features = load_features(rna_id, seq_emb_file, struct_feat_file)
            
            print(f"  Sequence embedding: {features['seq_emb'].shape}")
            print(f"  Conservation: {features['conservation'].shape}")
            print(f"  Coevolution: {features['coevolution'].shape}")
            print(f"  Accessibility: {features['accessibility'].shape}")
            print(f"  MSA depth: {features['msa_depth']}")
            
            # Forward pass
            with torch.no_grad():
                output = model(
                    features['seq_emb'],
                    features['conservation'],
                    features['coevolution'],
                    features['accessibility']
                )
            
            print(f"  Output embedding: {output.shape}")
            print()
            
        except Exception as e:
            print(f"  Error: {e}\n")
    
    print("Example completed successfully!")
    print("\nNext steps:")
    print("1. Integrate this encoder into your RNA-protein interaction model")
    print("2. Replace the RNA encoder in your existing model with CombinedRNAEncoder")
    print("3. Train the model and compare performance with sequence-only baseline")


def create_dataset_with_features(pairs_csv, seq_emb_file, struct_feat_file):
    """
    Create a dataset that includes both sequence and structural features
    
    Args:
        pairs_csv: CSV file with RNA-protein pairs
        seq_emb_file: Sequence embeddings NPZ
        struct_feat_file: Structural features NPZ
    
    Returns:
        Dataset object
    """
    import pandas as pd
    
    # Load pairs
    df = pd.read_csv(pairs_csv)
    
    # Load embeddings
    seq_data = np.load(seq_emb_file)
    struct_data = np.load(struct_feat_file, allow_pickle=True)
    
    dataset = []
    
    for _, row in df.iterrows():
        rna_id = row['rna_id']
        prot_id = row['prot_id']
        
        # Check if features exist
        if rna_id not in seq_data.files or rna_id not in struct_data.files:
            continue
        
        # Load features
        try:
            features = load_features(rna_id, seq_emb_file, struct_feat_file)
            features['rna_id'] = rna_id
            features['prot_id'] = prot_id
            dataset.append(features)
        except Exception as e:
            print(f"Warning: Could not load features for {rna_id}: {e}")
    
    print(f"Created dataset with {len(dataset)} RNA-protein pairs")
    
    return dataset


if __name__ == "__main__":
    print("=" * 80)
    print("Example: Integrating 2D RNA Structural Features")
    print("=" * 80)
    print()
    
    example_usage()
