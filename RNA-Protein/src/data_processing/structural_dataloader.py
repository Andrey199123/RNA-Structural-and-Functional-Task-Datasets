#!/usr/bin/env python3
"""
Data loader that combines sequence embeddings with structural features.

Loads:
- RNA-FM sequence embeddings
- ESM2 protein sequence embeddings  
- RNA 2D structural features (distance maps, contact maps)
- Protein 2D structural features (distance maps, contact maps)
- RNA base interaction features (base pairs, stacking)

Provides unified interface for training with both sequence and structure.
"""

import numpy as np
import torch
from pathlib import Path
from typing import Dict, Tuple, Optional, List
import json


class StructuralFeatureLoader:
    """Loads and manages structural features for RNA-protein pairs."""
    
    def __init__(
        self,
        rna_2d_dir: Path,
        protein_2d_dir: Path,
        rna_interactions_dir: Path,
        max_rna_len: int = 800,
        max_prot_len: int = 1000,
    ):
        self.rna_2d_dir = Path(rna_2d_dir)
        self.protein_2d_dir = Path(protein_2d_dir)
        self.rna_interactions_dir = Path(rna_interactions_dir)
        self.max_rna_len = max_rna_len
        self.max_prot_len = max_prot_len
        
        # Cache for loaded features
        self._rna_2d_cache = {}
        self._prot_2d_cache = {}
        self._rna_interact_cache = {}
    
    def _pad_or_crop_matrix(self, matrix: np.ndarray, target_size: int) -> np.ndarray:
        """Pad or crop a square matrix to target size."""
        n = matrix.shape[0]
        if n == target_size:
            return matrix
        elif n < target_size:
            # Pad with zeros
            padded = np.zeros((target_size, target_size), dtype=matrix.dtype)
            padded[:n, :n] = matrix
            return padded
        else:
            # Crop (center crop)
            start = (n - target_size) // 2
            return matrix[start:start+target_size, start:start+target_size]

    
    def load_rna_2d_features(self, rna_id: str) -> Optional[Dict[str, np.ndarray]]:
        """Load 2D structural features for an RNA."""
        if rna_id in self._rna_2d_cache:
            return self._rna_2d_cache[rna_id]
        
        npz_file = self.rna_2d_dir / f"{rna_id}_2d.npz"
        if not npz_file.exists():
            return None
        
        try:
            data = np.load(npz_file)
            features = {
                'distance_matrix': data['distance_matrix'],
                'contact_map': data['contact_map'],
                'normalized_distance': data['normalized_distance'],
                'sequence_length': int(data['sequence_length']),
            }
            self._rna_2d_cache[rna_id] = features
            return features
        except Exception:
            return None
    
    def load_protein_2d_features(self, prot_id: str) -> Optional[Dict[str, np.ndarray]]:
        """Load 2D structural features for a protein."""
        if prot_id in self._prot_2d_cache:
            return self._prot_2d_cache[prot_id]
        
        npz_file = self.protein_2d_dir / f"{prot_id}_2d.npz"
        if not npz_file.exists():
            return None
        
        try:
            data = np.load(npz_file)
            features = {
                'distance_matrix': data['distance_matrix'],
                'contact_map': data['contact_map'],
                'normalized_distance': data['normalized_distance'],
                'sequence_length': int(data['sequence_length']),
            }
            self._prot_2d_cache[prot_id] = features
            return features
        except Exception:
            return None
    
    def load_rna_interactions(self, rna_id: str) -> Optional[Dict[str, np.ndarray]]:
        """Load base pairing and stacking features for an RNA."""
        if rna_id in self._rna_interact_cache:
            return self._rna_interact_cache[rna_id]
        
        npz_file = self.rna_interactions_dir / f"{rna_id}_interactions.npz"
        if not npz_file.exists():
            return None
        
        try:
            data = np.load(npz_file)
            features = {
                'base_pair_matrix': data['base_pair_matrix'],
                'stacking_matrix': data['stacking_matrix'],
                'sequence_length': int(data['sequence_length']),
            }
            self._rna_interact_cache[rna_id] = features
            return features
        except Exception:
            return None

    
    def get_rna_structural_tensor(
        self, 
        rna_id: str, 
        target_size: Optional[int] = None
    ) -> Optional[torch.Tensor]:
        """
        Get combined RNA structural features as a tensor.
        
        Returns tensor of shape (C, H, W) where C=4 channels:
        - Channel 0: normalized distance map
        - Channel 1: contact map
        - Channel 2: base pair matrix
        - Channel 3: stacking matrix
        """
        target_size = target_size or self.max_rna_len
        
        feat_2d = self.load_rna_2d_features(rna_id)
        feat_interact = self.load_rna_interactions(rna_id)
        
        if feat_2d is None:
            return None
        
        # Get matrices
        norm_dist = self._pad_or_crop_matrix(feat_2d['normalized_distance'], target_size)
        contact = self._pad_or_crop_matrix(feat_2d['contact_map'], target_size)
        
        if feat_interact is not None:
            base_pair = self._pad_or_crop_matrix(feat_interact['base_pair_matrix'], target_size)
            stacking = self._pad_or_crop_matrix(feat_interact['stacking_matrix'], target_size)
        else:
            base_pair = np.zeros((target_size, target_size), dtype=np.float32)
            stacking = np.zeros((target_size, target_size), dtype=np.float32)
        
        # Stack into channels and ensure float32
        tensor = np.stack([norm_dist, contact, base_pair, stacking], axis=0).astype(np.float32)
        return torch.from_numpy(tensor)
    
    def get_protein_structural_tensor(
        self, 
        prot_id: str, 
        target_size: Optional[int] = None
    ) -> Optional[torch.Tensor]:
        """
        Get protein structural features as a tensor.
        
        Returns tensor of shape (C, H, W) where C=2 channels:
        - Channel 0: normalized distance map
        - Channel 1: contact map
        """
        target_size = target_size or self.max_prot_len
        
        feat_2d = self.load_protein_2d_features(prot_id)
        
        if feat_2d is None:
            return None
        
        norm_dist = self._pad_or_crop_matrix(feat_2d['normalized_distance'], target_size)
        contact = self._pad_or_crop_matrix(feat_2d['contact_map'], target_size)
        
        # Stack into channels and ensure float32
        tensor = np.stack([norm_dist, contact], axis=0).astype(np.float32)
        return torch.from_numpy(tensor)



class RPIDataset(torch.utils.data.Dataset):
    """
    Dataset for RNA-protein interaction prediction with structural features.
    
    Combines:
    - Sequence embeddings (RNA-FM, ESM2)
    - 2D structural features (distance/contact maps)
    - Base interaction features (pairing, stacking)
    """
    
    def __init__(
        self,
        pairs_csv: Path,
        rna_embeddings: Dict[str, np.ndarray],
        prot_embeddings: Dict[str, np.ndarray],
        structural_loader: Optional[StructuralFeatureLoader] = None,
        use_structural: bool = True,
        rna_struct_size: int = 256,  # Resize structural maps for efficiency
        prot_struct_size: int = 256,
    ):
        import pandas as pd
        
        self.pairs_df = pd.read_csv(pairs_csv, encoding='utf-8-sig')
        self.rna_embeddings = rna_embeddings
        self.prot_embeddings = prot_embeddings
        self.structural_loader = structural_loader
        self.use_structural = use_structural and structural_loader is not None
        self.rna_struct_size = rna_struct_size
        self.prot_struct_size = prot_struct_size
        
        # Filter to pairs where we have embeddings
        valid_pairs = []
        for idx, row in self.pairs_df.iterrows():
            rna_id = row['RNA names']
            prot_id = row['Protein names']
            if rna_id in rna_embeddings and prot_id in prot_embeddings:
                valid_pairs.append(idx)
        
        self.valid_indices = valid_pairs
        print(f"Dataset: {len(self.valid_indices)} valid pairs out of {len(self.pairs_df)}")
    
    def __len__(self):
        return len(self.valid_indices)
    
    def __getitem__(self, idx):
        row = self.pairs_df.iloc[self.valid_indices[idx]]
        rna_id = row['RNA names']
        prot_id = row['Protein names']
        label = float(row['Labels'])
        
        # Sequence embeddings
        rna_emb = torch.from_numpy(self.rna_embeddings[rna_id].astype(np.float32))
        prot_emb = torch.from_numpy(self.prot_embeddings[prot_id].astype(np.float32))
        
        sample = {
            'rna_id': rna_id,
            'prot_id': prot_id,
            'rna_embedding': rna_emb,
            'prot_embedding': prot_emb,
            'label': torch.tensor(label, dtype=torch.float32),
        }
        
        # Add structural features if available
        if self.use_structural:
            rna_struct = self.structural_loader.get_rna_structural_tensor(
                rna_id, self.rna_struct_size
            )
            prot_struct = self.structural_loader.get_protein_structural_tensor(
                prot_id, self.prot_struct_size
            )
            
            if rna_struct is not None:
                sample['rna_structural'] = rna_struct
            if prot_struct is not None:
                sample['prot_structural'] = prot_struct
        
        return sample



def load_embeddings_as_dict(npz_path: Path) -> Dict[str, np.ndarray]:
    """Load embeddings NPZ file as a dictionary.
    
    The NPZ files use IDs as keys directly (e.g., npz['BCAR4'] = embedding).
    """
    data = np.load(npz_path, allow_pickle=True)
    return {key: data[key] for key in data.files}


def create_dataloaders(
    pairs_csv: Path,
    rna_npz: Path,
    prot_npz: Path,
    structural_features_dir: Path,
    batch_size: int = 32,
    use_structural: bool = True,
    rna_struct_size: int = 256,
    prot_struct_size: int = 256,
    num_workers: int = 0,
) -> torch.utils.data.DataLoader:
    """
    Create a DataLoader with both sequence and structural features.
    """
    # Load sequence embeddings
    rna_emb = load_embeddings_as_dict(rna_npz)
    prot_emb = load_embeddings_as_dict(prot_npz)
    
    print(f"Loaded {len(rna_emb)} RNA embeddings")
    print(f"Loaded {len(prot_emb)} protein embeddings")
    
    # Create structural loader
    structural_loader = None
    if use_structural:
        structural_loader = StructuralFeatureLoader(
            rna_2d_dir=structural_features_dir / 'rna_2d_features',
            protein_2d_dir=structural_features_dir / 'protein_2d_features',
            rna_interactions_dir=structural_features_dir / 'rna_interactions',
        )
    
    # Create dataset
    dataset = RPIDataset(
        pairs_csv=pairs_csv,
        rna_embeddings=rna_emb,
        prot_embeddings=prot_emb,
        structural_loader=structural_loader,
        use_structural=use_structural,
        rna_struct_size=rna_struct_size,
        prot_struct_size=prot_struct_size,
    )
    
    # Create dataloader
    dataloader = torch.utils.data.DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=num_workers,
        collate_fn=collate_rpi_batch,
    )
    
    return dataloader


def collate_rpi_batch(batch: List[Dict]) -> Dict:
    """Custom collate function to handle variable presence of structural features."""
    result = {
        'rna_ids': [s['rna_id'] for s in batch],
        'prot_ids': [s['prot_id'] for s in batch],
        'rna_embeddings': torch.stack([s['rna_embedding'] for s in batch]),
        'prot_embeddings': torch.stack([s['prot_embedding'] for s in batch]),
        'labels': torch.stack([s['label'] for s in batch]),
    }
    
    # Handle structural features (may not be present for all samples)
    # Find a sample with structural features to get the shape
    rna_template = None
    prot_template = None
    for s in batch:
        if 'rna_structural' in s and rna_template is None:
            rna_template = s['rna_structural']
        if 'prot_structural' in s and prot_template is None:
            prot_template = s['prot_structural']
        if rna_template is not None and prot_template is not None:
            break
    
    # Stack RNA structural features
    if rna_template is not None:
        rna_structs = []
        for s in batch:
            if 'rna_structural' in s:
                rna_structs.append(s['rna_structural'].float())
            else:
                # Create zero tensor as placeholder with same dtype
                rna_structs.append(torch.zeros_like(rna_template).float())
        result['rna_structural'] = torch.stack(rna_structs)
    
    # Stack protein structural features
    if prot_template is not None:
        prot_structs = []
        for s in batch:
            if 'prot_structural' in s:
                prot_structs.append(s['prot_structural'].float())
            else:
                prot_structs.append(torch.zeros_like(prot_template).float())
        result['prot_structural'] = torch.stack(prot_structs)
    
    return result


if __name__ == "__main__":
    # Test the dataloader
    base_dir = Path(__file__).parent.parent.parent / 'data'
    
    loader = create_dataloaders(
        pairs_csv=base_dir / 'RPI7317/rpi7317_pairs_balanced.csv',
        rna_npz=base_dir / 'embeddings/rna_fm_rpi7317.npz',
        prot_npz=base_dir / 'embeddings/esm2_35m_rpi7317.npz',
        structural_features_dir=base_dir / 'structural_features',
        batch_size=4,
        use_structural=True,
    )
    
    # Test one batch
    batch = next(iter(loader))
    print("\nBatch contents:")
    for k, v in batch.items():
        if isinstance(v, torch.Tensor):
            print(f"  {k}: {v.shape}")
        else:
            print(f"  {k}: {type(v).__name__} len={len(v)}")
