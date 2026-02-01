#!/usr/bin/env python3
"""
RNA-Protein Interaction model with both sequence and structural features.

Combines:
- Sequence embeddings (RNA-FM, ESM2) processed through MLP
- 2D structural features (distance/contact maps) processed through CNN
- Fusion of both modalities for final prediction
"""

import torch
import torch.nn as nn
import torch.nn.functional as F

try:
    from .structural_encoder import StructuralCNN
except ImportError:
    from structural_encoder import StructuralCNN


class RPIStructural(nn.Module):
    """
    RNA-Protein Interaction prediction model with structural features.
    
    Architecture:
    1. Sequence branch: MLP to project embeddings
    2. Structure branch: CNN to encode 2D maps
    3. Fusion: Concatenate and predict interaction
    """
    
    def __init__(
        self,
        rna_seq_dim: int = 640,      # RNA-FM embedding dim
        prot_seq_dim: int = 480,     # ESM2 embedding dim
        rna_struct_channels: int = 4, # distance, contact, base_pair, stacking
        prot_struct_channels: int = 2, # distance, contact
        hidden_dim: int = 256,
        struct_dim: int = 128,        # Output dim of structural encoder
        dropout: float = 0.2,
        use_structural: bool = True,
    ):
        super().__init__()
        
        self.use_structural = use_structural
        
        # Sequence encoders (MLP)
        self.rna_seq_encoder = nn.Sequential(
            nn.Linear(rna_seq_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
        )
        
        self.prot_seq_encoder = nn.Sequential(
            nn.Linear(prot_seq_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
        )

        
        # Structural encoders (CNN)
        if use_structural:
            self.rna_struct_encoder = StructuralCNN(
                in_channels=rna_struct_channels,
                out_dim=struct_dim,
                hidden_channels=[32, 64, 128],
                dropout=dropout,
            )
            
            self.prot_struct_encoder = StructuralCNN(
                in_channels=prot_struct_channels,
                out_dim=struct_dim,
                hidden_channels=[32, 64, 128],
                dropout=dropout,
            )
        
        # Fusion and prediction head
        if use_structural:
            fusion_dim = (hidden_dim + struct_dim) * 2  # RNA + Protein, each with seq + struct
        else:
            fusion_dim = hidden_dim * 2  # RNA + Protein, seq only
        
        self.fusion = nn.Sequential(
            nn.Linear(fusion_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
        )
        
        self.classifier = nn.Linear(hidden_dim // 2, 1)
    
    def encode_rna(self, seq_emb, struct_feat=None):
        """Encode RNA using sequence and optionally structural features."""
        seq_out = self.rna_seq_encoder(seq_emb)
        
        if self.use_structural and struct_feat is not None:
            # Check if struct_feat is all zeros (placeholder for missing features)
            if struct_feat.abs().sum() > 0:
                struct_out = self.rna_struct_encoder(struct_feat)
                return torch.cat([seq_out, struct_out], dim=-1)
            else:
                # Return zeros for structural component
                batch_size = seq_emb.shape[0]
                device = seq_emb.device
                struct_out = torch.zeros(batch_size, self.rna_struct_encoder.fc.out_features, device=device)
                return torch.cat([seq_out, struct_out], dim=-1)
        return seq_out
    
    def encode_protein(self, seq_emb, struct_feat=None):
        """Encode protein using sequence and optionally structural features."""
        seq_out = self.prot_seq_encoder(seq_emb)
        
        if self.use_structural and struct_feat is not None:
            # Check if struct_feat is all zeros (placeholder for missing features)
            if struct_feat.abs().sum() > 0:
                struct_out = self.prot_struct_encoder(struct_feat)
                return torch.cat([seq_out, struct_out], dim=-1)
            else:
                # Return zeros for structural component
                batch_size = seq_emb.shape[0]
                device = seq_emb.device
                struct_out = torch.zeros(batch_size, self.prot_struct_encoder.fc.out_features, device=device)
                return torch.cat([seq_out, struct_out], dim=-1)
        return seq_out

    
    def forward(
        self,
        rna_seq_emb,
        prot_seq_emb,
        rna_struct=None,
        prot_struct=None,
    ):
        """
        Forward pass for interaction prediction.
        
        Args:
            rna_seq_emb: (B, rna_seq_dim) RNA sequence embeddings
            prot_seq_emb: (B, prot_seq_dim) Protein sequence embeddings
            rna_struct: (B, C, H, W) RNA structural features (optional)
            prot_struct: (B, C, H, W) Protein structural features (optional)
        
        Returns:
            (B, 1) logits for interaction prediction
        """
        # Encode RNA and protein
        rna_repr = self.encode_rna(rna_seq_emb, rna_struct)
        prot_repr = self.encode_protein(prot_seq_emb, prot_struct)
        
        # Concatenate representations
        combined = torch.cat([rna_repr, prot_repr], dim=-1)
        
        # Fusion and classification
        fused = self.fusion(combined)
        logits = self.classifier(fused)
        
        return logits.squeeze(-1)
    
    def predict_proba(self, *args, **kwargs):
        """Get probability of interaction."""
        logits = self.forward(*args, **kwargs)
        return torch.sigmoid(logits)


class RPIStructuralWithAttention(RPIStructural):
    """
    Extended model with cross-attention between RNA and protein representations.
    """
    
    def __init__(self, *args, num_heads: int = 4, **kwargs):
        super().__init__(*args, **kwargs)
        
        # Determine representation dimension
        if self.use_structural:
            repr_dim = kwargs.get('hidden_dim', 256) + kwargs.get('struct_dim', 128)
        else:
            repr_dim = kwargs.get('hidden_dim', 256)
        
        # Cross-attention layer
        self.cross_attention = nn.MultiheadAttention(
            embed_dim=repr_dim,
            num_heads=num_heads,
            dropout=kwargs.get('dropout', 0.2),
            batch_first=True,
        )
        
        # Update fusion layer for attention output
        fusion_dim = repr_dim * 2
        self.fusion = nn.Sequential(
            nn.Linear(fusion_dim, kwargs.get('hidden_dim', 256)),
            nn.ReLU(),
            nn.Dropout(kwargs.get('dropout', 0.2)),
            nn.Linear(kwargs.get('hidden_dim', 256), kwargs.get('hidden_dim', 256) // 2),
            nn.ReLU(),
            nn.Dropout(kwargs.get('dropout', 0.2)),
        )
        
        self.classifier = nn.Linear(kwargs.get('hidden_dim', 256) // 2, 1)

    
    def forward(
        self,
        rna_seq_emb,
        prot_seq_emb,
        rna_struct=None,
        prot_struct=None,
    ):
        """Forward with cross-attention between RNA and protein."""
        # Encode RNA and protein
        rna_repr = self.encode_rna(rna_seq_emb, rna_struct)
        prot_repr = self.encode_protein(prot_seq_emb, prot_struct)
        
        # Add sequence dimension for attention (B, 1, D)
        rna_repr = rna_repr.unsqueeze(1)
        prot_repr = prot_repr.unsqueeze(1)
        
        # Cross-attention: RNA attends to protein
        rna_attended, _ = self.cross_attention(rna_repr, prot_repr, prot_repr)
        
        # Cross-attention: Protein attends to RNA
        prot_attended, _ = self.cross_attention(prot_repr, rna_repr, rna_repr)
        
        # Remove sequence dimension
        rna_attended = rna_attended.squeeze(1)
        prot_attended = prot_attended.squeeze(1)
        
        # Concatenate and classify
        combined = torch.cat([rna_attended, prot_attended], dim=-1)
        fused = self.fusion(combined)
        logits = self.classifier(fused)
        
        return logits.squeeze(-1)


if __name__ == "__main__":
    # Test the model
    batch_size = 4
    
    # Create dummy inputs
    rna_seq = torch.randn(batch_size, 640)
    prot_seq = torch.randn(batch_size, 480)
    rna_struct = torch.randn(batch_size, 4, 256, 256)
    prot_struct = torch.randn(batch_size, 2, 256, 256)
    
    # Test basic model
    model = RPIStructural(use_structural=True)
    out = model(rna_seq, prot_seq, rna_struct, prot_struct)
    print(f"RPIStructural output shape: {out.shape}")
    
    # Test without structural features
    model_seq_only = RPIStructural(use_structural=False)
    out_seq = model_seq_only(rna_seq, prot_seq)
    print(f"Sequence-only output shape: {out_seq.shape}")
    
    # Test attention model
    model_attn = RPIStructuralWithAttention(use_structural=True, hidden_dim=256, struct_dim=128)
    out_attn = model_attn(rna_seq, prot_seq, rna_struct, prot_struct)
    print(f"Attention model output shape: {out_attn.shape}")
    
    # Count parameters
    total_params = sum(p.numel() for p in model.parameters())
    print(f"\nTotal parameters (with structural): {total_params:,}")
