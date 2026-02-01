#!/usr/bin/env python3
"""
Training script for RNA-protein interaction prediction with structural features.

Uses:
- Sequence embeddings (RNA-FM, ESM2)
- 2D structural features (distance maps, contact maps, base interactions)
"""

import argparse
import json
import random
from pathlib import Path
import numpy as np
import torch
import torch.nn as nn
from torch.optim import Adam
from torch.optim.lr_scheduler import ReduceLROnPlateau
from sklearn.metrics import (
    accuracy_score, matthews_corrcoef, roc_auc_score, average_precision_score
)
from tqdm import tqdm

import sys
sys.path.insert(0, str(Path(__file__).parent))

from models.rpi_structural import RPIStructural, RPIStructuralWithAttention
from data_processing.structural_dataloader import (
    create_dataloaders, load_embeddings_as_dict, 
    StructuralFeatureLoader, RPIDataset, collate_rpi_batch
)


def unknown_entity_split(df, entity_col, train_frac=0.8, val_frac=0.1, seed=42):
    """Split by holding out entities (RNAs or proteins)."""
    random.seed(seed)
    entities = df[entity_col].unique().tolist()
    random.shuffle(entities)
    n = len(entities)
    n_train = int(train_frac * n)
    n_val = int(val_frac * n)
    train_e = set(entities[:n_train])
    val_e = set(entities[n_train:n_train + n_val])
    test_e = set(entities[n_train + n_val:])
    tr = df[df[entity_col].isin(train_e)]
    va = df[df[entity_col].isin(val_e)]
    te = df[df[entity_col].isin(test_e)]
    return tr, va, te


def unknown_nodes_split(df, seed=42):
    """Split where both RNA and protein are unseen in test."""
    random.seed(seed)
    rnas = df['rna_id'].unique().tolist()
    prots = df['prot_id'].unique().tolist()
    random.shuffle(rnas)
    random.shuffle(prots)

    def thirds(xs):
        n = len(xs); a = n // 3; b = 2 * n // 3
        return set(xs[:a]), set(xs[a:b]), set(xs[b:])

    tr_r, va_r, te_r = thirds(rnas)
    tr_p, va_p, te_p = thirds(prots)
    tr = df[df['rna_id'].isin(tr_r) & df['prot_id'].isin(tr_p)]
    va = df[df['rna_id'].isin(va_r) & df['prot_id'].isin(va_p)]
    te = df[df['rna_id'].isin(te_r) & df['prot_id'].isin(te_p)]
    return tr, va, te


def evaluate(model, dataloader, device, threshold=0.5):
    """Evaluate model on a dataloader."""
    model.eval()
    all_labels = []
    all_probs = []
    
    with torch.no_grad():
        for batch in dataloader:
            rna_emb = batch['rna_embeddings'].to(device)
            prot_emb = batch['prot_embeddings'].to(device)
            labels = batch['labels'].to(device)
            
            rna_struct = batch.get('rna_structural')
            prot_struct = batch.get('prot_structural')
            if rna_struct is not None:
                rna_struct = rna_struct.to(device)
            if prot_struct is not None:
                prot_struct = prot_struct.to(device)
            
            logits = model(rna_emb, prot_emb, rna_struct, prot_struct)
            probs = torch.sigmoid(logits)
            
            all_labels.extend(labels.cpu().numpy())
            all_probs.extend(probs.cpu().numpy())
    
    all_labels = np.array(all_labels)
    all_probs = np.array(all_probs)
    all_preds = (all_probs >= threshold).astype(int)

    
    metrics = {
        'accuracy': accuracy_score(all_labels, all_preds),
        'mcc': matthews_corrcoef(all_labels, all_preds),
        'auroc': roc_auc_score(all_labels, all_probs),
        'auprc': average_precision_score(all_labels, all_probs),
    }
    
    return metrics, all_labels, all_probs


def find_best_threshold(labels, probs):
    """Find threshold that maximizes MCC."""
    best_thr, best_mcc = 0.5, -1
    for thr in np.linspace(0.1, 0.9, 17):
        preds = (probs >= thr).astype(int)
        mcc = matthews_corrcoef(labels, preds)
        if mcc > best_mcc:
            best_mcc = mcc
            best_thr = thr
    return best_thr, best_mcc


def train_epoch(model, dataloader, optimizer, criterion, device):
    """Train for one epoch."""
    model.train()
    total_loss = 0
    num_batches = 0
    
    for batch in dataloader:
        rna_emb = batch['rna_embeddings'].to(device)
        prot_emb = batch['prot_embeddings'].to(device)
        labels = batch['labels'].to(device)
        
        rna_struct = batch.get('rna_structural')
        prot_struct = batch.get('prot_structural')
        if rna_struct is not None:
            rna_struct = rna_struct.to(device)
        if prot_struct is not None:
            prot_struct = prot_struct.to(device)
        
        optimizer.zero_grad()
        logits = model(rna_emb, prot_emb, rna_struct, prot_struct)
        loss = criterion(logits, labels)
        loss.backward()
        optimizer.step()
        
        total_loss += loss.item()
        num_batches += 1
    
    return total_loss / num_batches


def main():
    # Get script directory for relative paths
    script_dir = Path(__file__).parent.parent  # RNA-Protein/
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-dir', type=Path, default=None,
                       help='Data directory (default: RNA-Protein/data)')
    parser.add_argument('--split', type=str, default='random', 
                       choices=['random', 'uRNA', 'uPROT', 'uNODES'],
                       help='Split type: random, uRNA, uPROT, or uNODES')
    parser.add_argument('--epochs', type=int, default=50)
    parser.add_argument('--batch-size', type=int, default=32)
    parser.add_argument('--lr', type=float, default=1e-3)
    parser.add_argument('--hidden-dim', type=int, default=256)
    parser.add_argument('--struct-dim', type=int, default=128)
    parser.add_argument('--dropout', type=float, default=0.2)
    parser.add_argument('--use-structural', action='store_true', default=True)
    parser.add_argument('--no-structural', dest='use_structural', action='store_false')
    parser.add_argument('--use-attention', action='store_true')
    parser.add_argument('--struct-size', type=int, default=128, help='Resize structural maps')
    parser.add_argument('--patience', type=int, default=10)
    parser.add_argument('--device', type=str, default='auto')
    parser.add_argument('--save-dir', type=Path, default=None,
                       help='Checkpoint directory (default: RNA-Protein/checkpoints)')
    args = parser.parse_args()
    
    # Set default paths relative to script location
    if args.data_dir is None:
        args.data_dir = script_dir / 'data'
    if args.save_dir is None:
        args.save_dir = script_dir / 'checkpoints'
    
    # Device
    if args.device == 'auto':
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    else:
        device = torch.device(args.device)
    print(f"Using device: {device}")
    
    # Paths
    data_dir = args.data_dir
    pairs_csv = data_dir / 'RPI7317/rpi7317_pairs_balanced.csv'
    rna_npz = data_dir / 'embeddings/rna_fm_rpi7317.npz'
    prot_npz = data_dir / 'embeddings/esm2_35m_rpi7317.npz'
    struct_dir = data_dir / 'structural_features'
    
    # Check if balanced pairs file exists
    if not pairs_csv.exists():
        print(f"ERROR: Balanced pairs file not found: {pairs_csv}")
        print("Please run generate_negative_samples.py first to create the balanced dataset.")
        return
    
    # Load embeddings
    print("Loading embeddings...")
    rna_emb = load_embeddings_as_dict(rna_npz)
    prot_emb = load_embeddings_as_dict(prot_npz)
    print(f"  RNA embeddings: {len(rna_emb)}")
    print(f"  Protein embeddings: {len(prot_emb)}")
    
    # Create structural loader
    struct_loader = None
    if args.use_structural:
        struct_loader = StructuralFeatureLoader(
            rna_2d_dir=struct_dir / 'rna_2d_features',
            protein_2d_dir=struct_dir / 'protein_2d_features',
            rna_interactions_dir=struct_dir / 'rna_interactions',
        )

    
    # Create dataset and split
    print("Creating dataset...")
    import pandas as pd
    
    pairs_csv = data_dir / 'RPI7317/rpi7317_pairs_balanced.csv'
    splits_dir = data_dir / 'interim'
    
    # Check if pre-made splits exist, otherwise create them
    if args.split == 'random':
        # Use simple random split
        full_df = pd.read_csv(pairs_csv, encoding='utf-8-sig')
        full_df = full_df.sample(frac=1, random_state=42).reset_index(drop=True)
        n = len(full_df)
        train_df = full_df[:int(0.8*n)]
        val_df = full_df[int(0.8*n):int(0.9*n)]
        test_df = full_df[int(0.9*n):]
    else:
        # Generate splits using make_splits logic
        full_df = pd.read_csv(pairs_csv, encoding='utf-8-sig')
        full_df = full_df.rename(columns={'RNA names': 'rna_id', 'Protein names': 'prot_id'})
        
        if args.split == 'uRNA':
            train_df, val_df, test_df = unknown_entity_split(full_df, 'rna_id')
        elif args.split == 'uPROT':
            train_df, val_df, test_df = unknown_entity_split(full_df, 'prot_id')
        elif args.split == 'uNODES':
            train_df, val_df, test_df = unknown_nodes_split(full_df)
        
        # Rename back for compatibility
        train_df = train_df.rename(columns={'rna_id': 'RNA names', 'prot_id': 'Protein names'})
        val_df = val_df.rename(columns={'rna_id': 'RNA names', 'prot_id': 'Protein names'})
        test_df = test_df.rename(columns={'rna_id': 'RNA names', 'prot_id': 'Protein names'})
    
    print(f"  Split: {args.split}")
    print(f"  Train: {len(train_df)}, Val: {len(val_df)}, Test: {len(test_df)}")
    
    # Save splits temporarily
    train_df.to_csv(data_dir / 'RPI7317/train_split.csv', index=False)
    val_df.to_csv(data_dir / 'RPI7317/val_split.csv', index=False)
    test_df.to_csv(data_dir / 'RPI7317/test_split.csv', index=False)
    
    # Create datasets
    train_dataset = RPIDataset(
        data_dir / 'RPI7317/train_split.csv', rna_emb, prot_emb,
        struct_loader, args.use_structural, args.struct_size, args.struct_size
    )
    val_dataset = RPIDataset(
        data_dir / 'RPI7317/val_split.csv', rna_emb, prot_emb,
        struct_loader, args.use_structural, args.struct_size, args.struct_size
    )
    test_dataset = RPIDataset(
        data_dir / 'RPI7317/test_split.csv', rna_emb, prot_emb,
        struct_loader, args.use_structural, args.struct_size, args.struct_size
    )
    
    # Create dataloaders
    train_loader = torch.utils.data.DataLoader(
        train_dataset, batch_size=args.batch_size, shuffle=True, collate_fn=collate_rpi_batch
    )
    val_loader = torch.utils.data.DataLoader(
        val_dataset, batch_size=args.batch_size, shuffle=False, collate_fn=collate_rpi_batch
    )
    test_loader = torch.utils.data.DataLoader(
        test_dataset, batch_size=args.batch_size, shuffle=False, collate_fn=collate_rpi_batch
    )

    
    # Get embedding dimensions from data
    sample_batch = next(iter(train_loader))
    rna_seq_dim = sample_batch['rna_embeddings'].shape[1]
    prot_seq_dim = sample_batch['prot_embeddings'].shape[1]
    
    # Create model
    print("Creating model...")
    ModelClass = RPIStructuralWithAttention if args.use_attention else RPIStructural
    model = ModelClass(
        rna_seq_dim=rna_seq_dim,
        prot_seq_dim=prot_seq_dim,
        rna_struct_channels=4,
        prot_struct_channels=2,
        hidden_dim=args.hidden_dim,
        struct_dim=args.struct_dim,
        dropout=args.dropout,
        use_structural=args.use_structural,
    ).to(device)
    
    total_params = sum(p.numel() for p in model.parameters())
    print(f"  Model parameters: {total_params:,}")
    print(f"  Using structural features: {args.use_structural}")
    print(f"  Using attention: {args.use_attention}")
    
    # Training setup
    optimizer = Adam(model.parameters(), lr=args.lr, weight_decay=1e-4)
    scheduler = ReduceLROnPlateau(optimizer, mode='max', factor=0.5, patience=5)
    criterion = nn.BCEWithLogitsLoss()
    
    # Training loop
    args.save_dir.mkdir(parents=True, exist_ok=True)
    print(f"Checkpoints will be saved to: {args.save_dir.absolute()}")
    best_val_mcc = -1
    patience_counter = 0
    
    print("\nTraining...")
    for epoch in range(1, args.epochs + 1):
        # Train
        train_loss = train_epoch(model, train_loader, optimizer, criterion, device)
        
        # Evaluate
        val_metrics, val_labels, val_probs = evaluate(model, val_loader, device)
        best_thr, _ = find_best_threshold(val_labels, val_probs)
        val_metrics_thr, _, _ = evaluate(model, val_loader, device, threshold=best_thr)
        
        scheduler.step(val_metrics_thr['mcc'])
        
        print(f"Epoch {epoch:3d} | Loss: {train_loss:.4f} | "
              f"Val ACC: {val_metrics_thr['accuracy']:.4f} MCC: {val_metrics_thr['mcc']:.4f} "
              f"AUROC: {val_metrics_thr['auroc']:.4f} | thr: {best_thr:.2f}")

        
        # Save best model
        if val_metrics_thr['mcc'] > best_val_mcc:
            best_val_mcc = val_metrics_thr['mcc']
            patience_counter = 0
            
            torch.save({
                'epoch': epoch,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'val_metrics': val_metrics_thr,
                'best_threshold': best_thr,
                'args': vars(args),
            }, args.save_dir / 'best_model.pt')
        else:
            patience_counter += 1
            if patience_counter >= args.patience:
                print(f"Early stopping at epoch {epoch}")
                break
    
    # Final evaluation on test set
    print("\nLoading best model for test evaluation...")
    checkpoint = torch.load(args.save_dir / 'best_model.pt', weights_only=False)
    model.load_state_dict(checkpoint['model_state_dict'])
    best_thr = checkpoint['best_threshold']
    
    test_metrics, _, _ = evaluate(model, test_loader, device, threshold=best_thr)
    
    print("\n" + "="*60)
    print("FINAL TEST RESULTS")
    print("="*60)
    print(f"Accuracy: {test_metrics['accuracy']:.4f}")
    print(f"MCC:      {test_metrics['mcc']:.4f}")
    print(f"AUROC:    {test_metrics['auroc']:.4f}")
    print(f"AUPRC:    {test_metrics['auprc']:.4f}")
    print(f"Threshold: {best_thr:.2f}")
    
    # Save results
    results = {
        'test_metrics': test_metrics,
        'best_val_mcc': best_val_mcc,
        'best_threshold': best_thr,
        'args': vars(args),
    }
    with open(args.save_dir / 'results.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print(f"\nResults saved to {args.save_dir}")


if __name__ == "__main__":
    main()
