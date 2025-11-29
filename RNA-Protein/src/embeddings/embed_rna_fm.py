# This is from Phase 1 (ZHMolGraph Replication, not part of RNA MSA pipeline)
import argparse
from pathlib import Path
import numpy as np
import torch
from tqdm import tqdm

from multimolecule import RnaTokenizer, RnaFmModel
from utils.io import read_fasta_to_dict, normalize_rna_id


def mean_pool(last_hidden_state, attention_mask):
    mask = attention_mask.unsqueeze(-1).float()
    s = (last_hidden_state * mask).sum(1)
    d = mask.sum(1).clamp(min=1e-9)
    return s / d


def encode_sequence_in_chunks(tokenizer, model, seq, device, max_len_tokens=900, add_special_tokens=True):
    model.eval()
    total_weighted = None
    total_tokens = 0.0
    with torch.no_grad():
        for i in range(0, len(seq), max_len_tokens):
            chunk = seq[i:i + max_len_tokens]
            enc = tokenizer(chunk, return_tensors="pt",
                            add_special_tokens=add_special_tokens, padding=False, truncation=False)
            enc = {k: v.to(device) for k, v in enc.items()}
            out = model(**enc).last_hidden_state
            pooled = mean_pool(out, enc["attention_mask"]).squeeze(0)
            n_tokens = float(enc["attention_mask"].sum().item())
            vec = pooled * n_tokens
            total_weighted = vec if total_weighted is None else (total_weighted + vec)
            total_tokens += n_tokens
    return (total_weighted / max(total_tokens, 1.0)).cpu().numpy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="RNA FASTA file")
    ap.add_argument("--out_npz", required=True)
    ap.add_argument("--model_name", default="multimolecule/rnafm")
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--chunk_len", type=int, default=900, help="chunk length (nucleotides) for long RNAs")
    args = ap.parse_args()

    device = torch.device(args.device)
    tok = RnaTokenizer.from_pretrained(args.model_name)
    mdl = RnaFmModel.from_pretrained(args.model_name).to(device)

    raw = read_fasta_to_dict(args.fasta)
    ids = []
    seqs = {}
    for hdr, s in raw.items():
        rid = normalize_rna_id(hdr)
        seqs[rid] = s
        ids.append(rid)

    out = {}
    for rid in tqdm(ids, desc="RNA-FM"):
        pooled = encode_sequence_in_chunks(tok, mdl, seqs[rid], device, max_len_tokens=args.chunk_len)
        out[rid] = pooled

    Path(args.out_npz).parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.out_npz, **out)
    print(f"✅ Saved RNA embeddings to {args.out_npz} (n={len(out)})")


if __name__ == "__main__":
    main()
