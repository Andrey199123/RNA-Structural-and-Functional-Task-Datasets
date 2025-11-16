#!/usr/bin/env python3

import argparse, re
from pathlib import Path
import numpy as np
import torch
from tqdm import tqdm
from transformers import AutoTokenizer, AutoModel, T5Tokenizer, T5EncoderModel

def read_fasta_to_dict(path) -> dict:
    d, header, seq = {}, None, []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line: continue
        if line.startswith(">"):
            if header is not None: d[header] = "".join(seq)
            header, seq = line[1:], []
        else:
            seq.append(line)
    if header is not None: d[header] = "".join(seq)
    return d

_UNIPROT_RE = re.compile(r"^[^|]*\|([A-Z0-9]+)\|")  # sp|Q01844|...

def normalize_prot_id(header: str) -> str:
    m = _UNIPROT_RE.match(header)
    return m.group(1) if m else header.split()[0]

# ---------- tokenization/pooling ----------
def space_separate(seq: str) -> str:
    return " ".join(list(seq.strip().upper()))

def mean_pool(last_hidden_state, attention_mask):
    mask = attention_mask.unsqueeze(-1).float()
    s = (last_hidden_state * mask).sum(1)
    d = mask.sum(1).clamp(min=1e-9)
    return s / d  # [B, H]

def encode_in_chunks(tokenizer, model, seq, device, use_space, max_len=700, add_special_tokens=True):
    model.eval()
    total, n_tokens_total = None, 0.0
    with torch.inference_mode():
        for i in range(0, len(seq), max_len):
            chunk = seq[i:i+max_len]
            text  = space_separate(chunk) if use_space else chunk
            enc = tokenizer(text, return_tensors="pt", add_special_tokens=add_special_tokens,
                            padding=False, truncation=False)
            enc = {k: v.to(device) for k, v in enc.items()}
            out = model(**enc).last_hidden_state  # [1, T, H]
            pooled = mean_pool(out, enc["attention_mask"]).squeeze(0)  # [H]
            n_tokens = float(enc["attention_mask"].sum().item())
            vec = pooled * n_tokens
            total = vec if total is None else total + vec
            n_tokens_total += n_tokens
            del enc, out, pooled
    return (total / max(n_tokens_total, 1.0)).cpu().numpy()

def load_model(model_name, backend, device, dtype):
    if backend == "t5":
        tok = T5Tokenizer.from_pretrained(model_name, do_lower_case=False, use_fast=False, legacy=True)
        mdl = T5EncoderModel.from_pretrained(model_name, torch_dtype=dtype, low_cpu_mem_usage=True)
    else:
        tok = AutoTokenizer.from_pretrained(model_name)
        mdl = AutoModel.from_pretrained(model_name, torch_dtype=dtype, low_cpu_mem_usage=True)
    mdl.to(device).eval()
    return tok, mdl

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out_npz", required=True)
    ap.add_argument("--model_name", default="Rostlab/prot_bert")
    ap.add_argument("--backend", choices=["auto","t5"], default="auto",
                    help="use 'auto' for ProtBert/ESM2; 't5' for ProtT5")
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--dtype", choices=["fp32","fp16","bf16"], default="fp32")
    ap.add_argument("--chunk_len", type=int, default=700)
    ap.add_argument("--space_separated", action="store_true",
                    help="set for models expecting spaced tokens (ProtBert/ProtT5)")
    args = ap.parse_args()

    dtype = {"fp32": torch.float32, "fp16": torch.float16, "bf16": torch.bfloat16}[args.dtype]
    device = torch.device(args.device)

    tok, mdl = load_model(args.model_name, args.backend, device, dtype)

    raw = read_fasta_to_dict(args.fasta)
    ids, seqs = [], {}
    for hdr, s in raw.items():
        pid = normalize_prot_id(hdr)
        ids.append(pid); seqs[pid] = s

    out = {}
    for pid in tqdm(ids, desc=f"Embed {args.model_name}"):
        out[pid] = encode_in_chunks(tok, mdl, seqs[pid], device,
                                    use_space=args.space_separated, max_len=args.chunk_len)

    Path(args.out_npz).parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.out_npz, **out)
    print(f"✅ Saved protein embeddings to {args.out_npz} (n={len(out)})")

if __name__ == "__main__":
    main()
