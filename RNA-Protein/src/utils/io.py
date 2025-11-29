import random
import re
from pathlib import Path
import numpy as np
import torch
import pandas as pd

def set_seed(seed: int):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)



def read_fasta_to_dict(path) -> dict:
    d = {}
    header = None
    seq_chunks = []
    text = Path(path).read_text()
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                d[header] = "".join(seq_chunks)
            header = line[1:]
            seq_chunks = []
        else:
            seq_chunks.append(line)
    if header is not None:
        d[header] = "".join(seq_chunks)
    return d


def normalize_rna_id(header: str) -> str:
    return header.split()[0]


_UNIPROT_RE = re.compile(r"^[^|]*\|([A-Z0-9]+)\|")


def normalize_prot_id(header: str) -> str:
    m = _UNIPROT_RE.match(header)
    if m:
        return m.group(1)
    return header.split()[0]



def load_npz_embeddings(rna_npz_path, prot_npz_path):
    rna_npz = np.load(rna_npz_path)
    prot_npz = np.load(prot_npz_path)
    rna_ids = list(rna_npz.files)
    prot_ids = list(prot_npz.files)
    Xr = np.stack([rna_npz[k] for k in rna_ids]).astype(np.float32)
    Xp = np.stack([prot_npz[k] for k in prot_ids]).astype(np.float32)
    return (rna_ids, torch.from_numpy(Xr)), (prot_ids, torch.from_numpy(Xp))


def pairs_to_index_tensors(df: pd.DataFrame, rna_ids, prot_ids, device="cpu"):
    rmap = {k: i for i, k in enumerate(rna_ids)}
    pmap = {k: i for i, k in enumerate(prot_ids)}
    offset = len(rna_ids)
    uu, vv = [], []
    missing_rnas, missing_prots = set(), set()
    for rna, prot in zip(df['rna_id'], df['prot_id']):
        if rna not in rmap:
            missing_rnas.add(rna)
            continue
        if prot not in pmap:
            missing_prots.add(prot)
            continue
        uu.append(rmap[rna])
        vv.append(offset + pmap[prot])
    if missing_rnas or missing_prots:
        print(f"Missing IDs when mapping pairs -> indices. "
              f"Missing RNAs: {len(missing_rnas)}, Missing Proteins: {len(missing_prots)}")
    if len(uu) == 0:
        return torch.empty(0, dtype=torch.long, device=device), torch.empty(0, dtype=torch.long, device=device)
    return torch.tensor(uu, dtype=torch.long, device=device), torch.tensor(vv, dtype=torch.long, device=device)
