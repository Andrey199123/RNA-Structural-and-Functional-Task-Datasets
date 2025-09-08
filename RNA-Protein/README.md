# rpi-baseline (ZHMolGraph-style baseline)

End-to-end baseline for RNA–Protein Interaction (RPI) link prediction inspired by the ZHMolGraph paper:
- Inputs: RNA-FM (RNA) + ProtTrans/ProtT5 (protein) pooled embeddings
- Graph encoder: 2-layer GraphSAGE over an RPI bipartite graph
- Head: VecNN-like 1D-conv + MLP pair scorer
- Splits: random, unknown-RNA, unknown-Protein, unknown-Nodes
- Metrics: Accuracy, MCC, AUROC, AUPRC

> After this baseline is stable, swap ProtTrans → **ESM2** and compare metrics (your first ablation).

---

## 0) Environment (Fedora Linux)

```bash
sudo dnf groupinstall "Development Tools" -y
sudo dnf install -y python3 python3-devel gcc-c++ cmake openssl-devel zlib-devel libffi-devel git
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -U pip wheel

# Torch (choose one)
pip install --index-url https://download.pytorch.org/whl/cpu torch torchvision torchaudio
# or, if you have a compatible NVIDIA setup (CUDA 12.1 wheels):
# pip install --index-url https://download.pytorch.org/whl/cu121 torch torchvision torchaudio

# PyG
pip install torch-geometric

# NLP + utils
pip install transformers tokenizers sentencepiece
pip install pandas numpy scikit-learn tqdm pyyaml biopython
