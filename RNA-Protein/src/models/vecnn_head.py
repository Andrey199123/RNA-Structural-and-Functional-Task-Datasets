import torch
from torch import nn


class VecNNHead(nn.Module):
    """
    Simple pair-scoring head (VecNN flavor):
      input: z_u, z_v  (node embeddings for RNA and Protein)
      features: [z_u, z_v, |z_u - z_v|, z_u * z_v]
      conv1d (kernel=3, groups=1) over the concatenated feature vector (as length L)
      -> MLP -> logit
    """
    def __init__(self, z_dim, hidden=256):
        super().__init__()
        self.feat_dim = 4 * z_dim
        self.conv = nn.Conv1d(1, 1, kernel_size=3, padding=1)
        self.mlp = nn.Sequential(
            nn.Linear(self.feat_dim, hidden),
            nn.ReLU(),
            nn.Linear(hidden, 1)
        )

    def forward(self, z_u, z_v):
        # z_u, z_v: [B, D]
        pair = torch.cat([z_u, z_v, torch.abs(z_u - z_v), z_u * z_v], dim=-1)  # [B, 4D]
        x = pair.unsqueeze(1)               # [B, 1, 4D]
        x = self.conv(x).squeeze(1)         # [B, 4D]
        logit = self.mlp(x).squeeze(-1)     # [B]
        return logit
