import torch
from torch import nn


class VecNNHead(nn.Module):
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
        pair = torch.cat([z_u, z_v, torch.abs(z_u - z_v), z_u * z_v], dim=-1)
        x = pair.unsqueeze(1)
        x = self.conv(x).squeeze(1)
        logit = self.mlp(x).squeeze(-1)
        return logit
