import torch
from torch import nn
from torch_geometric.nn import SAGEConv


class RPISAGE(nn.Module):
    def __init__(self, in_rna, in_prot, hidden_dim=256, out_dim=256, num_layers=2, dropout=0.1):
        super().__init__()
        self.proj_rna = nn.Linear(in_rna, hidden_dim)
        self.proj_prot = nn.Linear(in_prot, hidden_dim)
        self.convs = nn.ModuleList([SAGEConv(hidden_dim, hidden_dim) for _ in range(num_layers)])
        self.act = nn.ReLU()
        self.drop = nn.Dropout(dropout)
        self.out_proj = nn.Linear(hidden_dim, out_dim)

    def forward(self, x_rna, x_prot, edge_index):
        x_rna = self.act(self.proj_rna(x_rna))
        x_prot = self.act(self.proj_prot(x_prot))
        x = torch.cat([x_rna, x_prot], dim=0)
        for conv in self.convs:
            x = conv(x, edge_index)
            x = self.act(x)
            x = self.drop(x)
        z = self.out_proj(x)
        return z
