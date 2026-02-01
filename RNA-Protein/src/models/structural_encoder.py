#!/usr/bin/env python3
"""
CNN encoder for 2D structural features (distance maps, contact maps).

Processes image-like structural features and outputs a fixed-size embedding.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F


class StructuralCNN(nn.Module):
    """
    CNN encoder for 2D structural feature maps.
    
    Takes (B, C, H, W) input and outputs (B, out_dim) embedding.
    """
    
    def __init__(
        self,
        in_channels: int,
        out_dim: int = 256,
        hidden_channels: list = [32, 64, 128],
        dropout: float = 0.1,
    ):
        super().__init__()
        
        layers = []
        prev_ch = in_channels
        
        for ch in hidden_channels:
            layers.extend([
                nn.Conv2d(prev_ch, ch, kernel_size=3, padding=1),
                nn.BatchNorm2d(ch),
                nn.ReLU(inplace=True),
                nn.MaxPool2d(2),  # Halves spatial dimensions
                nn.Dropout2d(dropout),
            ])
            prev_ch = ch
        
        self.conv_layers = nn.Sequential(*layers)
        
        # Global average pooling + projection
        self.global_pool = nn.AdaptiveAvgPool2d(1)
        self.fc = nn.Linear(hidden_channels[-1], out_dim)
    
    def forward(self, x):
        """
        Args:
            x: (B, C, H, W) structural feature maps
        Returns:
            (B, out_dim) embedding
        """
        x = self.conv_layers(x)
        x = self.global_pool(x)  # (B, C, 1, 1)
        x = x.view(x.size(0), -1)  # (B, C)
        x = self.fc(x)
        return x
