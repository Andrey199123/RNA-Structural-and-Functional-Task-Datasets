"""
RNA Multiple Sequence Alignment (MSA) Pipeline
Generates structure-aware MSAs and 2D feature maps for RNA sequences
"""

from .pipeline import RNAMSAPipeline
from .feature_extraction import extract_2d_features

__all__ = ['RNAMSAPipeline', 'extract_2d_features']
