"""
st_odf — Structure Tensor ODF utilities.

Provides three functions used by the cutout_st_odf_glyph notebook:
    odf_from_vectors_simple       — build a directional ODF on a sphere
    generalized_fractional_anisotropy — compute GFA from an ODF
    apply_anatomical_rotation     — rotate sphere vertices into anatomical space
"""

from .utils import (
    odf_from_vectors_simple,
    generalized_fractional_anisotropy,
    apply_anatomical_rotation,
)

__all__ = [
    "odf_from_vectors_simple",
    "generalized_fractional_anisotropy",
    "apply_anatomical_rotation",
]
