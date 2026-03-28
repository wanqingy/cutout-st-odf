"""
st_odf — Structure Tensor ODF utilities.

Provides the full pipeline without depending on fiberorient or dipy:
    StructureTensor               — compute voxel-wise orientation vectors
    make_sphere                   — Fibonacci-sampled unit sphere mesh
    odf_from_vectors_simple       — build a directional ODF on a sphere
    generalized_fractional_anisotropy — compute GFA from an ODF
    apply_anatomical_rotation     — rotate sphere vertices into anatomical space
"""

from .structuretensor import StructureTensor
from .utils import (
    make_sphere,
    odf_from_vectors_simple,
    generalized_fractional_anisotropy,
    apply_anatomical_rotation,
)

__all__ = [
    "StructureTensor",
    "make_sphere",
    "odf_from_vectors_simple",
    "generalized_fractional_anisotropy",
    "apply_anatomical_rotation",
]
