"""
Structure Tensor ODF utilities.

Three self-contained functions for building and rendering an Orientation
Distribution Function (ODF) glyph from structure-tensor eigenvectors.
"""

from __future__ import annotations

import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
# ODF from vectors
# ─────────────────────────────────────────────────────────────────────────────

def odf_from_vectors_simple(
    vectors,
    sphere_vertices,
    method: str = "vonmises",
    normalize: bool = True,
    sigma: float | None = None,
    kappa: float | None = None,
    batch_size: int = 100_000,
) -> np.ndarray:
    """
    Compute an ODF from a set of orientation vectors by projecting onto a sphere.

    Each vector votes for the sphere vertices it is most aligned with, using
    one of three kernels.  Because structure-tensor eigenvectors have arbitrary
    sign, **pass antipodally-symmetrised vectors** — i.e.
    ``np.concatenate([vectors, -vectors])`` — to obtain a symmetric ODF.

    Parameters
    ----------
    vectors : array_like, shape (N, 3)
        Input 3-D orientation vectors (need not be unit vectors; normalised
        internally).
    sphere_vertices : array_like, shape (M, 3)
        Points on the unit sphere (e.g. from ``fiberorient.util.make_sphere``).
    method : {'vonmises', 'hist', 'kde'}
        Kernel type:
        * ``'vonmises'`` — von Mises–Fisher exponential kernel, concentration
          controlled by *kappa* (recommended).
        * ``'hist'``     — hard assignment: votes for the nearest hemisphere
          (fast but blocky).
        * ``'kde'``      — Gaussian-like soft kernel controlled by *sigma*.
    normalize : bool
        If ``True``, scale the output so its maximum value is 1.
    sigma : float, optional
        Bandwidth for the ``'kde'`` kernel (default 0.2).
    kappa : float, optional
        Concentration for the ``'vonmises'`` kernel (default 20).
        Higher values → sharper, more peaked ODF.
    batch_size : int
        Number of vectors processed per batch to limit memory usage.

    Returns
    -------
    odf_on_sphere : ndarray, shape (M,)
        ODF amplitude at each sphere vertex.
    """
    vectors = np.asarray(vectors, dtype=np.float64)
    vectors = vectors / (np.linalg.norm(vectors, axis=1, keepdims=True) + 1e-10)

    sphere_vertices = np.asarray(sphere_vertices, dtype=np.float64)
    sphere_vertices = sphere_vertices / (
        np.linalg.norm(sphere_vertices, axis=1, keepdims=True) + 1e-10
    )

    n_sphere = sphere_vertices.shape[0]
    n_vectors = vectors.shape[0]
    odf_on_sphere = np.zeros(n_sphere, dtype=np.float64)

    for start in range(0, n_vectors, batch_size):
        batch = vectors[start : start + batch_size]
        dots = np.dot(sphere_vertices, batch.T)   # (M, batch)

        if method == "hist":
            odf_on_sphere += np.sum(dots > 0.95, axis=1)
        elif method == "kde":
            s = sigma if sigma is not None else 0.2
            odf_on_sphere += np.sum(np.exp((dots - 1) / s**2), axis=1)
        elif method == "vonmises":
            k = kappa if kappa is not None else 20
            odf_on_sphere += np.sum(np.exp(k * (dots - 1)), axis=1)
        else:
            raise ValueError(f"Unknown method '{method}'. Choose 'vonmises', 'hist', or 'kde'.")

    if normalize and odf_on_sphere.max() > 0:
        odf_on_sphere /= odf_on_sphere.max()

    return odf_on_sphere


# ─────────────────────────────────────────────────────────────────────────────
# Generalized Fractional Anisotropy
# ─────────────────────────────────────────────────────────────────────────────

def generalized_fractional_anisotropy(odf_on_sphere: np.ndarray) -> float:
    """
    Compute Generalized Fractional Anisotropy (GFA) from an ODF.

    GFA measures how much the ODF deviates from isotropy:
    * 0 → perfectly isotropic (uniform sphere)
    * 1 → perfectly anisotropic (single spike)

    Parameters
    ----------
    odf_on_sphere : array_like, shape (M,)
        ODF amplitude values at each sphere vertex.

    Returns
    -------
    gfa : float
        GFA in [0, 1].
    """
    odf = np.asarray(odf_on_sphere, dtype=np.float64)
    n = len(odf)
    mean = odf.mean()
    numerator = np.sum((odf - mean) ** 2)
    denominator = np.sum(odf ** 2)
    if denominator == 0:
        return 0.0
    gfa = np.sqrt((n / (n - 1)) * (numerator / denominator))
    return float(np.clip(gfa, 0.0, 1.0))


# ─────────────────────────────────────────────────────────────────────────────
# Anatomical rotation
# ─────────────────────────────────────────────────────────────────────────────

def apply_anatomical_rotation(
    sphere_vertices: np.ndarray,
    R: np.ndarray,
) -> np.ndarray:
    """
    Rotate sphere vertices from image (z, y, x) space into anatomical space.

    Structure-tensor eigenvectors are in image axis order ``(v_z, v_y, v_x)``.
    *R* is a 3×3 matrix whose rows define how each output component is formed
    from those image-space components::

        v_out = R @ v_image

    The rotation matrix is dataset-specific and stored in
    ``datasets_config.json`` under ``"anatomical_rotation"``.

    Parameters
    ----------
    sphere_vertices : ndarray, shape (N, 3)
        Sphere vertex coordinates in image (z, y, x) space.
    R : ndarray, shape (3, 3)
        Rotation / transform matrix.

    Returns
    -------
    ndarray, shape (N, 3)
        Rotated vertex coordinates.
    """
    return (np.asarray(R) @ np.asarray(sphere_vertices).T).T
