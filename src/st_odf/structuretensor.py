"""
3-D Structure Tensor extracted from fiberorient (MIT licence).

Only depends on numpy and scipy.ndimage — no pkg_resources or dipy required.
"""

from __future__ import annotations

import logging
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.ndimage import gaussian_filter


class StructureTensor:
    """Compute the 3-D structure tensor at every voxel of an image.

    Parameters
    ----------
    d_sigma : float
        Sigma for the Gaussian derivative filters.
    n_sigma : float
        Sigma for the Gaussian neighbourhood filters.
    gaussargs : dict, optional
        Keyword arguments forwarded to ``scipy.ndimage.gaussian_filter``.
        Defaults to ``{'mode': 'nearest', 'cval': 0}``.
    n_jobs : int, optional
        Number of CPU processes for eigenanalysis.  ``None`` → all cores.
    verbose : bool
        Enable progress logging.  Default ``False``.
    """

    def __init__(
        self,
        d_sigma: float,
        n_sigma: float,
        gaussargs: dict | None = None,
        n_jobs: int | None = None,
        verbose: bool = False,
    ) -> None:
        self.d_sigma = d_sigma
        self.n_sigma = n_sigma
        self.gaussargs = gaussargs if gaussargs is not None else {"mode": "nearest", "cval": 0}
        self.n_jobs = n_jobs

        if verbose:
            logging.basicConfig(level=logging.INFO, format="%(message)s")
        self.logger = logging.getLogger()

        self.S: np.ndarray | None = None
        self.evals: np.ndarray | None = None
        self.vectors: np.ndarray | None = None
        self.confidence: np.ndarray | None = None

    # ------------------------------------------------------------------
    def fit(self, img: np.ndarray) -> "StructureTensor":
        """Compute structure tensor and perform eigenanalysis.

        Parameters
        ----------
        img : ndarray
            3-D image array, shape ``(Z, Y, X)``.

        Returns
        -------
        self
        """
        self.S = self._structure_tensor(img)
        logging.info("Calculating eigenvectors/values")

        with ProcessPoolExecutor(self.n_jobs) as pool:
            evals, evectors = zip(*[eig for eig in pool.map(np.linalg.eigh, self.S)])

        self.evals = np.array(evals)
        self.vectors = np.array(evectors)[..., 0]
        self.confidence = self._get_westin()
        logging.info("Done!")
        return self

    # ------------------------------------------------------------------
    def _structure_tensor(self, img: np.ndarray) -> np.ndarray:
        img = np.squeeze(img).astype("float32")

        logging.info("Computing gradient")
        imz, imy, imx = self._compute_gradient(img)

        logging.info("Forming ST elements")
        Szz = gaussian_filter(imz * imz, self.n_sigma, **self.gaussargs)
        Szy = gaussian_filter(imz * imy, self.n_sigma, **self.gaussargs)
        Szx = gaussian_filter(imz * imx, self.n_sigma, **self.gaussargs)
        Syy = gaussian_filter(imy * imy, self.n_sigma, **self.gaussargs)
        Syx = gaussian_filter(imy * imx, self.n_sigma, **self.gaussargs)
        Sxx = gaussian_filter(imx * imx, self.n_sigma, **self.gaussargs)

        S = np.array(
            [[Szz, Szy, Szx],
             [Szy, Syy, Syx],
             [Szx, Syx, Sxx]]
        )
        # np.linalg.eigh requires shape (..., 3, 3)
        return np.moveaxis(S, [0, 1], [3, 4])

    def _compute_gradient(self, img: np.ndarray):
        imz = gaussian_filter(img, self.d_sigma, order=[1, 0, 0], **self.gaussargs)
        imy = gaussian_filter(img, self.d_sigma, order=[0, 1, 0], **self.gaussargs)
        imx = gaussian_filter(img, self.d_sigma, order=[0, 0, 1], **self.gaussargs)
        return imz, imy, imx

    def _get_westin(self) -> np.ndarray:
        t2 = self.evals[..., 2]  # largest
        t1 = self.evals[..., 1]  # middle
        t0 = self.evals[..., 0]  # smallest
        with np.errstate(invalid="ignore"):
            westin = np.where(t2 != 0, (t1 - t0) / t2, np.zeros_like(t2))
        return westin
