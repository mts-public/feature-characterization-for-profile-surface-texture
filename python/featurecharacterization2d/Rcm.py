import numpy as np


def inverse_material_ratio(z, p):
    """
    Inverse material ratio Rcm(p) in µm according to ISO 21920-2, 4.5.1.4:
    level of intersection at the material ratio p in % relative to the
    maximum height (Rcm(0) = 0). The material ratio curve is given by the
    pairs (k/n, c_k) of the heights c_k sorted in descending order
    (Annex C), between them the next value c_k is used.

    Parameters
    ----------
        z : nd.array, float
            vertical profile values in µm
        p : float or nd.array, float
            material ratio in %
    Returns
    -------
        rcm : float or nd.array, float
            inverse material ratio in µm
    """
    p = np.asarray(p, dtype=float)
    if np.any((p < 0.0) | (p > 100.0)):
        raise ValueError("The material ratio p has to be within 0 % and 100 %.")

    z = np.asarray(z, dtype=float).reshape(-1)
    n = z.size
    # material ratio curve (Abbott Firestone Curve)
    heightintersection = np.sort(z)[::-1]
    # index of the smallest material ratio k/n >= p (1-based)
    # (tolerance for rounding errors of p*n/100 at the sampling points)
    k = np.maximum(1, np.ceil(p * n / 100.0 - 1e-9)).astype(int)
    rcm = heightintersection[k - 1] - heightintersection[0]
    return rcm[()] if rcm.ndim == 0 else rcm
