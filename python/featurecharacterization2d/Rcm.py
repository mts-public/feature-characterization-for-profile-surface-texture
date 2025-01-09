import numpy as np
from scipy.interpolate import interp1d


def inverse_material_ratio(z, p):
    """
    Material ratio curve (Abbott Firestone Curve)
    """

    try:
        z =z.flatten()
        heightintersection = np.sort(z)[::-1]
        material = np.linspace(100.0/len(z), 100, len(z))

        # Level of intersection at given material ratio p (p = 0% => Rcm = 0)
        interpolation_function = interp1d(material, heightintersection, kind='next')
        rcm = interpolation_function(p) - np.max(z)
    except Exception:
        raise RuntimeError("Rcm computation failed.")

    return rcm