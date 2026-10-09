from copy import deepcopy

import numpy as np

from .watershed import Watershed
from .parameter import feature_parameter
from .Rz import maximum_height


def default_fc_parameters(z: np.ndarray, dx: float) -> dict:
    """
    Calculates all named parameters of feature characterization defined by
    ISO 21920-2 with the default settings according to ISO 21920-3.

    Parameters
    ----------
        z : nd.array, float
            vertical profile values in µm
        dx : float
            step size in x-direction in mm
    Returns
    -------
        xFC : dict
            named feature parameters {"Rpd", "Rvd", "Rmpc", "Rmvc", "R5p",
            "R5v", "R10z"}
    """
    if z.ndim == 2:
        z = z.reshape(-1)

    # watershed segmentation with default pruning setting
    TH = (5 / 100) * maximum_height(z, dx)
    Mp = Watershed(z, dx, "P", "Wolfprune", TH).motifs()
    Mv = Watershed(z, dx, "V", "Wolfprune", TH).motifs()

    def parameter(M, Fsig, NIsig, AT, Astats):
        # feature_parameter changes the significance of the motifs
        return feature_parameter(z, dx, deepcopy(M), Fsig, NIsig, AT, Astats, np.nan)[0]

    # default feature parameters according to ISO 21920-2
    xFC = {
        "Rpd": parameter(Mp, "All", 1, "Count", "Density"),
        "Rvd": parameter(Mv, "All", 1, "Count", "Density"),
        "Rmpc": parameter(Mp, "All", 1, "Curvature", "Mean"),
        "Rmvc": parameter(Mv, "All", 1, "Curvature", "Mean"),
        "R5p": parameter(Mp, "Top", 5, "PVh", "Mean"),
        "R5v": parameter(Mv, "Bot", 5, "PVh", "Mean"),
    }
    xFC["R10z"] = xFC["R5p"] + xFC["R5v"]
    return xFC
