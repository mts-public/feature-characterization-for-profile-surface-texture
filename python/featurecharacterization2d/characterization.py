import numpy as np

from .watershed import Watershed
from .parameter import feature_parameter

from .Rz import maximum_height
from .Rcm import inverse_material_ratio


def feature_characterization(
    z: np.ndarray,
    dx: float,
    FT: str,
    pruning: str,
    significant: str,
    AT: str,
    stats: str,
):
    """
    Parameters
    ----------
        z : nd.array, float
            vertical profile values in µm
        dx : float
            step size in x-direction in mm
        FT : str
            feature type: {'D', 'V', 'H', 'P'}
        pruning : str
            'None', 'Wolfprune TH/X#', 'Width TH/X#', 'VolS TH', 'DevLength TH'
            Threshold TH in units of corresponding attribute.
            For Wolfprune or Width, the threshold value can also be
            specified as a percentage. In this case, X# of Rz is used
            as the threshold in the case of Wolfprune and X# of Le in
            the case of Width.
        significant : str
            All', 'Closed c', 'Open c', 'Bot N', 'Top N'
            c can be an absolute value in µm or if the value is given
            as a percentage it is interpreted as a material ratio
            from which the height is then determined.
            N specifies the number of top or bot values.
        AT : str
            attribute type {'HDh', 'HDv', 'HDw', 'HDl', 'PVh', 'Curvature', 'Count'}
        stats : str
            'Mean', 'Max', 'Min', 'StdDev', 'Perc X', 'Hist X', 'Sum', 'Density'
            for "Hist", x specifies the number of bins in the histogram.
            for "Perc", x specifies the threshold in the units of the
            corresponding attribute
    Returns
    -------
        xFC : np.ndarray, float
            parameter based on feature characterization
        M : Motif
            structured array of motifs
        meta :
            meta data for further processing (e.g. plotting)
    """

    # parse pruning
    pruning = pruning.replace("%", " %")
    str = pruning.split(" ")
    PT = str[0]
    N = len(str)
    try:
        TH = float(str[np.min([2, N]) - 1])
    except ValueError:
        TH = np.nan

    if N == 2 and np.isnan(TH):
        TH = "opt"

    if N >= 3:
        if PT == "Wolfprune":
            TH = TH / 100.0 * maximum_height(z, dx)
        elif PT == "Width":
            TH = TH / 100.0 * len(z) * dx

    # parse significant
    significant = significant.replace("%", " %")
    str = significant.split(" ")
    Fsig = str[0]
    N = len(str)
    try:
        NIsig = float(str[np.min([2, N]) - 1])
    except ValueError:
        NIsig = np.nan

    if N >= 3:
        NIsig = np.max(z) + inverse_material_ratio(z, NIsig)

    # parse stats
    str = stats.split(" ")
    Astats = str[0]
    try:
        vstats = float(str[-1])
    except ValueError:
        vstats = np.nan

    # feature characterization
    watershed = Watershed(z, dx, FT, PT, TH)
    M = watershed.motifs()
    xFC, M, attr, _, _ = feature_parameter(z, dx, M, Fsig, NIsig, AT, Astats, vstats)
    meta = {
        "attr": attr,
        "nM": len(M),
        "PT": PT,
        "TH": TH,
        "Fsig": Fsig,
        "NIsig": NIsig,
        "AT": AT,
        "Astats": Astats,
        "vstats": vstats,
    }
    return xFC, M, meta
