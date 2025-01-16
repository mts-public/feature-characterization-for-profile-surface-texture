import warnings
import numpy as np

from .motif import Motif
from .attribute import FeatureAttribute


def feature_parameter(
    z: np.ndarray,
    dx: float,
    M: Motif,
    Fsig: str,
    NIsig: float,
    AT: str,
    Astats: str,
    vstats: float,
):
    if z.ndim == 2:
        z = z.reshape(
            -1,
        )
    else:
        z = z
    # step 4: determine significant_features
    I_Nsig = []
    nM = len(M)
    # error handling if M is empty
    if nM == 0:
        warnings.warn("No features detected. Check pruning configuration.")
        xFC = np.nan
        attr = np.nan
        return xFC, M, attr, Fsig, NIsig
    match Fsig:
        case "Open" | "Closed":
            # feature type indicator (FTI=-1 if hills/peaks, FTI=1 if dales/pits)
            FTI = np.sign(z[M[0].ilp.astype(int)] - z[M[0].iv.astype(int)])
            # determine z-values of low-peaks and pits
            zlp = z[M.ilp.astype(int)]
            # determine indices of not significant features
            if Fsig == "Open":
                I_Nsig = np.where(FTI * zlp > FTI * NIsig)[0]
            else:
                zv = z[M.iv.astype(int)]
                I_Nsig = np.where(
                    np.logical_or(FTI * zlp < FTI * NIsig, FTI * zv > FTI * NIsig)
                )[0]
        case "Top" | "Bot":
            # determine attribute values
            attr = FeatureAttribute.compute(z, dx, M, "PVh")
            # determine indices (I) of sorted zv-values in zv
            I_sort = np.argsort(-attr)
            # case NIsig is higher than nM use nM
            NIsig = np.amin([NIsig, nM])
            I_Nsig = I_sort[int(NIsig + 0.5) :]

    # set indicator M.sig zero for not significant motifs
    for i in range(len(I_Nsig)):
        M.sig[I_Nsig[i]] = 0

    # error handling if there are no significant features
    if len(I_Nsig) == nM:
        warnings.warn("All features are declared as not significant.")
        xFC = np.nan
        attr = np.nan
        return xFC, M, attr, Fsig, NIsig

    # step 5: determine attibrute-values of significant features
    attr = FeatureAttribute.compute(z, dx, M, AT)

    ## step 6: attribute statistics
    match Astats:
        case "Mean":
            xFC = np.mean(attr)
        case "Max":
            xFC = np.amax(attr)
        case "Min":
            xFC = np.amin(attr)
        case "StdDev":
            xFC = np.std(attr, ddof=1)
        case "Perc":
            xFC = sum(attr > vstats) / len(attr)
        case "Hist":
            import matplotlib.pyplot as plt

            xFC = plt.hist(attr)
            plt.show()
        case "Sum":
            xFC = sum(attr)
        case "Density":
            xFC = sum(attr) / (dx * len(z) / 10)

    return xFC, M, attr, Fsig, NIsig
