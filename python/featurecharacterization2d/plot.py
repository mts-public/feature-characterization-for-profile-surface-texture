import numpy as np
import matplotlib.pyplot as plt


def interpolation(arr, idx_interpolated):

    idxfloor = np.floor(idx_interpolated).astype(int)
    idxceil = np.ceil(idx_interpolated).astype(int)

    if idxceil == idxfloor:
        return arr[idxceil]

    out = (arr[idxceil] - arr[idxfloor]) / 1.0 * (idx_interpolated- idxfloor) + arr[idxfloor]
    return out


def plot_motif(z, dx, M):
    """
    Parameters
    ----------
        z : nd.array, float
            vertical profile values in µm
        dx : float
            step size in x-direction in mm
        M : Motif
            structure for motifs  
    """
    if z.ndim == 2:
        z = z.reshape(-1,)
    x = dx * np.arange(z.size)
    # Plots
    plt.plot(x, z, color="black")
    plt.xlabel("profile length in mm")
    plt.ylabel("profile height in μm")

    for i in range(len(M)):
        iv = M[i].iv
        ilp = M[i].ilp
        ihp = M[i].ihp
        ihi = np.array(M[i].ihi[0])
        zv = interpolation(z, iv)
        xv = interpolation(x, iv)
        zlp = interpolation(z, ilp)
        xlp = interpolation(x, ilp)
        zhp = interpolation(z, ihp)
        xhp = interpolation(x, ihp)
        # Change tree
        plt.plot([xv, xlp, xlp], [zv, zv, zlp], c="red")
        plt.plot([xv, xhp, xhp], [zv, zv, zhp], c="red")
        # Waterlevel
        direction = np.sign(ihp - ilp)[0].astype(int)
        idxaux = np.hstack([ilp, ihi])
        for i in range(0, idxaux.size - 1, 2):
            i1 = np.abs(np.ceil(direction*idxaux[i]).astype(int))
            i2 = np.abs(np.floor(direction*idxaux[i+1]).astype(int))
            xwater = dx * np.hstack([idxaux[i], np.arange(i1, i2+1*direction, direction), idxaux[i+1]])
            zupper = np.ones(xwater.shape) * zlp
            zlower = np.hstack([zlp, z[i1:i2+1*direction:direction], zlp])
            plt.fill_between(xwater, zlower, zupper, alpha=0.2, linewidth=0, color="tab:blue")