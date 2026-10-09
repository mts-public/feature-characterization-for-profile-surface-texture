import numpy as np


def interpolation(arr, idx_interpolated):

    idxfloor = np.floor(idx_interpolated).astype(int)
    idxceil = np.ceil(idx_interpolated).astype(int)

    if idxceil == idxfloor:
        return arr[idxceil]

    out = (arr[idxceil] - arr[idxfloor]) / 1.0 * (idx_interpolated- idxfloor) + arr[idxfloor]
    return out

def patch_featureelement(z, x, dx, M, col):
    """
    Creates patch elements (shaded areas) for motifs.
    """
    dir = np.sign(M.ihp - M.ilp)[0].astype(int)
    # all height intersections incl. low-peak (number of intersections varies)
    ihi = np.hstack([M.ilp, np.asarray(M.ihi[0], dtype=float)])
    zlp = interpolation(z, M.ilp)

    patches = []
    for i in range(0, len(ihi) - 1, 2):
        i1 = int(np.abs(np.ceil(dir * ihi[i])))
        i2 = int(np.abs(np.floor(dir * ihi[i + 1])))
        idx = np.arange(i1, i2 + dir, dir)
        xf = np.hstack([ihi[i] * dx, x[idx], ihi[i + 1] * dx])
        zf = np.hstack([zlp, z[idx], zlp])

        patches.append((xf, zf, col))

    return patches

def plot_motifs(z, dx, M, Fsig="All", NIsig=None):
    """
    Parameters
    ----------
        z : nd.array, float
            vertical profile values in µm
        dx : float
            step size in x-direction in mm
        M : Motif
            structure for motifs  
        Fsig : str, optional
               significant features {”All”, ”Open”, ”Closed”, ”Top”, ”Bot”}
        NIsig : float, optional
                nesting index for significant features
    """
    import matplotlib.pyplot as plt

    if len(M) == 0:
        raise ValueError("No features detected. Check pruning configuration.")

    if z.ndim == 2:
        z = z.reshape(-1,)

    x = dx * np.arange(z.size)

    # Plot settings
    linewidth = 0.85
    alpha = 0.4  # Transparency for not significant motifs
    col = [0, 0.4470, 0.7410]  # Color of motifs
    frame_line_style = ['-', ':']
    
    # Plots
    fig, ax = plt.subplots(dpi=300)
    plt.xlabel("profile length in mm")
    plt.ylabel("profile height in μm")
    ax.grid(True)

    # Plot patches (shaded areas)
    for motif in M:
        patches = patch_featureelement(z, x, dx, motif, col)
        for xf, zf, color in patches:
            ax.fill_between(xf, zf, np.ones_like(zf) * interpolation(z, motif.ilp),
                            color=color, alpha=(1 - motif.sig)*alpha + motif.sig, linewidth=0, zorder=2)
        # Plot motif frames
        ilp, iv, ihp = motif.ilp, motif.iv, motif.ihp
        xvals = [(ilp) * dx, (ilp) * dx, (ihp) * dx, (ihp) * dx]
        zvals = [interpolation(z, ilp), interpolation(z, iv), interpolation(z, iv), interpolation(z, ihp)]
        ax.plot(xvals, zvals, color="red", linestyle=frame_line_style[int(motif.sig.item() != 1)],linewidth=0.5)


    # Threshold for Fsig = "Open" or "Closed"
    if Fsig in ["Open", "Closed"] and NIsig is not None and np.isfinite(NIsig):
        ax.axhline(NIsig, linestyle="--", linewidth=1, color="black")

    # Plot settings
    ax.plot(x, z, color="black", linewidth=linewidth, antialiased=True)
    plt.show()