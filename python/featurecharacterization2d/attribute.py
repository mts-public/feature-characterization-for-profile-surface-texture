from .motif import Motif
import numpy as np


class FeatureAttribute(object):
    @classmethod
    def compute(cls, z: np.ndarray, dx: float, M: Motif, AT: np.ndarray):
        """
        Parameters
        ----------
            z : nd.array, float
                vertical profile values in µm
            dx : float
                step size in x-direction in mm
            M : Motif
                Object containing motifs
            AT : str
                attribute type {"Wolfprune", "HDh", "Width", "HDw", "VolS",...
                                "HDv", "DevLength", "HDl", "PVh", "Curvature",...
                                "Count"}
        Returns
        -------
            ATTR : nd.array
                attribute value of given motifs
        """

        I_sig = np.where(M.sig == 1)[0]
        match AT:
            case "Wolfprune" | "HDh":
                ATTR = np.abs(z[M[I_sig].ilp.astype(int)] - z[M[I_sig].iv.astype(int)])
            case "Width" | "HDw":
                ATTR = np.zeros(I_sig.size)
                for i in range(I_sig.size):
                    ATTR[i] = dx * np.max(np.abs(M[I_sig[i]].ihi - M[I_sig[i]].ilp))
            case "VolS" | "HDv":
                ATTR = np.zeros(I_sig.size)
                for i in range(I_sig.size):
                    ATTR[i] = cls.HDvf(z, dx, M[I_sig[i]])
            case "DevLength" | "HDl":
                ATTR = np.zeros(I_sig.size)
                for i in range(I_sig.size):
                    ATTR[i] = cls.HDlf(z, dx, M[I_sig[i]])
            case "PVh":
                FTI = np.sign(z[M[0].ilp.astype(int)] - z[M[0].iv.astype(int)])
                ATTR = -FTI * z[M[I_sig].iv.astype(int)]
            case "Curvature":
                ATTR = np.zeros(I_sig.size)
                for i in range(I_sig.size):
                    ATTR[i] = cls.curvature(z, dx, M[I_sig[i]].iv)
            case "Count":
                ATTR = np.ones(I_sig.size)
        return ATTR

    @staticmethod
    def HDvf(z, dx, Mr):
        dx = dx*1000 # convert from mm to µm
        # all heightintersections incl. low-peak
        ihi = np.hstack([Mr.ilp, np.array(Mr.ihi[0])])
        zlp = z[Mr.ilp.astype(int)]
        A = 0
        i = 0
        direction = np.sign(Mr.ihp - Mr.ilp)
        while i < len(ihi) - 1:
            i1 = np.abs(np.ceil(direction * ihi[i])).astype(int)
            i2 = np.abs(np.floor(direction * ihi[i + 1])).astype(int)
            xf = (
                np.hstack(
                    [
                        ihi[i],
                        np.arange(i1[0], i2[0] + direction[0], direction[0]),
                        ihi[i + 1],
                    ]
                )
                * dx
            )
            zf = np.hstack(
                [
                    zlp,
                    z[np.arange(i1[0], i2[0] + direction[0], direction[0]).astype(int)],
                    zlp,
                ]
            )
            A = A + np.abs(np.trapz(xf, zf - zlp))
            i = i + 2
        HDv = A / (len(z) * dx)
        return HDv

    @staticmethod
    def HDlf(z, dx, Mr):
        z = z/1000 # convert from µm to mm
        zlp = z[Mr.ilp.astype(int)]
        direction = np.sign(Mr.ihp - Mr.ilp)
        ihi_end = Mr.ihi[-1][-1]
        i1 = np.abs(np.ceil(direction * Mr.ilp)).astype(int)
        i2 = np.abs(np.floor(direction * ihi_end)).astype(int)
        zf = z[np.arange(i1[0], i2[0] + direction[0], direction[0]).astype(int)]
        HDl = (
            sum(np.sqrt(1.0 + (np.diff(zf.flatten()) / dx) ** 2)) + np.mod(Mr.ilp, 1)
        ) * dx + np.sqrt((ihi_end - i2) ** 2 * dx**2 + (zlp - z[i2]) ** 2)
        return HDl

    @staticmethod
    def curvature(z, dx, ix):
        dx = dx*1000 # convert from mm to µm
        if np.mod(ix, 1) != 0:
            ix = np.hstack([np.floor(ix), np.ceil(ix)]).astype(int)
        else:
            ix = ix.astype(int)

        cx = np.zeros(len(ix))
        for n in range(len(ix)):
            i = ix[n]
            N = len(z)
            if i == 0:
                dz1 = (
                    - 147.0 * z[0]
                    + 360.0 * z[1]
                    - 450.0 * z[2]
                    + 400.0 * z[3]
                    - 225.0 * z[4]
                    + 72.0 * z[5]
                    - 10.0 * z[6]
                ) / (60.0 * dx)
                dz2 = (
                    812.0 * z[0]
                    - 3132.0 * z[1]
                    + 5265.0 * z[2]
                    - 5080.0 * z[3]
                    + 2970.0 * z[4]
                    - 972.0 * z[5]
                    + 137.0 * z[6]
                ) / (180.0 * (dx) ** 2)
            elif i == 1:
                dz1 = (
                    - 10.0 * z[0]
                    - 77.0 * z[1]
                    + 150.0 * z[2]
                    - 100.0 * z[3]
                    + 50.0 * z[4]
                    - 15.0 * z[5]
                    + 2.0 * z[6]
                ) / (60.0 * dx)
                dz2 = (
                    137.0 * z[0]
                    - 147.0 * z[1]
                    - 255.0 * z[2]
                    + 470.0 * z[3]
                    - 285.0 * z[4]
                    + 93.0 * z[5]
                    - 13.0 * z[6]
                ) / (180.0 * (dx) ** 2)
            elif i == 2:
                dz1 = (
                    2.0 * z[0]
                    - 24.0 * z[1]
                    - 35.0 * z[2]
                    + 80.0 * z[3]
                    - 30.0 * z[4]
                    + 8.0 * z[5]
                    - 1.0 * z[6]
                ) / (60.0 * dx)
                dz2 = (
                    - 13.0 * z[0]
                    + 288.0 * z[1]
                    - 420.0 * z[2]
                    + 200.0 * z[3]
                    + 15.0 * z[4]
                    + 12.0 * z[5]
                    + 2.0 * z[6]
                ) / (180.0 * (dx) ** 2)
            elif i == N - 3:
                dz1 = (
                    z[N - 7]
                    - 8.0 * z[N - 6]
                    + 30.0 * z[N - 5]
                    - 80.0 * z[N - 4]
                    + 35.0 * z[N - 3]
                    + 24.0 * z[N - 2]
                    - 2.0 * z[N - 1]
                ) / (60.0 * dx)
                dz2 = (
                    2.0 * z[N - 7]
                    - 12.0 * z[N - 6]
                    + 15.0 * z[N - 5]
                    + 200.0 * z[N - 4]
                    - 420.0 * z[N - 3]
                    + 228.0 * z[N - 2]
                    - 13.0 * z[N - 1]
                ) / (180.0 * (dx) ** 2)
            elif i == N - 2:
                dz1 = (
                    - 2.0 * z[N - 7]
                    + 15.0 * z[N - 6]
                    - 50.0 * z[N - 5]
                    + 100.0 * z[N - 4]
                    - 150.0 * z[N - 3]
                    + 77.0 * z[N - 2]
                    + 10.0 * z[N - 1]
                ) / (60.0 * dx)
                dz2 = (
                    - 13.0 * z[N - 7]
                    + 93.0 * z[N - 6]
                    - 285.0 * z[N - 5]
                    + 470.0 * z[N - 4]
                    - 255.0 * z[N - 3]
                    - 147.0 * z[N - 2]
                    + 137.0 * z[N - 1]
                ) / (180.0 * (dx) ** 2)
            elif i == N - 1:
                dz1 = (
                    10.0 * z[N - 7]
                    - 72.0 * z[N - 6]
                    + 225.0 * z[N - 5]
                    - 400.0 * z[N - 4]
                    + 450.0 * z[N - 3]
                    - 360.0 * z[N - 2]
                    + 147.0 * z[N - 1]
                ) / (60.0 * dx)
                dz2 = (
                    137.0 * z[N - 7]
                    - 972.0 * z[N - 6]
                    + 2970.0 * z[N - 5]
                    - 5080.0 * z[N - 4]
                    + 5265.0 * z[N - 3]
                    - 3132.0 * z[N - 2]
                    + 812.0 * z[N - 1]
                ) / (180.0 * (dx) ** 2)
            else:
                dz1 = (
                    - z[i - 3]
                    + 9.0 * z[i - 2]
                    - 45.0 * z[i - 1]
                    + 45.0 * z[i + 1]
                    - 9.0 * z[i + 2]
                    + z[i + 3]
                ) / (60.0 * dx)
                dz2 = (
                    2.0 * z[i - 3]
                    - 27.0 * z[i - 2]
                    + 270.0 * z[i - 1]
                    - 490.0 * z[i]
                    + 270.0 * z[i + 1]
                    - 27.0 * z[i + 2]
                    + 2.0 * z[i + 3]
                ) / (180.0 * (dx) ** 2)

            cx[n] = dz2 / ((1.0 + dz1**2) ** (3.0 / 2.0))
        cx = np.mean(cx)
        return cx
