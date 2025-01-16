from copy import deepcopy
import numpy as np

from .motif import Motif
from .attribute import FeatureAttribute


class Watershed(object):
    def __init__(self, z, dx, FT, PT=None, TH=0.0) -> None:
        """
        Parameters
        ----------
            z : nd.array, float
                vertical profile values in µm
            dx : float
                step size in x-direction in mm
            FT : str
                feature type: {'D', 'V', 'H', 'P'}
            PT : str, optional(default=None)
                pruning type: {'None', 'Wolfprune', 'Width', 'VolS', 'DevLength'}
            TH : float or str
                threshold for pruning (not needed if PT = 'None').
        """

        if z.ndim ==2:
            self.z = z.reshape(-1,)
        else:
            self.z = z
        self.dx = dx
        self.FT = FT
        self.PT = PT
        self.TH = TH

    def motifs(self) -> Motif:
        """
        Returns
        -------
            M : Motif
                structure array with motifs with four members
                (referring to Dale-motif):
                M.iv    - (interpolated) index of pit
                M.ilp   - (interpolated) index of low-peak
                M.ihp   - (interpolated) index of high-peak
                M.ihi   - (interpolated) index of heightintersection
                M.sig   - indicator for significant features
        """
        # step 1: determine indices of all peaks and pits

        # invert z-values if hill-motifs are searched. Allows the rest of the code
        # refers to dale-motifs
        if self.FT == "H" or self.FT == "P":
            z = -1.0 * self.z
        else:
            z = self.z
        # consider only the first of any adjacent pairs of equal values so that
        # plateaus are single points
        mask = z[:-1] != z[1:]
        iNeq = np.hstack([np.array([0]), np.where(mask)[0] + 1])
        # determine the slope for each point. s=-1 neg. slope; s=1 pos. slope
        s = np.sign(np.diff(z[iNeq]))
        # change from 1 to -1 corresponds to a peak. -1 to 1 to pit
        ipNeq = np.where(np.diff(s) == -2)[0] + 1
        ivNeq = np.where(np.diff(s) == 2)[0] + 1
        # peaks and pits indices into the original index vector
        ipv = np.hstack([iNeq[ipNeq], iNeq[ivNeq]], dtype=np.float64)
        # examine peaks and pits that are plateaus
        vp_idx = np.where(z[ipv.astype(int)] == z[ipv.astype(int) + 1])[0]
        for j in vp_idx:
            # number of equal high values on the plateau
            n_plateau = np.where(z[ipv[j].astype(int) :] != z[ipv[j].astype(int)])[0][0]
            # replace index with (interpolated) index of middle of plateau
            ipv[j] = ipv[j] + (n_plateau - 1.0) / 2.0
        # (interpolated) indices of peaks and pits (plateaus taken into account)
        ip = ipv[: ipNeq.size]
        iv = ipv[ipNeq.size :]

        # step 2: determine motifs

        # just keep indices of pits that are enclosed with peaks
        iv = iv[(iv > ip[0]) & (iv < ip[-1])]
        # enrich structure array M with information for each motif such as pit
        # (iv), low-peak (lp) and high-peak (hp) and the height intersection (ihi)
        nM = iv.size
        M = Motif()
        for k in range(nM):
            ilp, ihp = self.get_ilp_ihp(z, [ip[k], ip[k + 1]])
            ihi = self.height_intersections(z, ilp, ihp)
            M.append(values={"iv": iv[k], "ilp": ilp, "ihp": ihp, "ihi": ihi, "sig": 1})

        # step 3: pruning (see pruning cases in readme.md)

        # skip pruning if pruning type (PT) is "None"
        if (self.PT is not None) and (self.PT != "None"):
            # determine attribute-values of each motif
            ATTR = FeatureAttribute.compute(z, self.dx, M, self.PT)
            # find optimal limit for maximum periodicity (if requested)
            if self.TH == "opt":
                self.TH = self.optimal_periodicity(
                    z, self.dx, deepcopy(M), nM, ATTR, self.PT
                )
            # prune aslong minimal attribute value is lower than given threshold
            while min(ATTR) < self.TH:
                M, nM, ATTR = self.prune_min_motif(z, self.dx, M, nM, ATTR, self.PT)
        return M

    def get_ilp_ihp(self, z, ip_surr) -> tuple[float, float]:
        """
        determine indices of low-peak and high-peak based of 2 given indices

        Parameters
        ----------
            z : nd.array, float
                vertical profile values

            ip_surr : np.array or list
                (interpolated) indices of the two surrounding peaks of examined pit
        Returns
        -------
            ilp : float
                (interpolated) index of low-peak
            ihp : float
                (interpolated) index of high-peak
        """
        I = np.argmin([z[int(ip_surr[0])], z[int(ip_surr[1])]])
        ilp = ip_surr[I]
        ihp = ip_surr[1 - I]
        return ilp, ihp

    def height_intersections(self, z, ilp, ihp) -> list:
        """
        height intersection function

        Parameters
        ----------
            z : nd.array, float
                vertical profile values
            ilp : float
                (interpolated) index of low-peak
            ihp : float
                (interpolated) index of high-peak
        Returns
        -------
            ihi : nd.array
                (interpolated) index of height intersection
        """
        # direction in which to search ihi outgoing form low-peak
        direction = np.sign(ihp - ilp).astype(int)
        zlp = z[int(np.round(ilp))]
        ihi = []
        j = int(np.round(ilp)) + direction*(np.where(z[int(np.round(ilp)):int(np.round(ihp)):direction] != zlp)[0][0])
        # getting height-intesections (based on crossing-the-line-segmentation)
        while j != (ihp).astype(int):
            if (z[j] < zlp and z[j + direction] >= zlp) or (z[j] >= zlp and z[j + direction] < zlp):
                ihi.append(j + direction * (zlp - z[j]) / (z[j + direction] - z[j]))
            j += direction
        return ihi


    def optimal_periodicity(self, z, dx, M, nM, ATTR, PT) -> float:
        """
        optimal threshold function
        """
        # minimal Q-Value
        Qmin = 3
        # set default threshold for the case that Qmin is never exceeded
        TH = 100
        # prune until just two motifs are left
        while nM > 2:
            # parameter Q as a measure for periodicity
            Q = np.mean(ATTR) / np.std(ATTR, ddof=1)
            # if Q is greater than Qmin then overwrite Qmin with the current
            # Q-value and TH with minimal attribute value
            if Q > Qmin:
                Qmin = Q
                TH = min(ATTR)
            M, nM, ATTR = self.prune_min_motif(z, dx, M, nM, ATTR, PT)
        return TH

    def prune_min_motif(self, z, dx, M, nM, ATTR, PT) -> tuple[Motif, int, np.ndarray]:
        """
        prune min-motif function
        """
        # row-index of minimal attribute value
        rmin = np.argmin(ATTR)
        # save motif with minimal attribute-value temporarily, delete
        # corresponding entry in motif-array and attribute-vector. update nM
        Mmin = M[rmin]
        del M[rmin]
        ATTR = np.delete(ATTR, rmin)
        nM = nM - 1
        # determine row-index of motif which is to update (rU). dir=-1: left of
        # min-motif. dir=1 right of min-motif. rU = rmin if dir=1 because
        # min-motif was deleted in motif array. else rU = rmin - 1.
        direction = np.sign(Mmin.ilp - Mmin.iv).astype(int)[0]
        rU = rmin - int(direction == -1)
        # case 1: if in that direction is border no further steps are required
        if (rU == -1) or (rU > nM - 1):
            return M, nM, ATTR
        # case 2: if low-peak of min-motif and motif to update is the same then
        # determine low-peak and high-peak of motif to update. Update height
        # intersection and attribute-value (see behind if-query)
        if M[rU].ilp == Mmin.ilp:
            M.ilp[rU], M.ihp[rU] = self.get_ilp_ihp(z, [M[rU].ihp[0], Mmin.ihp[0]])
        # case 3: low-peak of min-motif is equal the high-peak of motif to
        # update. replace high-peak of motif to merge with high-peak of
        # min-motif
        else:
            M.ihp[rU] = Mmin.ihp[0]
            # case 3.1: if low-peak of motif to merge is lower or equal than
            # pit of min-motif no further steps are required
            if z[M[rU].ilp.astype(int)] <= z[Mmin.iv.astype(int)]:
                return M, nM, ATTR
        # case 3.2: low-peak of motif to merge is higher than pit of
        # min-motif: refresh height intersection and attribute-value
        # (see behind if-query)

        # update height intersection and attribute-value of motif to update for
        # case 2 and 3.2
        M.ihi[rU] = self.height_intersections(z, M[rU].ilp[0], M[rU].ihp[0])
        ATTR[rU] = FeatureAttribute.compute(z, dx, M[rU], PT)[0]

        return M, nM, ATTR