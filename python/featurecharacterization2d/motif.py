from typing import Union
import numpy as np


class Motif:
    """
    Struct for a motif
    """

    def __init__(self) -> None:
        # Placeholder
        self.iv = None
        self.ilp = None
        self.ihp = None
        # ihi can be a list of different sizes.
        # Here it is a list of lists.
        self.ihi = []
        self.sig = None

        self.size = 0

    def append(self, values: Union[dir]) -> None:
        """
        Adds only a single motif.
        """
        if self.iv is None:
            self.iv = np.array([values["iv"]])
        else:
            self.iv = np.hstack([self.iv, values["iv"]])

        if self.ilp is None:
            self.ilp = np.array([values["ilp"]])
        else:
            self.ilp = np.hstack([self.ilp, values["ilp"]])

        if self.ihp is None:
            self.ihp = np.array([values["ihp"]])
        else:
            self.ihp = np.hstack([self.ihp, values["ihp"]])

        self.ihi.append(values["ihi"])

        if self.sig is None:
            self.sig = np.array([values["sig"]])
        else:
            self.sig = np.hstack([self.sig, values["sig"]])

        self.size += 1
        return None

    def __len__(self) -> int:
        return self.size

    def __delitem__(self, index: Union[list, int, np.ndarray, np.integer]) -> None:
        """
        Removes motifs by index.
        """
        if isinstance(index, (list, np.ndarray)):
            self.size -= len(index)
            for id in index:
                del self.ihi[id]
        else:
            self.size -= 1
            del self.ihi[index]
        self.iv = np.delete(self.iv, index)
        self.ilp = np.delete(self.ilp, index)
        self.ihp = np.delete(self.ihp, index)
        self.sig = np.delete(self.sig, index)
        return None

    def __getitem__(self, index: Union[list, int, np.int_, np.ndarray]):
        """
        Returns a single motif or multiple motifs by index.
        """
        if isinstance(index, (int, np.integer)):
            motif_object = Motif()
            values = {
                "iv": self.iv[index],
                "ilp": self.ilp[index],
                "ihp": self.ihp[index],
                "ihi": self.ihi[index],
                "sig": self.sig[index],
            }
            motif_object.append(values)
            return motif_object
        else:
            motif_object = Motif()
            for id in index:
                values = {
                    "iv": self.iv[id],
                    "ilp": self.ilp[id],
                    "ihp": self.ihp[id],
                    "ihi": self.ihi[id],
                    "sig": self.sig[id],
                }
                motif_object.append(values)
            return motif_object