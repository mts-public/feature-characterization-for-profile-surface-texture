#!/usr/bin/env python3

import unittest

import numpy as np
from featurecharacterization2d import Watershed


class TestWatershedSegementation(unittest.TestCase):
    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)

        # Reference data
        self.z = np.array(
            [
                3.3,
                2,
                1,
                5,
                3.8,
                4,
                1.5,
                1.5,
                3.5,
                2.5,
                2,
                -1,
                0,
                3,
                1.2,
                2,
                -1.2,
                -5,
                -4,
                -4.5,
                -2,
                -2.3,
                1,
                3,
                3,
                3,
                4,
                4.5,
                4.5,
                4,
                1.5,
                1.5,
                3.5,
                4,
                9,
                8,
                -1,
                -1,
                -1,
                -1,
                7,
                7,
                7,
                0,
                0.5,
                1,
                5,
                2,
                6,
                4.5,
                0.5,
                1,
                2,
                -1,
                0,
                3,
                5.2,
                5,
                5.5,
                4,
                7,
            ]
        )
        self.z = self.z - np.mean(self.z)
        self.dx = 0.5

    def test_wolfprune(self):
        watershed = Watershed(self.z, self.dx, "D", PT="Wolfprune", TH=4.1)
        motifs = watershed.motifs()

        # Verfy with matlab solutions
        self.assertTrue(np.all(motifs.iv + 1.0 == np.array([18.0, 38.5, 44.0, 54.0])))
        self.assertTrue(np.all(motifs.ilp + 1.0 == np.array([4.0, 42.0, 49.0, 59.0])))
        self.assertTrue(np.all(motifs.ihp + 1.0 == np.array([35.0, 35.0, 42.0, 49.0])))
        self.assertTrue(
            np.allclose(
                np.array(motifs.ihi) + 1.0,
                np.array([[34.2], [36.1111], [43.1429], [49.3333]]),
            )
        )


if __name__ == "__main__":
    unittest.main()
