#!/usr/bin/env python3

import unittest

import numpy as np
from scipy.io import loadmat
from featurecharacterization2d import inverse_material_ratio
from featurecharacterization2d import maximum_height


class TestFC(unittest.TestCase):
    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.z = loadmat("data/profiles/Bu_1_56_ak.mat")["z"]
        self.z = self.z - np.mean(self.z)
        self.dx = 0.5

        # random generated c values for rcm by matlab
        self.c = 100 * np.array(
            [
                0.162182308193243,
                0.794284540683907,
                0.311215042044805,
                0.528533135506213,
                0.165648729499781,
            ]
        )

    def test_rcm(self):
        rcm = inverse_material_ratio(self.z, self.c)

        # verify with matlab results
        res = np.allclose(
            rcm,
            np.array(
                [
                    -12.0236795096572,
                    -30.7178301888893,
                    -16.8151745219686,
                    -29.1302084334353,
                    -12.0687978035015,
                ]
            ),
        )

        self.assertTrue(res)

    def test_rz(self):
        rz = maximum_height(self.z, self.dx)

        # verify with matlab result
        self.assertAlmostEqual(rz, 29.4112166847386)


if __name__ == "__main__":
    unittest.main()
