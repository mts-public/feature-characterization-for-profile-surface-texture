#!/usr/bin/env python3

import unittest

import numpy as np
from scipy.io import loadmat
from featurecharacterization2d import feature_characterization


class TestFC(unittest.TestCase):
    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.z = loadmat("data/profiles/Bu_1_56_ak.mat")["z"]
        self.z = self.z - np.mean(self.z)
        self.dx = 0.5

    def test_wolfprune5_open_max(self):
        xFC, _, _ = feature_characterization(
            self.z, self.dx, "D", "Wolfprune 5 %", "Open 10 %", "HDv", "Max"
        )

        # Verify with matlab solutions
        self.assertAlmostEqual(xFC, 0.0216294517099762)

    def test_wolfprune5_top_max(self):
        xFC, _, _ = feature_characterization(
            self.z, self.dx, "D", "Wolfprune 5 %", "Top 10 %", "HDv", "Max"
        )

        # Verfy with matlab solutions
        self.assertAlmostEqual(xFC, 0.0134717823931572)

    def test_wolfprune3(self):
        xFC, _, _ = feature_characterization(
            self.z, self.dx, "D", "Wolfprune 3", "All", "HDh", "Mean"
        )

        # Verify with matlab solutions
        self.assertAlmostEqual(xFC, 4.312365774565255)

    def test_Volume1(self):
        xFC, _, _ = feature_characterization(
            self.z, self.dx, "D", "VolS 1", "All", "HDh", "Mean"
        )

        # Verify with matlab solutions
        self.assertAlmostEqual(xFC, 27.938628855042230)

    def test_Width50(self):
        xFC, _, _ = feature_characterization(
            self.z, self.dx, "D", "Width 50", "All", "HDh", "Mean"
        )

        # Verify with matlab solutions
        self.assertAlmostEqual(xFC, 26.461823860088230)

    def test_DevLength50(self):
        xFC, _, _ = feature_characterization(
            self.z, self.dx, "D", "DevLength 50", "All", "HDh", "Mean"
        )

        # Verify with matlab solutions
        self.assertAlmostEqual(xFC, 27.519549401462786)


if __name__ == "__main__":
    unittest.main()
