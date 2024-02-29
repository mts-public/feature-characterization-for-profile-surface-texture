#!/usr/bin/env python3
import unittest
import numpy as np
from scipy.io import loadmat
from featurecharacterization2d import Watershed
from featurecharacterization2d import feature_characterization


class TestFC(unittest.TestCase):
    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.z = loadmat('data\profiles for case studies\Bu_1_56_ak.mat')['z']
        self.z = self.z - np.mean(self.z)
        self.dx = 0.5

    def test_wolfprune3(self):
        xFC,M = feature_characterization(self.z, self.dx, "D", "Wolfprune 3", "All",  "HDh", "Mean")
        self.assertAlmostEqual(xFC, 4.312365774565255)
        
    def test_Volume1(self):
        xFC,M = feature_characterization(self.z, self.dx, "D", "VolS 1", "All",  "HDh", "Mean")
        self.assertAlmostEqual(xFC, 27.938628855042230)


if __name__ == "__main__":
    unittest.main()