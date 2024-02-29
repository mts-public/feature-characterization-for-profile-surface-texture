#!/usr/bin/env python3

import numpy as np
from scipy.io import loadmat
from featurecharacterization2d import Watershed
from featurecharacterization2d import feature_characterization

z = loadmat('data\profiles for case studies\Bu_1_56_ak.mat')['z']
z = z - np.mean(z)
dx = 0.5

M = feature_characterization(z, dx, "D", "None", "All",  "HDh", "Mean")
print(M)

