import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from featurecharacterization2d import feature_characterization
from featurecharacterization2d import plot_motifs

# Simple self-generated profile to illustrate handling and functionality
z = np.array([
    3.3, 2, 1, 5, 3.8, 4, 1.5, 1.5, 3.5, 2.5, 2, -1, 0, 3, 1.2, 2, -1.2, -5, -4, -4.5, -2,
    -2.3, 1, 3, 3, 3, 4, 4.5, 4.5, 4, 1.5, 1.5, 3.5, 4, 9, 8, -1, -1, -1, -1, 7, 7, 7, 0, 0.5,
    3, 5, 4, 5, 4.5, 0.5, 1, 2, -1, 0, 3, 5.2, 5, 5.5, 4, 7
])
z = (z - np.mean(z))

# Alternatively, load profile from data folder (uncomment if needed)
# z = loadmat('..\data\profiles\Bu_1_56_ak.mat')['z']

# Feature characterization
dx = 0.5e-3
xFC, motifs, meta = feature_characterization(z, dx, "D", "Wolfprune 5 %", "All", "HDh", "Mean")
print(xFC)
plot_motifs(z, dx, motifs, meta["Fsig"], meta["NIsig"])