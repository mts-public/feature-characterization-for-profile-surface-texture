import numpy as np
from scipy.interpolate import interp1d

from scipy.io import loadmat
from featurecharacterization2d import feature_characterization
  



z = loadmat("data/profiles/Bu_1_56_ak.mat")["z"]
dx=0.5
p=10

heightintersection = np.sort(z[:, 0])[::-1]
material = np.linspace(100.0/len(z), 100, len(z))

# Level of intersection at given material ratio p (p = 0% => Rcm = 0)
interpolation_function = interp1d(material, heightintersection, kind='next')
rcm = interpolation_function(p) - np.max(z)

print(rcm)

R,_,_ = feature_characterization(z, dx, 'D', "Wolfprune 5 %", "Open 10 %", "HDl", "Mean")
print(R)