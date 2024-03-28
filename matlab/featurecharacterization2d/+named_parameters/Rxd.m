function [Rpd, Rvd] = Rxd(z,dx,tree)
Rpd = feature_characterization(z, dx, tree, "P", "Wolfprune 5 %", "All", "Count", "Density");
Rvd = feature_characterization(z, dx, tree, "V", "Wolfprune 5 %", "All", "Count", "Density");
end