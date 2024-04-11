function [Rmvc,Rmpc] = Rmxc(z,dx)
Rmvc = feature_characterization(z, dx, "V", "Wolfprune 5 %", "All", "Curvature", "Mean");
Rmpc = feature_characterization(z, dx, "P", "Wolfprune 5 %", "All", "Curvature", "Mean");
end