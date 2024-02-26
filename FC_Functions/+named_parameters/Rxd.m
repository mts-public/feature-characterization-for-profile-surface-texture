function [Rpd, Rvd] = Rxd(z,dx,tree)
Rpd = sammtop_functions.feature_parameters.feature_based_parameters.feature_parameter(z, dx, tree, "P", "All", 1, "Count", "Density");
Rvd = sammtop_functions.feature_parameters.feature_based_parameters.feature_parameter(z, dx, tree, "V", "All", 1, "Count", "Density");
end