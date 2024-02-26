function [Rmvc,Rmpc] = Rmxc(z,dx,tree)
Rmvc = sammtop_functions.feature_parameters.feature_based_parameters.feature_parameter(z, dx, tree, "V", "All", 1, "Curvature", "Mean");
Rmpc = sammtop_functions.feature_parameters.feature_based_parameters.feature_parameter(z, dx, tree, "P", "All", 1, "Curvature", "Mean");
end