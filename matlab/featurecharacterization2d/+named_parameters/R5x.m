function [R5p, R5v, R10z] = R5x(z, dx, tree)
    R5p = sammtop_functions.feature_parameters.feature_based_parameters.feature_parameter(z, dx, tree, "P", "Top", 5, "PVh", "Mean");
    R5v = sammtop_functions.feature_parameters.feature_based_parameters.feature_parameter(z,  dx, tree, "V", "Bot", 5, "PVh", "Mean");
    R10z = R5p + R5v;
end