function [R5p, R5v, R10z] = R5x(z, dx)
    R5p = feature_characterization(z, dx, "P", "Wolfprune 5 %","Top 5", "PVh", "Mean");
    R5v = feature_characterization(z, dx, "V", "Wolfprune 5 %", "Bot 5", "PVh", "Mean");
    R10z = R5p + R5v;
end