clc;clear;close all;

% simple self-generated profile to illustrate handling and functionality
z = [3.3 2 1 5 3.8 4 1.5 1.5 3.5 2.5 2 -1 0 3 1.2 2 -1.2 -5 -4 -4.5 -2 -2.3 1 3 3 3 4 4.5 4.5 4 1.5 1.5 3.5 4 9 8 -1 -1 -1 -1 7 7 7 0 0.5 3 5 4 5 4.5 0.5 1 2 -1 0 3 5.2 5 5.5 4 7];
z = (z-mean(z))';
dx = 0.5;

[xFC, M, META] = feature_characterization(z, dx, "D", "Wolfprune 5 %", "Open 10 %", "HDh", "Perc 2.5");
plot_motifs(z, dx, M, META.Fsig, META.NIsig);