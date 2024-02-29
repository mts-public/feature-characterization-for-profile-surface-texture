clc;clear;close all;
%% simple self-generated profile to illustrate handling and functionality
z = [3.3 2 1 5 3.8 4 1.5 1.5 3.5 2.5 2 -1 0 3 1.2 2 -1.2 -5 -4 -4.5 -2 -2.3 1 3 3 3 4 4.5 4.5 4 1.5 1.5 3.5 4 9 8 -1 -1 -1 -1 7 7 7 0 0.5 1 5 2 6 4.5 0.5 1 2 -1 0 3 5.2 5 5.5 4 7]';
z = z - mean(z);
dx = 0.5;

%% feature characterization
[xFC, M, meta] = feature_characterization(z, dx, "H", "Wolfprune 1", "All", "HDh", "Mean");

%% plot
plot_motifs(z, dx, M, meta.Fsig, meta.NIsig);