clc;clear;close all;

% simple self-generated profile to illustrate handling and functionality
z = [3.3 2 1 5 3.8 4 1.5 1.5 3.5 2.5 2 -1 0 3 1.2 2 -1.2 -5 -4 -4.5 -2 -2.3 1 3 3 3 4 4.5 4.5 4 1.5 1.5 3.5 4 9 8 -1 -1 -1 -1 7 7 7 0 0.5 1 5 2 6 4.5 0.5 1 2 -1 0 3 5.2 5 5.5 4 7];
z = (z-mean(z))';
dx = 0.5;
% load("case studies\profile\Bu_1_56_ak.mat")


tic
[xFC, M, meta] = feature_characterization(z, dx, "D", "Wolfprune 4.1", "Top 3", "HDv", "Mean")
toc
figure(1)
plot_motifs(z, dx, M, meta.Fsig, meta.NIsig);

% %% Feature characterization
% feature_type = 'D';                 % 'D' or 'V' for dale tree; 'H' or 'P' for hill tree
% pruning_type = 'Wolfprune';         % prunging type: 'Wolfprune', 'Width', 'VolS' or 'DevLength'
% TH_pruning = 0.71;                  % Treshold for pruning
% Fsig = "All";                       % signicant Feature "All", "Open", "Closed", "Top", "Bot"
% NIsig = -1;                         % Nesting index for significant features
% ATsig = [];
% attribute ="PVh";                   % "HDh","HDv","HDw","HDl","PVh","Curvature" or "Count"
% attribute_stats = "Mean";           % "Mean","Median","Max","Min","StdDev","Perc","Hist","Sum","Density"
% vstats = 1;                         % extra value for Hist or Perc
% 
% tic
% M = based_on_feature_characterization.segmentation.watershed_segmentation(z, dx, feature_type, pruning_type, TH_pruning);
% [fp, M, ATTR, Fsig, NIsig] = based_on_feature_characterization.feature_parameter(z, dx, M, Fsig, NIsig, ATsig, attribute, attribute_stats);
% toc
% figure(2)
% based_on_feature_characterization.visualisation.plot_motifs(z, dx, M, Fsig, NIsig);