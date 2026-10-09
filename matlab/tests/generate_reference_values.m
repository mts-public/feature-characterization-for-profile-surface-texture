% Generates MATLAB reference values for the regression tests of the Python
% package (python/tests). Run from any folder:
%   matlab -batch "run('matlab/tests/generate_reference_values.m')"
% Outputs:
%   python/tests/data/matlab_reference_fc.csv     - feature characterization
%   python/tests/data/matlab_reference_rz_rcm.csv - Rz and Rcm
%   python/tests/data/matlab_reference_default.csv - default_FC_parameters
%   python/tests/data/matlab_reference_smd.csv     - smd2mat
clear; warning('off', 'all');

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'matlab', 'featurecharacterization2d'));
addpath(fullfile(root, 'matlab', 'softgauge'));
outdir = fullfile(root, 'python', 'tests', 'data');
if ~exist(outdir, 'dir'), mkdir(outdir); end

profiles = ["Bu_1_56_ak", "Luftpresser_1_56_ak"];
dx = 0.5e-3; % mm

% FT, pruning, significant, AT, stats
configs = [
    "D", "None",           "All",         "HDh",       "Mean"
    "D", "Wolfprune 3",    "All",         "HDh",       "Mean"
    "D", "VolS 1",         "All",         "HDh",       "Mean"
    "D", "Width 0.05",     "All",         "HDh",       "Mean"
    "D", "DevLength 0.05", "All",         "HDh",       "Mean"
    "D", "Wolfprune 5 %",  "All",         "HDh",       "Mean"
    "D", "Width 5 %",      "All",         "HDw",       "Mean"
    "D", "Wolfprune opt",  "All",         "HDh",       "Mean"
    "D", "Wolfprune 1000", "All",         "HDh",       "Mean"
    "P", "Wolfprune 5 %",  "Top 5",       "PVh",       "Mean"
    "V", "Wolfprune 5 %",  "Bot 5",       "PVh",       "Mean"
    "D", "Wolfprune 5 %",  "Open 0",      "HDh",       "Mean"
    "H", "Wolfprune 5 %",  "Open 0",      "HDh",       "Mean"
    "D", "Wolfprune 5 %",  "Closed 50 %", "HDh",       "Mean"
    "D", "Wolfprune 5 %",  "Closed 0 %",  "HDh",       "Mean"
    "H", "Wolfprune 5 %",  "Closed 50 %", "HDv",       "Mean"
    "D", "Wolfprune 5 %",  "All",         "HDv",       "Sum"
    "D", "Wolfprune 5 %",  "All",         "HDl",       "Max"
    "D", "Wolfprune 5 %",  "All",         "HDw",       "Min"
    "D", "Wolfprune 5 %",  "All",         "HDh",       "StdDev"
    "D", "Wolfprune 5 %",  "All",         "HDh",       "Perc 2"
    "V", "Wolfprune 5 %",  "All",         "Curvature", "Mean"
    "P", "Wolfprune 5 %",  "All",         "Curvature", "Mean"
    "P", "Wolfprune 5 %",  "All",         "Count",     "Density"
];

fid = fopen(fullfile(outdir, 'matlab_reference_fc.csv'), 'w');
fprintf(fid, 'profile,FT,pruning,significant,AT,stats,xFC,nM\n');
for p = profiles
    S = load(fullfile(root, 'data', 'profiles', p + ".mat"));
    z = S.z - mean(S.z);
    for k = 1:size(configs, 1)
        c = configs(k, :);
        [xFC, ~, META] = feature_characterization(z, dx, c(1), c(2), c(3), ...
            c(4), c(5));
        fprintf(fid, '%s,%s,%s,%s,%s,%s,%.17g,%d\n', p, c, xFC, META.nM);
    end
end
fclose(fid);

fid = fopen(fullfile(outdir, 'matlab_reference_rz_rcm.csv'), 'w');
fprintf(fid, 'profile,quantity,p,value\n');
for p = profiles
    S = load(fullfile(root, 'data', 'profiles', p + ".mat"));
    z = S.z - mean(S.z);
    fprintf(fid, '%s,Rz,,%.17g\n', p, Rz(z, dx));
    n = length(z);
    % incl. material ratios exactly at sampling points k/n (k = 1, n/2, n-1)
    for mr = [0 10 50 90 100 [1 floor(n/2) n-1]/n*100]
        fprintf(fid, '%s,Rcm,%.17g,%.17g\n', p, mr, Rcm(z, mr));
    end
end
fclose(fid);

% named feature parameters with default settings (ISO 21920-2/-3)
names = ["Rpd", "Rvd", "Rmpc", "Rmvc", "R5p", "R5v", "R10z"];
fid = fopen(fullfile(outdir, 'matlab_reference_default.csv'), 'w');
fprintf(fid, 'profile,parameter,value\n');
for p = profiles
    S = load(fullfile(root, 'data', 'profiles', p + ".mat"));
    z = S.z - mean(S.z);
    xFC = default_FC_parameters(z, dx);
    for name = names
        fprintf(fid, '%s,%s,%.17g\n', p, name, xFC.(name));
    end
end
fclose(fid);

% reading of softgauge files (*.smd)
files = dir(fullfile(root, 'data', 'profiles', '*.smd'));
fid = fopen(fullfile(outdir, 'matlab_reference_smd.csv'), 'w');
fprintf(fid, 'file,n,dx,L,x_end,z_first,z_last,z_sum\n');
for k = 1:numel(files)
    [z, L, x, dxs] = smd2mat(fullfile(files(k).folder, files(k).name));
    fprintf(fid, '%s,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n', files(k).name, ...
        numel(z), dxs, L, x(end), z(1), z(end), sum(z));
end
fclose(fid);
