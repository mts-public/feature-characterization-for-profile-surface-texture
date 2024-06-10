function [xFC, M, META] = feature_characterization(z, dx, FT, pruning,...
    significant, AT, stats)
% INPUTS:
%   z           - vertical profile values
%   dx          - step size in x-direction
%   FT          - 'D', 'V', 'H', 'P'
%   pruning     - 'None', 'Wolfprune TH/X%', 'Width TH/X%', 
%                 'VolS TH', 'DevLength TH'
%                 Threshold TH in units of corresponding attribute.
%                 For Wolfprune or Width, the threshold value can also be 
%                 specified as a percentage. In this case, X% of Rz is used
%                 as the threshold in the case of Wolfprune and X% of Le in
%                 the case of Width.
%   significant - 'All', 'Closed c', 'Open c', 'Bot N', 'Top N'
%                 c can be an absolute value in µm or if the value is given 
%                 as a percentage it is interpreted as a material ratio 
%                 from which the height is then determined.
%                 N specifies the number of top or bot values.
%   AT          - 'HDh', 'HDv', 'HDw', 'HDl', 'PVh', 'Curvature', 'Count'
%   stats       - 'Mean', 'Max', 'Min', 'StdDev', 'Perc X',
%                 'Hist', 'Sum', 'Density'
%                 for "Perc", x specifies the threshold in the units of the
%                 corresponding attribute
% OUTPUTS:
%   xFC         - parameter based on feature characterization
%   M           - structured array of motifs
%   META        - meta data for further processing (e.g. plotting)

%% parse pruning
pruning = strrep(pruning,"%"," %"); % add blank before "%"
str = split(pruning, " ");
PT = str(1);
N = length(str);
TH = str2double(str(min(2,N)));     % TH=NaN when no number
if N == 2 && isnan(TH)              % 
    TH = "opt";
end
if N >= 3                           % N=4 if by strrep 2 blanks
    if PT=="Wolfprune"
        TH = (TH/100)*Rz(z, dx);
    elseif PT=="Width"
        TH = (TH/100)*length(z)*dx;
    end
end

%% parse significant
significant = strrep(significant,"%"," %");
str = split(significant, " ");
Fsig = str(1);
N = length(str);
NIsig = str2double(str(min(2,N)));
if N >= 3
    NIsig = max(z) + Rcm(z, dx, NIsig);
end

%% parse stats
str = split(stats, " ");
Astats = str(1);
vstats = str2double(str(end));

%% feature characterization
M = watershed_segmentation(z, dx, FT, PT, TH);
[xFC, M, attr] = feature_parameter(z, dx, M,...
    Fsig, NIsig, AT, Astats, vstats);
META = struct('attr', attr, 'nM', length(M), 'PT', PT, 'TH', TH, 'Fsig',...
    Fsig, 'NIsig', NIsig,'AT', AT, 'Astats', Astats,'vstats', vstats);
end