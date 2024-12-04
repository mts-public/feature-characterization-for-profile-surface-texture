function TH = TH_for_optimal_periodicity(z, dx, FT, PT)
% INPUTS:
%  z  - vertical profile values
%  dx - step size in x-direction
%  FT - feature type: {'D', 'V', 'H', 'P'}
%  PT - pruning type: {'None', 'Wolfprune', 'Width', 'VolS', 'DevLength'}
% OUTPUTS:
%  M    - structure array with motifs with four members 
%         (referring to Dale-motif):
%         M.iv    - (interpolated) index of pit
%         M.ilp   - (interpolated) index of low-peak
%         M.ihp   - (interpolated) index of high-peak
%         M.ihi   - (interpolated) index of heightintersection
%         M.sig   - indicator for significant features

%% step 1: determine indices of all peaks and pits
% invert z-values if hill-motifs are searched. Allows the rest of the code
% refers to dale-motifs
if FT == 'H' || FT == 'P'
    z = -z;
end
% consider only the first of any adjacent pairs of equal values so that
% plateaus are single points
iNeq = [1; 1 + find(z(1:end-1) ~= z(2:end))];
% determine the slope for each point. s=-1 neg. slope; s=1 pos. slope
s = sign(diff(z(iNeq)));
% change from 1 to -1 corresponds to a peak. -1 to 1 to pit
ipNeq = 1 + find(diff(s) == -2);
ivNeq = 1 + find(diff(s) == 2);
% peaks and pits indices into the original index vector
ipv = [iNeq(ipNeq); iNeq(ivNeq)];
% examine peaks and pits that are plateaus
for j = find(z(ipv) == z(ipv+1))'
    % number of equal high values on the plateau
    n_plateau = find(z(ipv(j):end) ~= z(ipv(j)), 1) - 1;
    % replace index with (interpolated) index of middle of plateau
    ipv(j) = ipv(j) + (n_plateau - 1)/2;                      
end
% (interpolated) indices of peaks and pits (plateaus taken into account)
ip = ipv(1:length(ipNeq));
iv = ipv(length(ipNeq) + 1:end);

%% step 2: determine motifs
% just keep indices of pits that are enclosed with peaks
iv = iv(iv > ip(1) & iv < ip(end));
% enrich structure array M with information for each motif such as pit
% (iv), low-peak (lp) and high-peak (hp) and the height intersection (ihi)
nM = length(iv);
for k = 1:nM
    [ilp, ihp] = get_ilp_ihp(z, [ip(k), ip(k+1)]);
    ihi = height_intersections(z, ilp, ihp);
    M(k) = struct('iv', iv(k), 'ilp', ilp, 'ihp', ihp, 'ihi', ihi,'sig',1);
end

%% step 3: pruning (see pruning cases in readme.md)
% determine attribute-values of each motif
attr = feature_attribute(z, dx, M, PT);
% minimal Q-Value
Qmin = 3;
% set default threshold for the case that Qmin is never exceeded
TH = 100;
% prune until just two motifs are left
while nM > 3
    % parameter Q as a measure for periodicity
    Q = mean(attr) / std(attr);
    % if Q is greater than Qmin then overwrite Qmin with the current 
    % Q-value and TH with minimal attribute value
    if Q > Qmin
        TH = min(attr);
        Qmin = Q;
    end
    % Inline prune_min_motif steps
    % row-index of minimal attribute value
    [~, rmin] = min(attr);
    % save motif with minimal attribute-value temporarily, delete 
    % corresponding entry in motif-array and attribute-vector. update nM
    Mmin = M(rmin);
    M(rmin) = []; attr(rmin) = []; 
    nM = nM - 1;
    % determine row-index of motif which is to update (rU). dir=-1: left of
    % min-motif. dir=1 right of min-motif. rU = rmin if dir=1 because
    % min-motif was deleted in motif array. else rU = rmin - 1.
    dir = sign(Mmin.ilp - Mmin.iv);
    rU = rmin - (dir == -1);
    % case 1: if in that direction is border no further steps are required
    if rU == 0 || rU > nM
        continue
    end
    % case 2: if low-peak of min-motif and motif to update is the same then 
    % determine low-peak and high-peak of motif to update. Update height
    % intersection and attribute-value (see behind if-query)
    if M(rU).ilp == Mmin.ilp
        [M(rU).ilp, M(rU).ihp] = get_ilp_ihp(z, [M(rU).ihp, Mmin.ihp]);
    % case 3: low-peak of min-motif is equal the high-peak of motif to
    % update. replace high-peak of motif to merge with high-peak of
    % min-motif
    else
        M(rU).ihp = Mmin.ihp;
        % case 3.1: if low-peak of motif to merge is lower or equal than
        % pit of min-motif no further steps are required
        if z(floor(M(rU).ilp)) <= z(floor(Mmin.iv))
            continue
        end
    end
    % case 3.2: low-peak of motif to merge is higher than pit of
    % min-motif: refresh height intersection and attribute-value
    % (see behind if-query)
    % update height intersection and attribute-value of motif to update for
    % case 2 and 3.2
    M(rU).ihi = height_intersections(z, M(rU).ilp, M(rU).ihp);
    attr(rU) = feature_attribute(z, dx, M(rU), PT);
    % End of prune_min_motif steps inlined
end
end
%% height intersection function
function ihi = height_intersections(z, ilp, ihp)
% INPUTS:
%   z    - vertical profile values
%   ilp  - (interpolated) index of low-peak
%   ihp  - (interpolated) index of high-peak
% OUTPUTS:
%   ihi  - (interpolated) index of height intersection

% direction in which to search ihi outgoing form low-peak
dir = sign(ihp - ilp);
ilp = round(ilp);ihp = round(ihp);
zlp = z(ilp);
ihi = [];
% starting index. if plateau: index of edge of plateau
j = ilp+dir*(find(z(ilp:dir:round(ihp))~=zlp,1)-1);
% getting height-intesections (based on crossing-the-line-segmentation)
while j ~= ihp
    if (z(j) < zlp && z(j+dir) >= zlp) || (z(j) >= zlp && z(j+dir) < zlp)
        ihi = [ihi; j + dir*(zlp - z(j))/(z(j + dir) - z(j))];
    end
    j = j + dir;
end
end

%% determine indices of low-peak and high-peak based of 2 given indices
function [ilp, ihp] = get_ilp_ihp(z, ip_surr)
% INPUTS:
%   z       - vertical profile values
%   ip_surr - (interpolated) indices of the two surrounding peaks of examined pit
% OUTPUTS:
%   ilp     - (interpolated) index of low-peak
%   ihp     - (interpolated) index of high-peak
[~, I] = min([z(floor(ip_surr(1))), z(floor(ip_surr(2)))]);
ilp  = ip_surr(I);
ihp = ip_surr(3 - I);
end

%% reduced find function (for translation to other languages)
function iN0 = find(A, firstOnly) 
% INPUTS
%   A         - Array 
%   firstOnly - optional: if set to 1 only first index of nonzeroelement is
%               returned. if set to 0 or not provided all indices of 
%               nonzeroelements are returned
% OUTPUTS
%   iN0       - indices of nonzero entries in A
if nargin < 2
    firstOnly = 0;
end
iN0 = [];
for i = 1:length(A)
    if A(i) ~= 0
        iN0 = [iN0; i];
        if firstOnly == 1
            return
        end
    end
end
end