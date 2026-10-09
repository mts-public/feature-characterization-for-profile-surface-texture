function [Rcm, material, hintersection] = Rcm(z, p)
% inverse material ratio Rcm(p) in µm according to ISO 21920-2, 4.5.1.4:
% level of intersection at the material ratio p in % relative to the
% maximum height (Rcm(0) = 0). The material ratio curve is given by the
% pairs (k/n, c_k) of the heights c_k sorted in descending order
% (Annex C), between them the next value c_k is used (as the ISO21920 dll).
% INPUTS:
%   z             - vertical profile values in µm
%   p             - material ratio in % (scalar or array)
% OUTPUTS:
%   Rcm           - inverse material ratio in µm
%   material      - material ratio k/n in % of the material ratio curve
%   hintersection - heights c_k of the material ratio curve in µm
if any(p < 0 | p > 100)
    error("The material ratio p has to be within 0 % and 100 %.")
end
n = length(z);
% material ratio curve (Abbott Firestone Curve)
hintersection = sort(z(:), 'descend');
material = (1:n)'/n*100;
% index of the smallest material ratio k/n >= p
% (tolerance for rounding errors of p*n/100 at the sampling points)
k = max(1, ceil(p*n/100 - 1e-9));
Rcm = reshape(hintersection(k) - hintersection(1), size(p));
end
