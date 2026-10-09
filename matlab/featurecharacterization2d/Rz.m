function Rz = Rz(z, dx, n_sc)
% maximum height Rz according to ISO 21920-2: mean of the sum of the
% largest peak height and the largest pit depth of each section. Peaks and
% pits are determined by an adapted crossing-the-line segmentation at
% z = 0 (zero crossings).
% INPUTS:
%   z    - vertical profile values in µm
%   dx   - step size in x-direction in mm
%   n_sc - number of sections (optional, default 5)
% OUTPUTS:
%   Rz   - maximum height in µm
if nargin < 3
    n_sc = 5; % default value
end
z = z(:);
[PF, nPF] = zero_crossing(z, dx);
% length of section
l_sc = floor(length(z)/n_sc)*dx;
% computing Rpi and Rvi
Rpi = zeros(n_sc, 1);
Rvi = zeros(n_sc, 1);
for i = 1:n_sc
    Rpsc = []; Rvsc = [];
    for j = 1:nPF
        if PF(j).xh >= (i-1)*l_sc && PF(j).xh < i*l_sc
            if PF(j).t == 1
                Rpsc = [Rpsc; PF(j).h];
            else
                Rvsc = [Rvsc; PF(j).h];
            end
        end
    end
    Rpi(i) = max([Rpsc; 0]);
    Rvi(i) = max([Rvsc; 0]);
end
Rz = mean(Rpi + Rvi);
end

%% adapted crossing-the-line algorithm
function [PF, nPF] = zero_crossing(z, dx)
% INPUTS:
%   z   - vertical profile values in µm
%   dx  - step size in x-direction in mm
% OUTPUTS:
%   PF  - profile elements with type t (1: peak, -1: pit), height h and
%         position xh of the largest absolute value
%   nPF - number of profile elements
PF = struct('t', {}, 'h', {}, 'xh', {});
n = length(z);
i = 1; j = 2; nPF = 0;
while j <= n
    if (z(j-1) <= 0 && z(j) > 0) || (z(j-1) >= 0 && z(j) < 0)
        nPF = nPF + 1;
        [~, kh] = max(abs(z(i:j-1)));
        kh = kh - 1 + i;
        PF(nPF) = struct('t', sign(z(kh)), 'h', abs(z(kh)), 'xh', (kh-1)*dx);
        i = j;
    end
    j = j + 1;
end
nPF = nPF + 1;
[~, kh] = max(abs(z(i:n)));
kh = kh - 1 + i;
PF(nPF) = struct('t', sign(z(kh)), 'h', abs(z(kh)), 'xh', (kh-1)*dx);
end
