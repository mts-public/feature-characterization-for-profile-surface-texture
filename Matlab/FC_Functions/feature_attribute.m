function ATTR = feature_attribute(z, dx, M, AT)
% INPUTS:
%   z   - vertical profile values
%   dx  - step size in x-direction
%   M   - motif array
%   AT  - attribute type {"Wolfprune", "HDh", "Width", "HDw", "VolS",...
%                         "HDv", "DevLength", "HDl", "PVh", "Curvature",...
%                         "Count"}
% Outputs
%   ATTR - attribute value of given motifs

I_sig = find([M.sig] == 1);
switch AT
    case {"Wolfprune", "HDh"}
        ATTR = abs(z(floor([M(I_sig).ilp])) - z(floor([M(I_sig).iv])));
    case {"Width", "HDw"}
        for i = 1:length(I_sig)
            ATTR(i) = dx*max(abs(M(I_sig(i)).ihi - M(I_sig(i)).ilp));
        end
    case {"VolS", "HDv"}
        for i = 1:length(I_sig)
            ATTR(i) = HDvf(z, dx, M(I_sig(i)));
        end
    case {"DevLength", "HDl"}
        for i = 1:length(I_sig)
            ATTR(i) = HDlf(z, dx, M(I_sig(i)));
        end
    case "PVh"
        FTI = sign(z(floor(M(1).ilp)) - z(floor(M(1).iv)));
        ATTR = -FTI*z(floor([M(I_sig).iv]));
    case "Curvature"
        for i = 1:length(I_sig)
            ATTR(i) = curvature(z, dx, M(I_sig(i)).iv);
        end
    case "Count"
        ATTR = ones(1, length(I_sig));
end
end

function HDv = HDvf(z, dx, Mr)
% all heightintersections incl. low-peak
ihi = [Mr.ilp; Mr.ihi];
zlp = z(floor(Mr.ilp));
A = 0;i=1;
dir = sign(Mr.ihp - Mr.ilp);
while i < length(ihi)
    i1 = abs(ceil(dir*ihi(i)));     % round toward pit of current area
    i2 = abs(floor(dir*ihi(i+1)));
    xf = [ihi(i); (i1:dir:i2)'; ihi(i+1)]*dx;
    zf = [zlp; z(i1:dir:i2); zlp];
    A = A + abs(trapz(xf,zf-zlp));
    i = i + 2;
end
HDv = A/(length(z)*dx);
end

function  HDl = HDlf(z, dx, Mr)
zlp = z(floor(Mr.ilp));
dir = sign(Mr.ihp - Mr.ilp);
ihi_end = Mr.ihi(end);
i1 = abs(ceil(dir*Mr.ilp));
i2 = abs(floor(dir*ihi_end));
zf = [z(i1:dir:i2)];
HDl = (sum(sqrt(1 + (diff(zf)./dx).^2))...
    + mod(Mr.ilp, 1))*dx...
    + sqrt((ihi_end - i2)^2*dx^2+(zlp-z(i2))^2);
end

function cx = curvature(z, dx, ix)
if mod(ix,1)~=0
    ix = [floor(ix) ceil(ix)];
end
for n=1:length(ix)
i = ix(n);
switch i
    case 1
    dz1 = (-147*z(1) + 360*z(2) - 450*z(3) + 400*z(4) - 225*z(5) + ...
        72*z(6) - 10*z(7))/(60*dx);
    dz2 = (812*z(1) - 3132*z(2) + 5265*z(3) - 5080*z(4) + 2970*z(5) ...
        - 972*z(6) + 137*z(7))/(180*(dx)^2);
    case 2
    dz1 = (-10*z(1) - 77*z(2) + 150*z(3) - 100*z(4) + 50*z(5) - ...
        15*z(6) + 2*z(7))/(60*dx);
    dz2 = (137*z(1) - 147*z(2) - 255*z(3) + 470*z(4) - 285*z(5) + ...
        93*z(6) - 13*z(7))/(180*(dx)^2);
    case 3
    dz1 = (2*z(1) - 24*z(2) - 35*z(3) + 80*z(4) - 30*z(5) + 8*z(6) ...
        - 1*z(7))/(60*dx);
    dz2 = (-13*z(1) + 288*z(2) - 420*z(3) + 200*z(4) + 15*z(5) + ...
        12*z(6) + 2*z(7))/(180*(dx)^2);
    case length(z) - 3
    dz1 = (z(end-7) - 8*z(end-6)+30*z(end-5) - 80*z(end-4) + ...
        35*z(end-3) + 24*z(end-2) - 2*z(end-1))/(60*dx);
    dz2 = (2*z(end-7) - 12*z(end-6) + 15*z(end-5) + 200*z(end-4) - ...
        420*z(end-3) + 228*z(end-2) - 13*z(end-1))/(180*(dx)^2);
    case length(z) - 2
    dz1 = (-2*z(end-7) + 15*z(end-6) - 50*z(end-5) + 100*z(end-4) - ...
        150*z(end-3) + 77*z(end-2) + 10*z(end-1))/(60*dx);
    dz2 = (-13*z(end-7) + 93*z(end-6) - 285*z(end-5) + 470*z(end-4) ...
        - 255*z(end-3) - 147*z(end-2) + 137*z(end-1))/(180*(dx)^2);
    case length(z)-1
    dz1 = (10*z(end-7) - 72*z(end-6) + 225*z(end-5) - 400*z(end-4) ...
        + 450*z(end-3) - 360*z(end-2) + 147*z(end-1))/(60*dx);
    dz2 = (137*z(end-7) + 93*z(end-6) - 285*z(end-5) + 470*z(end-4) ...
        - 255*z(end-3) - 147*z(end-2) + 137*z(end-1))/(180*(dx)^2);
    otherwise
    dz1 = (-z(i-3) + 9*z(i-2) - 45*z(i-1) + 45*z(i+1) - 9*z(i+2) + ...
        z(i+3))/(60*dx);
    dz2 = (2*z(i-3) - 27*z(i-2) + 270*z(i-1) - 490*z(i) + 270*z(i+1)...
        - 27*z(i+2) + 2*z(i+3))/(180*(dx)^2);
end
cx(n) = dz2/(1 + dz1^2)^(3/2);
end
cx=mean(cx);
end

% Herausforderung hier: wie umgehen mit den interpolierten indizes?
% z.B. eine Plateau bestehend aus zwei Punkten 
%% Curvature with spline
function K = curvature_spline(z, dx, c, ik)
i0 = floor(ik);
ir = mod(ik, 1);
% other coefficients
di = (1/(3*dx))*(c(i0+1) - c(i0));
bi = (z(i0+1) - z(i0))/dx - c(i0)*dx - di*dx^2;
% first derivative
d1 = bi+2*c(i0)*dx*ir+3*di*(dx*ir)^2;
% second derivative
d2 = 2*c(i0) + 6*di*dx*ir;
% curvature
K = d2 / (1 + d1^2)^(3/2);
end