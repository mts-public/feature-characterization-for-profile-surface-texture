function attr = feature_attribute(z, dx, M, AT)
% INPUTS:
%   z   - vertical profile values
%   dx  - step size in x-direction
%   M   - motif array
%   AT  - attribute type {"Wolfprune", "HDh", "Width", "HDw", "VolS",...
%                         "HDv", "DevLength", "HDl", "PVh", "Curvature",...
%                         "Count"}
% Outputs
%   attr - attribute value of given motifs

I_sig = find([M.sig] == 1);
switch AT
    case {"Wolfprune", "HDh"}
        attr = abs(z(floor([M(I_sig).ilp])) - z(floor([M(I_sig).iv])));
    case {"Width", "HDw"}
        for i = 1:length(I_sig)
            attr(i) = dx*max(abs(M(I_sig(i)).ihi - M(I_sig(i)).ilp));
        end
    case {"VolS", "HDv"}
        for i = 1:length(I_sig)
            attr(i) = HDvf(z, dx, M(I_sig(i)));
        end
    case {"DevLength", "HDl"}
        for i = 1:length(I_sig)
            attr(i) = HDlf(z, dx, M(I_sig(i)));
        end
    case "PVh"
        FTI = sign(z(floor(M(1).ilp)) - z(floor(M(1).iv)));
        attr = -FTI*z(floor([M(I_sig).iv]));
    case "Curvature"
        for i = 1:length(I_sig)
            attr(i) = curvature(z, dx, M(I_sig(i)).iv);
        end
    case "Count"
        attr = ones(1, length(I_sig));
end
end

function HDv = HDvf(z, dx, Mr)
ihi = [Mr.ilp; Mr.ihi];
zlp = z(floor(Mr.ilp));
A = 0;i=1;
dir = sign(Mr.ihp - Mr.ilp);
while i < length(ihi)
    i1 = abs(ceil(dir*ihi(i)));
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
if mod(ix,1) ~= 0
    ix = [floor(ix) ceil(ix)];
end
for n = 1:length(ix)
i = ix(n);
switch i
    case 1
        dz1(i) = (-147*z(1) + 360*z(2) - 450*z(3) + 400*z(4) - 225*z(5) + ...
            72*z(6) - 10*z(7))/(60*dx);
        dz2(i) = (812*z(1) - 3132*z(2) + 5265*z(3) - 5080*z(4) + 2970*z(5) ...
            - 972*z(6) + 137*z(7))/(180*(dx)^2);
    case 2
        dz1(i) = (-10*z(1) - 77*z(2) + 150*z(3) - 100*z(4) + 50*z(5) - ...
            15*z(6) + 2*z(7))/(60*dx);
        dz2(i) = (137*z(1) - 147*z(2) - 255*z(3) + 470*z(4) - 285*z(5) + ...
            93*z(6) - 13*z(7))/(180*(dx)^2);
    case 3
        dz1(i) = (2*z(1) - 24*z(2) - 35*z(3) + 80*z(4) - 30*z(5) + 8*z(6) ...
            - 1*z(7))/(60*dx);
        dz2(i) = (-13*z(1) + 288*z(2) - 420*z(3) + 200*z(4) + 15*z(5) + ...
            12*z(6) + 2*z(7))/(180*(dx)^2);
    case length(z) - 2
        dz1(i) = (z(end-6) - 8*z(end-5)+30*z(end-4) - 80*z(end-3) + ...
            35*z(end-2) + 24*z(end-1) - 2*z(end))/(60*dx);
        dz2(i) = (2*z(end-6) - 12*z(end-5) + 15*z(end-4) + 200*z(end-3) - ...
            420*z(end-2) + 228*z(end-1) - 13*z(end))/(180*(dx)^2);
    case length(z) - 1
        dz1(i) = (-2*z(end-6) + 15*z(end-5) - 50*z(end-4) + 100*z(end-3) - ...
            150*z(end-2) + 77*z(end-1) + 10*z(end))/(60*dx);
        dz2(i) = (-13*z(end-6) + 93*z(end-5) - 285*z(end-4) + 470*z(end-3) ...
            - 255*z(end-2) - 147*z(end-1) + 137*z(end))/(180*(dx)^2);
    case length(z)
        dz1(i) = (10*z(end-6) - 72*z(end-5) + 225*z(end-4) - 400*z(end-3) ...
            + 450*z(end-2) - 360*z(end-1) + 147*z(end))/(60*dx);
        dz2(i) = (137*z(end-6) + 93*z(end-5) - 285*z(end-4) + 470*z(end-3) ...
            - 255*z(end-2) - 147*z(end-1) + 137*z(end))/(180*(dx)^2);
    otherwise
        dz1(i) = (-z(i-3) + 9*z(i-2) - 45*z(i-1) + 45*z(i+1) - 9*z(i+2) + ...
            z(i+3))/(60*dx);
        dz2(i) = (2*z(i-3) - 27*z(i-2) + 270*z(i-1) - 490*z(i) + 270*z(i+1)...
            - 27*z(i+2) + 2*z(i+3))/(180*(dx)^2);
end
cx(n) = dz2/(1 + dz1^2)^(3/2);
end
cx = mean(cx);
end