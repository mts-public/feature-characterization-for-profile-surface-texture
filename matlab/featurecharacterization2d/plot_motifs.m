function plot_motifs(z, dx, M, Fsig, NIsig)
% INPUTS:
%   z      - vertical profile values in µm
%   dx     - step size in x-direction in mm
%   Fsig   - significant features {”All”, ”Open”, ”Closed”, ”Top”, ”Bot”}
%   NIsig  - nesting index for significant features

%% Initialization
if nargin < 5, NIsig = []; end
if isempty(M)
    error('No features detected. Maybe check pruning configuration')
end
if isscalar(dx)
    x = (0:dx:length(z)*dx - dx)';
else
    x = dx(:);
    dx = x(2) - x(1);
end
z = z(:); % ensure column vector


%% plot settings
linewidth=0.85;
% transparicy-value for not significant motifs
alpha = 0.5;   
% linestyle for motifs frames
frame_line_style = {'-',':'};
% color of motifs
col = [0 0.4470 0.7410];
% plot limits
xlim([min(x) max(x)])
ylim(1.1*[min([z; NIsig]) max([z; NIsig])])
% axis labels
xlabel('profile length / mm','FontSize',10,'FontName','Times New Roman')
ylabel('profile height / µm','FontSize',10,'FontName','Times New Roman')
grid on
box on
hold on

%% fill motifs and motif "frames"
for k = 1:length(M)
    Mt = M(k);
    pm{k} = plot(([Mt.ilp Mt.ilp Mt.ihp Mt.ihp]-1)*dx+x(1),...
        [z(ceil(Mt.ilp)) z(ceil(Mt.iv)) z(ceil(Mt.iv)) z(ceil(Mt.ihp))],...
        'color', [1 0 0],'LineStyle', frame_line_style((Mt.sig ~= 1) + 1));
    p{k} = patch_featureelement(z, x, dx, Mt, col, alpha);
end

%% plot
plot(x, z, 'k', 'LineWidth', linewidth)

%% Threshold for Fsig = "Open" or "Closed"
if nargin > 3
    switch Fsig
        case "Open"
            yline(NIsig, '--', 'LineWidth', 1);
        case "Closed"
            yline(NIsig, '--', 'LineWidth', 1); 
    end
end
end

function p = patch_featureelement(z, x, dx, Mr, col, alpha)
    dir = sign(Mr.ihp - Mr.ilp);
    ihi = [Mr.ilp; Mr.ihi];
    zlp = z(floor(Mr.ilp));
    for i = 1:2:length(ihi)-1
        i1 = abs(ceil(dir*ihi(i)));
        i2 = abs(floor(dir*ihi(i+1)));
        xf = ([(ihi(i)-1)*dx; x(i1:dir:i2); (ihi(i+1)-1)*dx]);
        zf = [zlp; z(i1:dir:i2); zlp];
        p{ceil(i/2)} = patch([xf; flip(xf)],...
            [zf; ones(length(zf), 1)*zlp],...
            col, 'LineStyle', 'none','FaceAlpha', ~Mr.sig*alpha + Mr.sig);
    end
end