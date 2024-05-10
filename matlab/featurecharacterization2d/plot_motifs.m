function plot_motifs(z, dx, M, Fsig, NIsig, qualitative, changetree)
%% preallocating
iend = length(z);
x = (0:dx:iend*dx - dx)'/1000;
z = z(:);                                              % make sure that vector is colum-vector
FT = sign(z(floor(M(1).ilp)) - z(floor(M(1).iv)));     % identify feature type. FT = -1 if Hills, FT = 1 if Dales
%% plot settings
xlabel('profile length / mm','FontSize',10,'FontName','Times New Roman')
ylabel('profile height / µm','FontSize',10,'FontName','Times New Roman')
linewidth=0.85;
alpha = 0.5;                                           % transparicy-value for not significant motifs
col = [0 0.4470 0.7410];                               % color of motifs
xlim([0 max(x)])
ylim(1.1*[min(z) max(z)])
grid on
box on
hold on

%% fill motifs and motif "frames"
if nargin ~= 7 || ~(nargin == 7 && changetree == 1)
for k = 1:length(M)
    p{k} = patch_featureelement(z, dx, M(k), col);
end

for k = 1:length(M)
    pm{k} = plot((dx*[(M(k).ilp) (M(k).ilp) (M(k).ihp) (M(k).ihp)]-dx)/1000,...
        [z(floor(M(k).ilp)) z(floor(M(k).iv)) z(floor(M(k).iv)) z(floor(M(k).ihp))],...
        'color', [1 0 0],'LineStyle','-');
end
end

%% plot
hold on
plot(x, z, 'k', 'LineWidth', linewidth)

%% not significant features
for n = find([M.sig]==0)
    for i = 1:length(p{1,n})             
        p{1,n}{i}.FaceAlpha = alpha;
        pm{n}.LineStyle =':';
    end
end
if nargin > 3
    switch Fsig
        case "Open"
            yline(NIsig, '--', 'LineWidth', 1);
        case "Closed"
            yline(NIsig, '--', 'LineWidth', 1); 
        case "Top"
    end
end

%% qualitativ
if nargin >= 6 && qualitative == 1
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    xlabel('profile length','Interpreter','latex')
    ylabel('profile height','Interpreter','latex')
end

%% change tree
if nargin == 7 && changetree == 1
treegreen=[0.4660 0.6740 0.1880];
ip = unique([[M(:).ilp] [M(:).ihp]])';
iv = [M(:).iv];
[zp, ip_s] = sort(FT*z(floor(ip)), 'descend');                              % sort peaks descending based on their heights(z-values)
ip = ip(ip_s);                                                              % rearrange the peak-indices in the same order
ipv = [ip(:); iv(:)];                                                       % indices of peaks followed from pits
tpv = [NaN(1, length(ip)) ones(1, length(iv))];                             % type: NaN for peak, 1 for motif (later replaced with respective height intersection)
tree = zeros(length(ip), 5);                                                % preallocating change tree
tree(:, 3) = ip;                                                            % third column equal descending peaks
for k=1:length(ip)
    % left arm
    iliml = max([ip(ip(1:k) < ip(k) & zp(1:k) >= zp(k)); 1]);               % "searchlimit" for left arm (if empty then 1)
    ixl = find(ipv > iliml & ipv < ipv(k), 1);                              % index of max peak or the only pit between current peak and "searchlimit".
    tree(k, 2) = max([1 ipv(ixl)]);                                         % if ixl empty then 1
    tree(k, 1) = max([NaN tpv(ixl)]);                                       % if ixl empty then NaN
    % right arm
    ilimr = min([ip(ip(1:k) > ip(k) & zp(1:k) >= zp(k)); iend]);            % "searchlimit" for right arm (if empty then iend)
    ixr = find(ipv < ilimr & ipv > ipv(k), 1);                              % index of max peak or the only pit between current peak and "searchlimit"
    tree(k, 4) = min([iend ipv(ixr)]);                                      % if ixr empty then length(x)
    tree(k, 5) = max([NaN tpv(ixr)]);                                       % if ixr empty then NaN
end
hold on
for k = 1:size(tree,1)
if tree(k,2) == 1           % left border case         
    plot(([dx*(tree(k, 2)) dx*(tree(k, 2)) dx*(tree(k, 4)) dx*(tree(k, 4))]-dx)/1000,...
        [z(floor(tree(k, 3))) z(floor(tree(k, 3))) z(floor(tree(k, 3))) z(floor(tree(k, 4)))],...
        'color', treegreen)
elseif tree(k,4) == iend    % right border case
    plot(([dx*(tree(k, 2)) dx*(tree(k, 2)) dx*(tree(k, 4)) dx*(tree(k, 4))]-dx)/1000,...
        [z(floor(tree(k, 2))) z(floor(tree(k, 3))) z(floor(tree(k, 3))) z(floor(tree(k, 3)))],...
        'color', treegreen)
else
    plot(([dx*(tree(k, 2)) dx*(tree(k, 2)) dx*(tree(k, 4)) dx*(tree(k, 4))]-dx)/1000,...
        [z(floor(tree(k, 2))) z(floor(tree(k, 3))) z(floor(tree(k, 3))) z(floor(tree(k, 4)))],...
        'color', treegreen)
end

end
end
hold off

function p = patch_featureelement(z, dx, Mr, col)
    dir = sign(Mr.ihp - Mr.ilp);
    ihi = [Mr.ilp; Mr.ihi];
    zlp = z(floor(Mr.ilp));
    for i=1:2:length(ihi)-1
        i1 = abs(ceil(dir*ihi(i)));     % round toward pit of current area
        i2 = abs(floor(dir*ihi(i+1)));
        xf = dx*([ihi(i); (i1:dir:i2)'; ihi(i+1)] - 1)/1000;
        zf = [zlp; z(i1:dir:i2); zlp];
        p{ceil(i/2)} = patch([xf; flip(xf)], [zf; ones(length(zf), 1)*zlp], col, 'LineStyle', 'none');
    end
end
end


