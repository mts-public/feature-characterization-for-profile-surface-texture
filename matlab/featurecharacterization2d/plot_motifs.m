function plot_motifs(z, dx, M, Fsig, NIsig, qualitative, changetree)
if nargin < 5
    NIsig = [];
end
if nargin < 6
    qualitative = 0;
end
%% preallocating
iend = length(z);
if isscalar(dx)
    x = (0:dx:iend*dx - dx)';
else
    x = dx(:);
    dx = x(2) - x(1);
end
% make sure that vector is colum-vector
z = z(:);                                             
if isempty(M)
    error('No features detected. Maybe check pruning configuration')
end
% feature type identifier. FTI = -1 if Hills, FTI = 1 if Dales
FTI = sign(z(floor(M(1).ilp)) - z(floor(M(1).iv)));     
%% plot settings
linewidth=0.85;
% transparicy-value for not significant motifs
alpha = 0.5;                                           
% color of motifs
col = [0 0.4470 0.7410];                               
xlim([min(x) max(x)])
ylim(1.1*[min([z; NIsig']) max([z; NIsig'])])
grid on
box on
hold on

%% fill motifs and motif "frames"
if nargin ~= 7 || ~(nargin == 7 && changetree == 1)
for k = 1:length(M)
    p{k} = patch_featureelement(z, x, dx, M(k), col);
end
for k = 1:length(M)
    Mt = M(k);
    pm{k} = plot(([Mt.ilp Mt.ilp Mt.ihp Mt.ihp]-1)*dx+x(1),...
        [z(ceil(Mt.ilp)) z(ceil(Mt.iv)) z(ceil(Mt.iv)) z(ceil(Mt.ihp))],...
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
if qualitative == 0
    xlabel('profile length / mm','FontSize',10,'FontName','Times New Roman')
    ylabel('profile height / µm','FontSize',10,'FontName','Times New Roman')
elseif qualitative == 1
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    xlabel('profile length','Interpreter','latex')
    ylabel('profile height','Interpreter','latex')
elseif qualitative == 2

elseif qualitative == 3
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
end

%% change tree
if nargin == 7 && changetree == 1
treegreen=[0.4660 0.6740 0.1880];
ip = unique([[M(:).ilp] [M(:).ihp]])';
iv = [M(:).iv];
% sort peaks descending based on their heights(z-values)
[zp, ip_s] = sort(FTI*z(floor(ip)), 'descend');    
% rearrange the peak-indices in the same order                        
ip = ip(ip_s);                                                              
% indices of peaks followed from pits
ipv = [ip(:); iv(:)];                                                       
% type: NaN for peak, 1 for motif (later replaced with height intersection)
tpv = [NaN(1, length(ip)) ones(1, length(iv))];                             
% preallocating change tree
tree = zeros(length(ip), 5);                                                
% third column equal descending peaks
tree(:, 3) = ip;                                                            
for k=1:length(ip)
    % "searchlimit" for left arm (if empty then 1)
    iliml = max([ip(ip(1:k) < ip(k) & zp(1:k) >= zp(k)); 1]);
    % index of max peak or the only pit between current peak and iliml
    ixl = find(ipv > iliml & ipv < ipv(k), 1);
    % if ixl empty then 1
    tree(k, 2) = max([1 ipv(ixl)]);
    % if ixl empty then NaN
    tree(k, 1) = max([NaN tpv(ixl)]);
    % "searchlimit" for right arm (if empty then iend)
    ilimr = min([ip(ip(1:k) > ip(k) & zp(1:k) >= zp(k)); iend]);
     % index of max peak or the only pit between current peak and ilimr
    ixr = find(ipv < ilimr & ipv > ipv(k), 1);
    % if ixr empty then length(x)
    tree(k, 4) = min([iend ipv(ixr)]);
     % if ixr empty then NaN
    tree(k, 5) = max([NaN tpv(ixr)]);                                      
end
hold on
for k = 1:size(tree,1)
% current arm
T = tree(k,:);
% left border case   
if tree(k, 2) == 1                 
    plot(...
        ([T(2) T(2) T(4) T(4)]-1)*dx,...
        [z(floor(T(3))) z(floor(T(3))) z(floor(T(3))) z(floor(T(4)))],...
        'color', treegreen)
% right border case
elseif tree(k, 4) == iend   
    plot(...
        ([T(2) T(2) T(4) T(4)]-1)*dx,...
        [z(floor(T(2))) z(floor(T(3))) z(floor(T(3))) z(floor(T(3)))],...
        'color', treegreen)
else
    plot(([dx*(T(2)) dx*(T(2)) dx*(T(4)) dx*(T(4))]-dx),...
        [z(floor(T(2))) z(floor(T(3))) z(floor(T(3))) z(floor(T(4)))],...
        'color', treegreen)
end
end
end
hold off

function p = patch_featureelement(z, x, dx, Mr, col)
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
            col, 'LineStyle', 'none');
    end
end
end