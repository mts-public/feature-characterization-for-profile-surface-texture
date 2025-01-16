function [xFC, M, attr] = feature_parameter(z, dx, M,...
                                           Fsig, NIsig, AT, Astats, vstats)
% INPUTS:
%   z      - vertical profile values in µm
%   dx     - step size in x-direction in mm
%   M      - motif array
%   Fsig   - significant features {”All”, ”Open”, ”Closed”, ”Top”, ”Bot”}
%   NIsig  - nesting index for significant features
%   AT     - attribute type {”HDh”, ”HDw”, ”HDv”, ”HDl”, ”PVh”, 
%                           ”Curvature”, ”Count”}
%   Astats - attribute statistics {”Mean”, ”Max”, ”Min”, ”StdDev”, ”Perc”,
%                                 ”Hist”, ”Sum”, ”Density” }
%   vstats - required value if Astats = ”Perc”
% OUTPUTS:
%   xFC    - feature parameter
%   M      - motif array with possibly changed significance flags
%   attr   - attribute values of significant motifs

%% step 4: determine significant_features
I_Nsig = [];
nM = length(M);
% error handling if M is empty
if isempty(M)
    warning("No features detected. Check pruning configuration.")
    xFC = NaN; attr = NaN;
    return
end
switch Fsig
    case {"Top", "Bot"}
        % if NIsig is higher than nM use nM
        NIsig = min(NIsig, nM);
        % determine attribute values
        attr = feature_attribute(z, dx, M, "PVh");
        % determine indices (I) of sorted zv-values in zv
        [~, I_sort] = sort(attr, 'descend');
        I_Nsig = I_sort(NIsig+1:end);
    case {"Open", "Closed"}
        % feature type indicator (FTI=-1: hills/peaks, FTI=1: dales/pits)
        FTI = sign(z(floor(M(1).ilp)) - z(floor(M(1).iv)));
        % determine z-values of low-peaks and pits
        zlp = z(floor([M.ilp]));
        % determine indices of not significant features
        if Fsig == "Open"
            I_Nsig = find(FTI*zlp > FTI*NIsig);
        else
            zv = z(floor([M.iv]));
            I_Nsig = find(FTI*zlp < FTI*NIsig | FTI*zv > FTI*NIsig);
        end
end
% set indicator M.sig zero for not significant motifs
for i = 1:length(I_Nsig)
    M(I_Nsig(i)).sig = 0;
end
% error handling if there are no significant features
if length(I_Nsig) == nM
    warning("All features are declared as not significant.")
    xFC = NaN; attr = NaN;
    return
end

%% step 5: determine attibrute-values of significant features
attr = feature_attribute(z, dx, M, AT);

%% step 6: attribute statistics
switch Astats
    case "Mean"
        xFC = mean(attr);
    case "Max"
        xFC = max(attr);
    case "Min"
        xFC = min(attr);
    case "StdDev"
        xFC = std(attr);
    case "Perc"
        xFC = sum(attr > vstats)/length(attr);
    case "Hist"
        figure
        xFC = histogram(attr,length(attr));
    case "Sum"
        xFC = sum(attr);
    case "Density"
        xFC = sum(attr)/(dx*length(z)/10); % unit: 1/cm
 end
end
