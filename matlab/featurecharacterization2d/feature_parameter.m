function [xFC, M, ATTR, Fsig, NIsig] = feature_parameter(z, dx, M, Fsig, NIsig, AT, Astats, vstats)
%% step 4: determine significant_features
I_Nsig = [];
switch Fsig
    case {"Open", "Closed"}
        % feature type indicator (FTI=-1 if hills/peaks, FTI=1 if dales/pits)
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
    case {"Top", "Bot"}
        % determine attribute values
        ATTR = feature_attribute(z, dx, M, "PVh");
        % determine indices (I) of sorted zv-values in zv
        [~, I_sort] = sort(ATTR, 'descend');
        % if NIsig is higher than nM use nM
        NIsig = min(NIsig, length(M));
        I_Nsig = I_sort(NIsig+1:end);
end
% set indicator M.sig zero for not significant motifs
for i = 1:length(I_Nsig)
    M(I_Nsig(i)).sig = 0;
end

%% step 5: determine attibrute-values of significant features
ATTR = feature_attribute(z, dx, M, AT);

%% step 6: attribute statistics
switch Astats
    case "Mean"
        xFC = mean(ATTR);
    case "Max"
        xFC = max(ATTR);
    case "Min"
        xFC = min(ATTR);
    case "StdDev"
        xFC = std(ATTR);
    case "Perc"
        xFC = sum(ATTR > vstats)/nMsig;
    case "Hist"
        figure
        xFC = histogram(ATTR,length(ATTR));
    case "Sum"
        xFC = sum(ATTR);
    case "Density"
        xFC = ATTR/(dx*(length(z)-1));
    case "Median"
        xFC = median(ATTR);
    case "Span"
        xFC = max(ATTR) - min(ATTR);
    case "RMS"  
        xFC = rms(ATTR);
end
end