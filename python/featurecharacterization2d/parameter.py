import numpy as np
from .attribute import FeatureAttribute


def feature_parameter(z = None, dx = None, M = None, Fsig = None, NIsig = None, AT = None, Astats = None, xstats = None): 
    ## step 4: determine significant_features
    I_Nsig = []
    match Fsig:
        case 'Open' | 'Closed':
            # feature type indicator (FTI=-1 if hills/peaks, FTI=1 if dales/pits)
            FTI = np.sign(z(int(np.floor(M[1].ilp))) - z(int(np.floor(M[1].iv))))
            # determine z-values of low-peaks and pits
            zlp = z(int(np.floor(np.array([M.ilp]))))
            # determine indices of not significant features
            if Fsig == 'Open':
                I_Nsig = find(FTI * zlp > FTI * NIsig)
            else:
                zv = z(int(np.floor(np.array([M.iv]))))
                I_Nsig = find(FTI * zlp < np.logical_or(FTI * NIsig,FTI * zv) > FTI * NIsig)
        case 'Top' | 'Bot':
            # determine attribute values
            ATTR = feature_attribute(z,dx,M,'PVh')
            # determine indices (I) of sorted zv-values in zv
            __,I_sort = __builtint__.sorted(ATTR,'descend')
            # case NIsig is higher than nM use nM
            NIsig = np.amin(NIsig,len(M))
            I_Nsig = I_sort(np.arange(NIsig + 1,end()+1))
    
    # set indicator M.sig zero for not significant motifs
    for i in np.arange(1,len(I_Nsig)+1).reshape(-1):
        M(I_Nsig[i]).sig = 0
    
    ## step 5: determine attibrute-values of significant features
    ATTR = FeatureAttribute.compute(z,dx,M,AT)
    ## step 6: attribute statistics
    match Astats:
        case 'Mean':
            xFC = np.mean(ATTR)
        case 'Max':
            xFC = np.amax(ATTR)
        case 'Min':
            xFC = np.amin(ATTR)
        case 'StdDev':
            xFC = np.std(ATTR)
        case 'Perc':
            xFC = np.sum(ATTR > xstats) / nMsig
        case 'Hist':
            xFC = histogram(ATTR,xstats)
        case 'Sum':
            xFC = sum(ATTR)
        case 'Density':
            xFC = ATTR / (dx * (len(z) - 1))
        case 'Median':
            xFC = median(ATTR)
        case 'Span':
            xFC = np.amax(ATTR) - np.amin(ATTR)
        case 'RMS':
            xFC = rms(ATTR)
    
    return xFC,M,ATTR,Fsig,NIsig