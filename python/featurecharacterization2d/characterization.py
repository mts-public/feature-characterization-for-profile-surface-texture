import numpy as np

from .watershed import Watershed
from .parameter import feature_parameter


def feature_characterization(z = None, dx = None, FT = None, pruning = None, significant = None, AT = None, stats = None): 
    """
    Parameters
    ----------
        z : nd.array, float
            vertical profile values
        dx : float
            step size in x-direction
        FT : str
            feature type: {'D', 'V', 'H', 'P'}
        pruning : str
            'None', 'Wolfprune TH/X#', 'Width TH/X#', 'VolS TH', 'DevLength TH'
            Threshold TH in units of corresponding attribute.
            For Wolfprune or Width, the threshold value can also be
            specified as a percentage. In this case, X# of Rz is used
            as the threshold in the case of Wolfprune and X# of Le in
            the case of Width.
        significant : str
            All', 'Closed c', 'Open c', 'Bot N', 'Top N'
            c can be an absolute value in µm or if the value is given
            as a percentage it is interpreted as a material ratio
            from which the height is then determined.
            N specifies the number of top or bot values.
        AT : str
            attribute type {'HDh', 'HDv', 'HDw', 'HDl', 'PVh', 'Curvature', 'Count'}
        stats : str
            'Mean', 'Max', 'Min', 'StdDev', 'Perc X', 'Hist X', 'Sum', 'Density'
            for "Hist", x specifies the number of bins in the histogram.
            for "Perc", x specifies the threshold in the units of the
            corresponding attribute
    Returns
    -------
        xFC : np.ndarray, float
            parameter based on feature characterization
        M : Motif
            structured array of motifs
        meta :
            meta data for further processing (e.g. plotting)
    """

    ## parse pruning
    pruning = pruning.replace('%',' %')
    str = pruning.split(' ')
    PT = str[0]
    N = len(str)
    try:
        TH = float(str[np.min([2,N])-1])
    except ValueError:
        TH = np.nan
    
    if N == 2 and np.isnan(TH):
        TH = 'opt'
    
    if N >= 3:
        TH = (TH / 100) * iso21920('fnciso21920_feature_parameters_peak_pit',z,len(z) / 5,5).xz
    
    ## parse significant
    significant = significant.replace('%',' %')
    str = significant.split(' ')
    Fsig = str[0]
    N = len(str)
    try:
        NIsig = float(str[np.min([2,N])-1])
    except ValueError:
        NIsig = np.nan
    
    if N >= 3:
        hintersection,material = iso21920('fnciso21920_material_ratio_functions_imp',z,len(z))
        NIsig = np.amax(z) + iso21920('fnciso21920_material_ratio_functions_xcm',material,hintersection,len(material),NIsig,1)
    
    ## parse stats
    str = stats.split(' ')
    Astats = str[0]
    try:
        vstats = float(str[-1])
    except ValueError:
        vstats = np.nan
        
    watershed = Watershed(z,dx,FT,PT,TH)
    M = watershed.motifs()
    xFC,M,ATTR,_,_ = feature_parameter(z,dx,M, Fsig, NIsig, AT, Astats, vstats)
    # meta = struct('ATTR',ATTR,'nM',len(M),'PT',PT,'TH',TH,'Fsig',Fsig,'NIsig',NIsig,'AT',AT,'Astats',Astats,'vstats',vstats)
    # return xFC,M,meta
    return xFC,M