import numpy as np


def zero_crossing(z, dx):
    PF = []
    n = len(z)
    i = 0
    j = 1
    nPF = 0
    
    while j < n:
        if (z[j - 1] <= 0 and z[j] > 0) or (z[j - 1] >= 0 and z[j] < 0):
            nPF += 1
            kh = np.argmax(np.abs(z[i:j])) + i
            PF.append({'t': np.sign(z[kh]), 'h': np.abs(z[kh]), 'xh': kh * dx})
            i = j
        j += 1
    
    nPF += 1
    kh = np.argmax(np.abs(z[i:n])) + i
    PF.append({'t': np.sign(z[kh]), 'h': np.abs(z[kh]), 'xh': kh * dx})
    
    return PF, nPF

def maximum_height(zinput, dx, n_sc=5):
    try: 
        # formatting
        z = zinput.copy()
        if z.ndim == 2:
            z = z.reshape(-1,)

        PF, nPF = zero_crossing(z.reshape(-1,), dx)
        
        # number of values per section
        l_sc = int(np.floor(len(z) / n_sc) * dx)

        # Initialize lists for Rpi and Rvi
        Rpi = np.zeros(n_sc)
        Rvi = np.zeros(n_sc)

        # computing Rpi and Rvi
        for i in range(n_sc):
            Rpsc = []
            Rvsc = []
            for j in range(nPF):
                if PF[j]['xh'] >= (i * l_sc) and PF[j]['xh'] < (i + 1) * l_sc:
                    if PF[j]['t'] == 1:
                        Rpsc.append(PF[j]['h'])
                    else:
                        Rvsc.append(PF[j]['h'])
            
            Rpi[i] = max(Rpsc, default=0)
            Rvi[i] = max(Rvsc, default=0)
        
        Rzi = Rpi + Rvi
        rz = np.mean(Rzi)

    except Exception:
        raise RuntimeError("Rz computation failed.")    

    return rz