import numpy as _np
import esutil as _es


def Completeness(sim=None, truth=None, bins=None):
    c = _np.zeros(len(bins)-1)
    for i in range(len(bins)-1):
        scut = (sim > bins[i]) & (sim < bins[i+1])
        tcut = (truth > bins[i]) & (truth < bins[i+1])
        print bins[i], bins[i+1], _np.sum(scut), _np.sum(tcut)
        c[i] = float(_np.sum(scut)) / _np.sum(tcut)
    return c
