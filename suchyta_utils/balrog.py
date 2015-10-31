import numpy as _np
import esutil as _es


def Completeness(sim=None, truth=None, bins=None):
    c = _np.zeros(len(bins)-1)
    for i in range(len(bins)-1):
        scut = (sim > bins[i]) & (sim < bins[i+1])
        tcut = (truth > bins[i]) & (truth < bins[i+1])
        den = _np.sum(tcut)

        if den > 0:
            c[i] = float(_np.sum(scut)) / den
        else:
            c[i] = _np.nan

    return c


def _nmatch(matched, nosim, id):
    kill = _np.in1d(matched[id], nosim[id])
    matched = matched[-kill]
    return matched

def _nnmatch(matched, nosim, id, version):
    versions = _np.unique(nosim[version])
    for v in versions:
        mv_cut = (matched[version]==v)
        nv_cut = (nosim[version]==v)
        bad = nosim[nv_cut][id]
        mb_cut = _np.in1d(matched[id],bad)
        both = (mv_cut & mb_cut)
        matched = matched[-both]
    return matched


def RemoveNosim(m, nosim, version=None, id='balrog_index'):
    matched = _np.copy(m)
    if version is None:
        matched = _nmatch(matched, nosim, id)
    else:
        matched = _nnmatch(matched, nosim, id, version)
    return matched
