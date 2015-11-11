import numpy as _np
import esutil as _es
import copy as _copy


def BinnedAvg(cat=None, bins=None, binon=None, avgon=None, kind='avg'):
    a = _np.zeros(len(bins)-1)
    for i in range(len(bins)-1):
        cut =  (cat[binon] > bins[i]) & (cat[binon] < bins[i+1])
        if _np.sum(cut) > 0:
            if kind=='avg':
                avg = _np.average(cat[cut][avgon])
            elif kind=='median':
                avg = _np.median(cat[cut][avgon])
            a[i] = avg
        else:
            a[i] = _np.nan
    return a


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
    """
    Do some stuff

    Parameters
    ----------
    m: catalog
        Something

    Returns
    -------
    thing: some type
        Explain

    """
    matched = _np.copy(m)
    if version is None:
        matched = _nmatch(matched, nosim, id)
    else:
        matched = _nnmatch(matched, nosim, id, version)
    return matched

def AddModestNeed(keys, release='sva1'):
    k = _copy.copy(keys)
    if release=='sva1':
        extra = ['flags_i', 'class_star_i', 'mag_auto_i', 'spread_model_i', 'spreaderr_model_i', 'mag_psf_i']
    elif release=='y1a1':
        extra = ['spread_model_i', 'spreaderr_model_i']

    for e in extra:
        if e not in keys:
            k.append(e)

    return k


def Modest(data, release='sva1'):
    """
    Find the modest class classifications for the data. The array must have all the required fields or you'll get an error.

    * ``data``: Equivalent to a numpy recarray
    * ``release``: Allowed values ['sva1', 'y1a1']

    .. warning::
        Currently, the Balrog setup is not able to completely recreate the modest classification for Y1A1.
        It depends on `wavg_spread_model_i`, a weighted average of the i-band single-epoch `spread_model` measurements.
        Balrog is only running on the coadd images for now.

    """

    modest = _np.zeros(len(data), dtype=_np.int32)

    if release=='sva1':
        galcut = (data['flags_i'] <=3) & -( ((data['class_star_i'] > 0.3) & (data['mag_auto_i'] < 18.0)) | ((data['spread_model_i'] + 3*data['spreaderr_model_i']) < 0.003) | ((data['mag_psf_i'] > 30.0) & (data['mag_auto_i'] < 21.0)))
        modest[galcut] = 1

        starcut = (data['flags_i'] <=3) & ((data['class_star_i'] > 0.3) & (data['mag_auto_i'] < 18.0) & (data['mag_psf_i'] < 30.0) | (((data['spread_model_i'] + 3*data['spreaderr_model_i']) < 0.003) & ((data['spread_model_i'] +3*data['spreaderr_model_i']) > -0.003)))
        modest[starcut] = 2

        neither = -(galcut | starcut)
        modest[neither] = 3

    elif release=='y1a1':
        ncut = (data['spread_model_i'] + (5.0/3.0)*data['spreaderr_model_i'] < -0.002)
        modest[ncut] = 0

        starcut = (_np.fabs(data['spread_model_i'] + (5.0/3.0)*data['spreaderr_model_i']) < 0.002)
        modest[starcut] = 2

        galcut = (data['spread_model_i'] + (5.0/3.0)*data['spreaderr_model_i'] > 0.005) #& (-(abs(data['wavg_spread_model_i']) < 0.002 and MAG_AUTO_I < 21.5) then 1)
        modest[galcut] = 1

        nncut = -(galcut | starcut | ncut)
        modest[nncut] = 3

    return modest
