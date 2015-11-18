"""
:mod:`suchyta_utils.balrog` is some functions for working with Balrog catalogs.

"""

import numpy as _np
import esutil as _es
import copy as _copy
import numpy.lib.recfunctions as _rec


def _DefineKwargs(kwargs):
    if 'version' not in kwargs:
        kwargs[version] = 'version'
    if 'id' not in kwargs:
        kwargs['id'] = 'balrog_index'
    if 'out' not in kwargs:
        kwargs['out'] = 'balrog_id'
    kwargs['truth'] = 0
    return kwargs


#def AddUniqueID(truth, *morecats, version='version', id='balrog_index', out='balrog_id'):
def AddUniqueID(*cats, **kwargs):
    """
    Add a column that makes the id which joins the Balrog truth, sim, and nosim tables (called `balrog_index` if you get it from me)
    be unique, even if multiple tables were stacked together such that `id` was no longer unique.
    The function just offsets the ids from different table versions.
    Effectively, it changes a 2D index into a 1D one.

    .. note::
        By "unique", I mean that `out` will be unique within the truth catalog. 
        It is possible for mulitiple detections to be associated with the same truth Balrog objects,
        such that are duplicate `id` entries in `sim` or `nosim` even with a unique identifier in the truth catalog.

    Parameters
    ----------
    cats (1 or more stuctured arrays)
        The leading arguments are non-keyworded, as many catalogs as you want to add the new column to.
        One of these MUST be the truth catalog, whose position in the argument ordering is given by `truth`.

    kwargs (Keyword arguments)
        | version (str) -- A column name that differentiates which table the data came from. Default = ``'version'``
        | id (str) -- Column name that joins truth/sim/nosim. Default = ``'balrog_index'``
        | out (str) -- Name of the new column. Default = ``'balrog_id'``
        | truth (int) -- Position in the arguments given of the truth catalog. Default = ``0``

    Returns
    -------
    cats (1 or more structured arrays)
        A list of the arrays with the new columns, ordered in the same order you gave the function arguments.

    """
    
    kwargs = _DefineKwargs(kwargs)
    version = kwargs['version']
    id = kwargs['id']
    out = kwargs['out']
   
    cats = list(cats)
    for i in range(len(cats)):
        cats[i] = _rec.append_fields(cats[i], out, _np.copy(cats[i][id]))

    versions = _np.unique(cats[kwargs['truth']][version])
    for i in range(len(versions)):
        if i==0:
            continue

        v = versions[i]
        v1 = versions[i-1]

        vcut = (cats[kwargs['truth']][version]==v1)
        amax = _np.amax(cats[kwargs['truth']][vcut][out])
        start = amax + 1

        for j in range(len(cats)):
            cut = (cats[j][version]==v)
            cats[j][out][cut] = cats[j][out][cut] + start

    return cats


def BinnedAvg(cat=None, bins=None, binon=None, avgon=None, kind='avg'):
    """
    Bin the data on some field and then tompute an average of another field in each bin.

    Parameters
    ----------
    cat (Structured array)
        The structured data array (e.g. numpy recarray)
    bins (float array)
        The bin edges to bin the data
    binon (str)
        A column name. The column which the data will be divided into `bins`.
    avgon (str)
        A column name. The column to take an average in each bin.
    kind (str)
        What kind of average to do. Available options are ['avg', 'median']

    Returns
    -------
    avg (float array)
        The averages in each bin

    """
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


def Completeness(sim=None, truth=None, binon='mag_i', bins=None, method='cut', matchon='balrog_index'):
    """
    Compute the completeness of a Balrog catalog.

    Parameters
    ----------
    sim (structured array)
       'sim' data array, i.e. the detected Balrog objects
    truth (stuctured array)
        'truth' data array, i.e. all the truth Balrog objects
    binon (str)
        Column name to use in both the `sim` and `truth` arrays (but it needs to be a truth column).
        What to compute the completeness as a function of.
    bins (array)
        Array of the bin edges to use.

    Returns
    -------
    comp (float array)
        The completeness in each bin
        
    """
    c = _np.zeros(len(bins)-1)

    for i in range(len(bins)-1):
        tcut = (truth[binon] > bins[i]) & (truth[binon] < bins[i+1])
        if method=='cut':
            scut = (sim[binon] > bins[i]) & (sim[binon] < bins[i+1])
        elif method=='match':
            scut = _np.in1d(sim[matchon], truth[tcut][matchon])

        den = float(_np.sum(tcut))
        num = float(_np.sum(scut))

        if den > 0:
            c[i] = num / den
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
    Remove entries from the Balrog measurement catalog ('sim') whose detection position is nearby a previously existing object ('nosim').
    This is done by matching the two catalogs on an index.
        
    Parameters
    ----------
    m (stuctured array)
        Measurement ('sim') data array
    nosim (structured array)
        'Nosim' data array
    id (str)
        Column name to match on
    version (str/None)
        None amounts to assuming that truth values for the `id` column are unique throughout the table.
        However, we often concatenate multiple Balrog tables together, so this is not always the case.
        `version` is a column name which differentiates which table the data came from, as a second indexing column
        to remove the 'nosim' entries from 'sim'

    Returns
    -------
    trimmed (structured array)
        The new data array

    """
    matched = _np.copy(m)
    if version is None:
        matched = _nmatch(matched, nosim, id)
    else:
        matched = _nnmatch(matched, nosim, id, version)
    return matched

def AddModestNeed(keys, release='sva1'):
    """
    Add any missing keys to a list that needed to compute `modest_class` for a dataset.

    Parameters
    ----------
    keys (str array)
        Some set of field names
    release (str)
        Year of the dataset. Allowed values are ['sva1', 'y1a1']

    Returns
    -------
    newkeys (str array)
        A new list which includes any keys which were missing

    """
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

    .. warning::
        Currently, the Balrog setup is not able to completely recreate the modest classification for Y1A1.
        It depends on `wavg_spread_model_i`, a weighted average of the i-band single-epoch `spread_model` measurements.
        Balrog is only running on the coadd images for now. The Y1A1 `modest_class` here ignores the `wavg_spread_model_i` condition.

    Parameters
    ----------
    data (structure array) 
        The data array to find `modest_class` for each object
    release (str)
        Year of the dataset. Allowed values are ['sva1', 'y1a1']


    Returns
    -------
    modest (int array)
        Array of the modest classification for each object

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
