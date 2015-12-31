"""
:mod:`suchyta_utils.balrog` is some functions for working with Balrog catalogs.

"""

import numpy as _np
import esutil as _es
import copy as _copy
import os as _os
import numpy.lib.recfunctions as _rec
from sklearn.neighbors import NearestNeighbors as _NN

import suchyta_utils.hp as _hp
import suchyta_utils.slr as _slr


def _DefineKwargs(kwargs):
    if 'version' not in kwargs:
        kwargs[version] = 'version'
    if 'id' not in kwargs:
        kwargs['id'] = 'balrog_index'
    if 'out' not in kwargs:
        kwargs['out'] = 'balrog_id'
    kwargs['truth'] = 0
    return kwargs


def ApplyCut(cat, key='objtype', val=1, cond='='):
    """
    Apply a cut on a single field to the catalog.

    Parameters
    ----------
    cat (structured array)
        The structured data array (e.g. numpy recarray)
    key (str)
        What field to make the cut on
    cond (str)
        How the cut is to be applied. Possible values: ``=``, ``<``, ``<=``, ``>``, ``>=``
    val (float)
        The value on the right side of `cond`.

    Returns
    -------
    cutcat (structured array)
        The cut version of the catalog

    """

    if cond=='=':
        use = (cat[key]==val)
    elif cond=='<':
        use = (cat[key] < val)
    elif cond=='<=':
        use = (cat[key] <= val)
    elif cond=='>':
        use = (cat[key] > val)
    elif cond=='>=':
        use = (cat[key] >= val)

    return cat[use]


def ColorCutGold(cat, key='mag_auto'):
    """
    Apply Gold "crazy color" cut

    Parameters
    ----------
    cat (structured array)
        The structured data array (e.g. numpy recarray)
    key (str)
        What field to compute the colors using. Strictly speaking DES defines this on `mag_auto` but I kept the field general.

    Returns
    -------
    cutcat (structued array)
        The cut version of the catalog

    """

    gr = cat['%s_g'%(key)] - cat['%s_r'%(key)]
    iz = cat['%s_i'%(key)] - cat['%s_z'%(key)]
    cut1 = (-1 < gr) & (gr < 4)
    cut3 = (-1 < iz) & (iz < 4)
    cut = (cut1 & cut3)
    return cat[cut]


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
        modest[neither] = 0

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



def AngularDistance(ra1, ra2, dec1, dec2):
    """
    Get the angular distance between points (in degrees)

    Parameters
    ----------
    ra1 (float array)
        The ra values of the first points
    ra2 (float array)
        The ra values of the second points
    dec1 (float array)
        The dec values of the first points
    dec2 (float array)
        The dec values of the second points

    Returns
    -------
    dist (float array)
        The angular distance between the points, in degrees

    """

    ddec = _np.radians(dec2-dec1)
    dra = _np.radians(ra2-ra1)
    a = _np.sin(ddec/2) * _np.sin(ddec/2) + _np.cos(_np.radians(dec1)) * _np.cos(_np.radians(dec2)) * _np.sin(dra/2) * _np.sin(dra/2)
    c = 2 * _np.arctan2(_np.sqrt(a), _np.sqrt(1-a))
    return c*180.0/_np.pi


def BadPos(cat):
    """
    Get objects of large astrometric color (windowed offsets between g and i band)

    Parameters
    ----------
    cat (structured array)
        The data array

    Returns
    -------
    cut (bool array)
        An array where the entries that are True are the bad ones
    """

    offset = 3600.0 * AngularDistance(cat['alphawin_j2000_g'], cat['alphawin_j2000_i'], cat['deltawin_j2000_g'], cat['deltawin_j2000_i'])
    c = (cat['fluxerr_auto_g'] > 0)
    s2n = _np.zeros(len(cat), dtype=_np.float32)
    s2n[c] = cat['flux_auto_g'][c] / cat['fluxerr_auto_g'][c]

    bad = _np.zeros(len(cat), dtype=_np.int16)
    cut = (s2n > 5) & (_np.fabs(offset) > 1.0)

    return cut


def UniformRandom(size, ramin=0, ramax=360, decmin=-90, decmax=90):
    """
    Populate uniform randoms on a sphere. The default generates over the whole sphere, and the ra/dec min/max parameters let you choose only a subregion.

    Parameters
    ----------
    size (int)
        Number of random points to generate
    ramin (float)
        Set a minimum RA (degress) for the region in which to generate points
    ramax (float)
        Set a maximum RA (degress) for the region in which to generate points
    decmin (float)
        Set a minimum DEC (degress) for the region in which to generate points
    decmax (float)
        Set a maximum DEC (degress) for the region in which to generate points

    Returns
    -------
    [ra, dec] (float arrays)
        The random positions
    """

    ra = _np.random.uniform(ramin,ramax, size)
    tmin = _np.cos( _np.radians(90.0 - decmax) )
    tmax = _np.cos( _np.radians(90.0 - decmin) )
    theta = _np.degrees( _np.arccos( _np.random.uniform(tmin,tmax, size) ) )
    dec = 90.0 - theta
    return [ra, dec]



def ReweightMatch(keys=None, matchto=None, reweight=None, nn=100):
    """
    Find the weights that match one sample to another, using a nearest neighbor calculation.

    Parameters
    ----------
    keys (str array)
        Sequence of what fields to perform the reweighting
    matchto (structured array)
        The data you are matching to (which will not change)
    reweight (structured array)
        The data you are going to reweight (which will be resampled)
    nn (int)
        Number of nearest neighbors to query

    Returns
    -------
    weight (float array)
        The weights for resampling
    """

    if keys is None:
        raise Exception('have to match something!')

    truthSample_arr = _np.zeros( (len(matchto), len(keys)) )
    matchSample_arr = _np.zeros( (len(reweight), len(keys)) )
    
    for i in range(len(keys)):
        truthSample_arr[:,i] =  matchto[keys[i]]
        matchSample_arr[:,i] =  reweight[keys[i]]

    NP = _calcNN(nn, truthSample_arr, matchSample_arr)
    bad = (NP <= 1)
    wts = NP * 1./nn
    wts[bad] = 0.

    return wts


def _calcNN(Nnei, magP, magT):

    # Find Nnei neighbors around each point
    # in the training sample.
    # Nnei+1 because [0] is the point itself.
    nbrT        = _NN(n_neighbors=Nnei+1, n_jobs=-1).fit(magT)
    distT, indT = nbrT.kneighbors(magT, n_neighbors=Nnei+1)
   
    # Find how many neighbors are there around 
    # each point in the photometric sample 
    # within the radius that contains Nnei in 
    # the training one. 
    nbrP        = _NN(radius=distT[:, Nnei], n_jobs=-1).fit(magP)
    distP, indP = nbrP.radius_neighbors(magT, radius=distT[:, Nnei])


    # Get the number of photometric neighbors
    NP = []
    for i in range(len(distP)):
        NP.append(len(distP[i])-1)

    NP = _np.asarray(NP)
    return NP




class Y1Dataset(object):

    def ApplyFootprint(self, ra='alphawin_j2000_i', dec='deltawin_j2000_i'):
        if self.footprint is not None:
            self.data = _hp.ApplyMask(mask=self.footprint[0], nest=self.footprint[1], cat=self.data, ra=ra, dec=dec, val=1, cond='>=')
        else:
            print 'No footprint to apply'

    def ApplyMask(self, ra='alphawin_j2000_i', dec='deltawin_j2000_i', cond='=', val=0):
        if self.badmask is not None:
            self.data = _hp.ApplyMask(mask=self.badmask[0], nest=self.badmask[1], cat=self.data, ra=ra, dec=dec, val=val, cond=cond)
        else:
            print 'No bad mask to apply'

    def ApplyDepth(self, ra='alphawin_j2000_i', dec='deltawin_j2000_i', cond='>', val=22):
        if self.depth is not None:
            self.data = _hp.ApplyMask(mask=self.depth[0], nest=self.depth[1], cat=self.data, ra=ra, dec=dec, val=val, cond=cond)
        else:
            print 'No depth mask to apply'


    def UsualMasking(self, ra='alphawin_j2000_i', dec='deltawin_j2000_i', cond='=', val=0):
        self.ApplyFootprint(ra=ra, dec=dec)
        self.ApplyMask(ra=ra, dec=dec, cond=cond, val=val)

    def BenchmarkMasking(self, ra='alphawin_j2000_i', dec='deltawin_j2000_i'):
        self.UsualMasking(ra=ra, dec=dec)
        self.ApplyDepth(cond='>', val=22)


    def ApplySLR(self, ra='alphawin_j2000_i', dec='deltawin_j2000_i', bands=['g','r','i','z','y'], key='mag_auto'):
        if self.slr is not None:
            for band in bands:
                shift = self.slr.GetMagShifts(band, self.data[ra], self.data[dec])
                self.data['{1}_{0}'.format(band,key)] = self.data['{1}_{0}'.format(band,key)] + shift
        else:
            print 'No SLR to apply'

    def ApplyModest(self, val):
        modest = Modest(self.data, release='y1a1')
        self.data = self.data[ modest==val ]


    def ModestGalaxies(self):
        self.ApplySLR(key='mag_auto')
        self.ApplyModest(1)

    def BenchmarkQualityCuts(self):
        # Crazy Colors
        gr = self.data['mag_auto_g'] - self.data['mag_auto_r']
        ri = self.data['mag_auto_r'] - self.data['mag_auto_i']
        iz = self.data['mag_auto_i'] - self.data['mag_auto_z']
        self.data = self.data[ (gr > -1) & (gr < 3) & (ri > -1) & (ri < 2.5) & (iz > -1) & (iz < 2) ]

        # Cuts to do the flags_gold that I can. (BM color cuts supercede the gold ones. Can't really do BBJ.)
        for band in ['g','r','i','z']:
            cut = (self.data['flags_%s'%(band)] > 3)
            self.data = self.data[-cut]
        bad = BadPos(self.data)
        self.data = self.data[-bad]


    def BenchmarkGalaxies(self, ra='alphawin_j2000_i', dec='deltawin_j2000_i'):
        self.BenchmarkMasking(ra=ra, dec=dec)
        self.ModestGalaxies()
        self.BenchmarkQualityCuts()

        # 'Complete'
        self.data = self.data[ (self.data['mag_auto_i'] > 17.5) & (self.data['mag_auto_i'] < 22) ]


    def PseudoBenchmarkGalaxies(self, ra='alphawin_j2000_i', dec='deltawin_j2000_i'):
        self.UsualMasking(ra=ra, dec=dec)
        self.ModestGalaxies()
        self.BenchmarkQualityCuts()


    def __init__(self, data, processing=None):
        self.data = data

        if processing is not None:
            done = self.__dict__.keys()
            for key in processing.__dict__.keys():
                if key not in done:
                    self.__setattr__(key, processing.__dict__[key])
        else:
            print 'WARNING: did not give an object for masking, SLR, etc.'


class Y1Processing(object):

    def __init__(self, dir=None, version='1.0.2', subset='wide'):
        if dir is None:
            print 'WARNING: must specify a directory where the the mask(s), SLR, etc. live. Nothing done.'

        else:
            self.version = version
            self.dir = dir
            self.subset = subset
            self.footprint_file = _os.path.join(self.dir, 'y1a1_gold_%s_%s_footprint_4096.fits.gz' %(self.version,self.subset))
            self.badmask_file = _os.path.join(self.dir, 'y1a1_gold_%s_%s_badmask_4096.fits.gz'%(self.version,self.subset))
            self.depth_file = _os.path.join(self.dir, 'y1a1_gold_%s_%s_auto_nside4096_i_10sigma.fits.gz'%(self.version,self.subset))
            self.slrcode = _os.path.join(self.dir, 'y1a1_slr_shiftmap.py')
            self.slrfits = _os.path.join(self.dir, 'y1a1_%s_slr_wavg_zpshift2.fit'%(self.subset))
            
            if not _os.path.exists(self.footprint_file):
                print 'WARNING: the footprint file does not exist for your given directory, subset, and version: %s'%(self.footprint_file)
                self.footprint = None
            else:
                self.footprint = _hp.GetHPMap(self.footprint_file)

            if not _os.path.exists(self.badmask_file):
                print 'WARNING: the bad mask file does not exist for your given directory, subset, and version: %s'%(self.badmask_file)
                self.badmask = None
            else:
                self.badmask = _hp.GetHPMap(self.badmask_file)

            if not _os.path.exists(self.depth_file):
                print 'WARNING: the depth mask file does not exist for your given directory, subset, and version: %s'%(self.depth_file)
                self.depth = None
            else:
                self.depth = _hp.GetHPMap(self.depth_file)

            noslr = False
            if not _os.path.exists(self.slrcode):
                print 'WARNING: SLR python file does not exist for your given directory and subset: %s'%(self.slrcode)
                noslr = True
            if not _os.path.exists(self.slrfits):
                print 'WARNING: SLR FITS file does not exist for your given directory and subset: %s'%(self.slrfits)
                noslr = True
           
            if not noslr:
                self.slr = _slr.SLR(release='y1a1', area=self.subset, slrdir=self.dir)
            else:
                self.slr = None

