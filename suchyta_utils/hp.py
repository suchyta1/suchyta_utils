"""
The :mod:`hp` module is for working with HEALPix maps in python.
Most of the functions are things I commonly do which make healpy calls underneath.
    
Examples:
::
    
    # Apply the Y1A1 footprint and 4% exclusion mask to Y1 DES and Balrog data
    sim_y1 = esutil.io.read('sim-y1.fits')
    des_y1 = esutil.io.read('des-y1.fits')
    y1fmask, y1fnest = es.hp.GetHPMap('sva1_gold_1.0.4_goodregions_04_equ_ring_4096.fits.gz')
    y1mask, y1nest = es.hp.GetHPMap('y1a1_gold_1.0.1_wide+d04_badmask_4096.fit.gz')
    sim_y1 = es.hp.ApplyMask(mask=y1fmask, nest=y1fnest, cat=sim_y1, ra='alphawin_j2000_i', dec='deltawin_j2000_i', val=1, cond='=')
    des_y1 = es.hp.ApplyMask(mask=y1fmask, nest=y1fnest, cat=des_y1, ra='alphawin_j2000_i', dec='deltawin_j2000_i', val=1, cond='=')
    sim_y1 = es.hp.ApplyMask(mask=y1mask, nest=y1nest, cat=sim_y1, ra='alphawin_j2000_i', dec='deltawin_j2000_i', val=0, cond='=')
    des_y1 = es.hp.ApplyMask(mask=y1mask, nest=y1nest, cat=des_y1, ra='alphawin_j2000_i', dec='deltawin_j2000_i', val=0, cond='=')

    # Apply the SVA1 4% exclusion mask to SVA1 DES and Balrog data
    sim_sv = esutil.io.read('sim-sv.fits')
    des_sv = esutil.io.read('des_sv.fits')
    svmask, svnest = es.hp.GetHPMap('sva1_gold_1.0.4_goodregions_04_equ_ring_4096.fits.gz')
    sim_sv = es.hp.ApplyMask(mask=svmask, nest=svnest, cat=sim_sv, ra='alphawin_j2000_i', dec='deltawin_j2000_i', val=1, cond='=')
    des_sv = es.hp.ApplyMask(mask=svmask, nest=svnest, cat=des_sv, ra='alphawin_j2000_i', dec='deltawin_j2000_i', val=1, cond='=')

"""

import numpy as _np
import healpy as _hp
import astropy.io.fits as _pyfits


class _HPException(Exception):
    pass

def _nsideExcept(nside):
    if nside is None:
        raise _HPException('must specify nside')


def _CatOrArrays(cat, ra, dec):
    if cat is not None:
        r = cat[ra]
        d = cat[dec]
    else:
        r = ra
        d = dec
    return r, d


def GetArea(cat=None, ra=None, dec=None, nside=4096, nest=False):
    """
    Return the area covered by the dataset, in square degrees.
    The function checks if any objects are each HEALPixel, then adds up how many pixels were found, and multiplies by the pixel area.
   
    * ``cat``: if not None, equivalent to a numpy recarray
    * ``ra``: if `cat` is None, an array of the RA values. Otherwise, the column name for the RA column in `cat`
    * ``dec``: if `cat` is None, an array of the DEC values. Otherwise, the column name for the DEC column in `cat`
    * ``nside``: HEALPix nside; determines the pixel size used in the computation
    * ``nest``: whether or not to use nested format, though this parameter is irrelevant (I just allow it for completeness)

    """
    hps = RaDec2Healpix(cat=cat, ra=ra, dec=dec, nside=nside, nest=nest)
    uhps = _np.unique(hps)
    area = len(uhps) * _hp.nside2pixarea(nside, degrees=True)
    return area

def Healpix2RaDec(pixels, nside=None, nest=False):
    """
    Return the RA/DEC for each HEALPixel index in `pixels`. A 2D list is returned. 
    Element 0 is the RA array, and Element 1 is the DEC array.
   
    * ``pixels``: Array of the pixel integers 
    * ``nside``: HEALPix nside
    * ``nest``: whether or not to use nested format

    """
    _nsideExcept(nside)
    theta, phi = _hp.pix2ang(nside, pixels, nest=nest)
    dec = 90.0 - _np.degrees(theta)
    ra = _np.degrees(phi)
    return [ra, dec]


def RaDec2Healpix(ra=None, dec=None, nside=None, nest=False, cat=None):
    """
    Return the HEALPix index for each RA/DEC. Numpy array is returned. 
   
    * ``cat``: if not None, equivalent to a numpy recarray
    * ``ra``: if `cat` is None, an array of the RA values. Otherwise, the column name for the RA column in `cat`
    * ``dec``: if `cat` is None, an array of the DEC values. Otherwise, the column name for the DEC column in `cat`
    * ``nside``: HEALPix nside
    * ``nest``: whether or not to use nested format

    """
    _nsideExcept(nside)
    r, d = _CatOrArrays(cat, ra, dec)

    phi = _np.radians(r)
    theta = _np.radians(90.0 - d)
    hpInd = _hp.ang2pix(nside, theta, phi, nest=nest)
    return hpInd


def _NestFromHeaderHP(mask, ext, nside=False):
    nest = False
    header = _pyfits.open(mask)[ext].header
    try:
        ordering = header['ORDERING']
        if ordering.upper()=='NEST':
            nest = True
    except:
        pass

    if nside:
        nside = header['NSIDE']
        return nest, nside
    else:
        return nest

def _BFromHeader(file, ext=-1, nest=None, nside=None):
    if nest is None:
        if nside is None:
            nest, nside = _NestFromHeaderHP(file, ext, nside=True)
        else:
            nest = _NestFromHeaderHP(file, ext, nside=False)
    return nest, nside


def ApplyMask(ra=None, dec=None, mask=None, ext=None, nest=False, cat=None, nocut=False, val=1, cond='='):
    """
    Apply a HEALPix mask. 

    * ``mask``: filename or an array (such as the first element returned by :func:`~suchyta_utils.hp.GetHPMap`).
      If you're going to be applying the same mask repeated times, it's faster to use the array so you don't have to read the file multiple times.
    * ``cat``: if not None, equivalent to a numpy recarray
    * ``ext``: Only relevant for `mask` which is a filename. Extension in the file where the map lives. None is equivalent to -1, the final extension.
    * ``ra``: if `cat` is None, an array of the RA values. Otherwise, the column name for the RA column in `cat`
    * ``dec``: if `cat` is None, an array of the DEC values. Otherwise, the column name for the DEC column in `cat`
    * ``nest``: whether or not to use nested format. If `mask` is a file name, the `nest` keyword is irrelevant (because it is automatically determined from the header).
    * ``nocut``: If True, do not apply the map to the data, and instead return the array of Booleans if the entry should be kept or not.
    * ``cond``: How the cut is to be applied. Possible values: ``=``, ``<``, ``<=``, ``>``, ``>=``
    * ``val``: The value on the right side of `cond`.


    """
    if ext is None:
        ext = -1

    if type(mask)==str:
        nest = _NestFromHeaderHP(mask, ext)
        map = _hp.read_map(mask, nest=nest)
    else:
        map = mask

    nside = _hp.npix2nside(map.size)
    pix = RaDec2Healpix(cat=cat, ra=ra, dec=dec, nside=nside, nest=nest)

    if cond=='=':
        use = (map[pix]==val)
    elif cond=='<':
        use = (map[pix] < val)
    elif cond=='<=':
        use = (map[pix] <= val)
    elif cond=='>':
        use = (map[pix] > val)
    elif cond=='>=':
        use = (map[pix] >= val)
    
    if nocut:
        return use
    elif cat is not None:
        return cat[use]
    else:
        return [ra[use], dec[use]]

def GetHPMap(mask):
    """
    Read and return a HEALPix mask file in FITS (.fz) format,
    returning a 2D tuple: (map [i.e. an array], nest)
    This function calls ``healpy.read_map(mask, nest=nest)``, automatically populating `nest` by reading the image header.

    * ``mask``: filename

    """
    nest = _NestFromHeaderHP(mask, -1)
    map = _hp.read_map(mask, nest=nest)
    return map, nest


def GetBorisMap(file, ext=-1, nside=None, nest=None):
    """
    Read one of Boris Leistedt's maps and return it as a HEALPix map, returning a 2D tuple: (map [i.e. an array], nest)

    * ``file``: Filename of the map
    * ``ext``: Extension of the file which the map lives in.
    * ``nside``: None means read from header. I would only use None.
    * ``nest``: None means read from header. I would only use None.

    """
    nest, nside = _BFromHeader(file, ext, nside=nside, nest=nest)
    data = _pyfits.open(file)[ext].data
    pix = data['PIXEL']
    value = data['SIGNAL']
    map = _np.zeros(_hp.nside2npix(nside))
    map[:] = _hp.UNSEEN
    map[pix] = value
    return map, nest

def RaDec2MapValue(map=None, nest=False, cat=None, ra=None, dec=None):
    """
    Return the value of a map for each entry in the catalog. Numpy array is returned. 
   
    * ``cat``: if not None, equivalent to a numpy recarray
    * ``ra``: if `cat` is None, an array of the RA values. Otherwise, the column name for the RA column in `cat`
    * ``dec``: if `cat` is None, an array of the DEC values. Otherwise, the column name for the DEC column in `cat`
    * ``nest``: whether or not the map is in nested format

    """
    nside = _hp.npix2nside(map.size)
    pix = RaDec2Healpix(ra=ra, dec=dec, nside=nside, nest=nest, cat=cat)
    return map[pix]
