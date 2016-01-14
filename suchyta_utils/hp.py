"""
The :mod:`suchyta_utils.hp` submodule is for working with HEALPix maps in python.
Most of the functions are things I commonly do which make healpy calls underneath.
The functions are generally designed to be able to intuitively handle either objects which amount to recarrays
or seperate RA/DEC arrays.

Here's a snippet from the git repository's `example.py <https://github.com/suchyta1/suchyta_utils/blob/master/examples/example.py>`_,
which demonstrates how to use some of the functions::

    # Put all the supplementary DES files in some directory
    supdir = './supplementary'

    # Y1 files
    y1_data_file = os.path.join(supdir, 'balrog-small-y1.fits')
    y1_footprint_file = os.path.join(supdir, 'y1a1_gold_1.0.1_wide+d04_footprint_4096.fit.gz')
    y1_badmask_file = os.path.join(supdir, 'y1a1_gold_1.0.1_wide+d04_badmask_4096.fit.gz')

    # Apply the Y1A1 footprint and 4% exclusion mask to the Balrog data. 
    # Do everything with objects, not file names. This is always at least as fast as all filen names. (Particularly if one uses things multiple times.)
    y1_data = esutil.io.read(y1_data_file)
    y1_footprint, y1_footprint_nest = es.hp.GetHPMap(y1_footprint_file)
    y1_data = es.hp.ApplyMask(mask=y1_footprint, nest=y1_footprint_nest, cat=y1_data, ra='alphawin_j2000_i', dec='deltawin_j2000_i', val=1, cond='>=')
    y1_badmask, y1_badmask_nest = es.hp.GetHPMap(y1_badmask_file)
    y1_data = es.hp.ApplyMask(mask=y1_badmask, nest=y1_badmask_nest, cat=y1_data, ra='alphawin_j2000_i', dec='deltawin_j2000_i', val=0, cond='=')

    # SV files
    sv_data_file = os.path.join(supdir, 'balrog-small-sv.fits')
    sv_goodmask_file = os.path.join(supdir, 'sva1_gold_1.0.4_goodregions_04_equ_ring_4096.fits.gz')

    # Apply the SVA1 4% exclusion mask to the Balrog data. Do everything with file names.
    sv_data = esutil.io.read(sv_data_file)
    sv_ra, sv_dec = es.hp.ApplyMask(mask=sv_goodmask_file, ra=sv_data['alphawin_j2000_i'], dec=sv_data['deltawin_j2000_i'], val=1, cond='=')
    
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


def GetArea(cat=None, ra=None, dec=None, nside=4096):
    """
    Compute the area covered by the dataset, in square degrees.
    The function checks if any objects are each in HEALPixel, then adds up how many pixels were found, and multiplies by the pixel area.

    Parameters
    ----------
    cat (None/structured array)
        If not None, the structured data array (e.g. numpy recarray)
    ra (float array/str)
        If `cat` is None, an array of the RA values. Otherwise, the column name for the RA column in `cat`.
    dec (float array/str)
        if `cat` is None, an array of the DEC values. Otherwise, the column name for the DEC column in `cat`.
    nside (int)
        HEALPix nside, which determines the pixel size used in the computation

    Returns
    -------
    area (float)
        The area covered by the dataset (in square degrees).

    """
    hps = RaDec2Healpix(cat=cat, ra=ra, dec=dec, nside=nside, nest=False)
    uhps = _np.unique(hps)
    area = len(uhps) * _hp.nside2pixarea(nside, degrees=True)
    return area

def Healpix2RaDec(pixels, nside=None, nest=False):
    """
    Compute the RA/DEC for each HEALPixel index in `pixels`.
  
    Parameters
    ----------
    pixels (int array)
        Array of the pixel integers 
    nside (int)
        HEALPix nside
    nest (bool)
        Whether or not to use nested format

    Returns
    -------
    ra (float array)
        The RA coordinates of the entries (in degrees)
    dec (float array)
        The DEC coordinates of the entries (in degrees)

    """
    _nsideExcept(nside)
    theta, phi = _hp.pix2ang(nside, pixels, nest=nest)
    dec = 90.0 - _np.degrees(theta)
    ra = _np.degrees(phi)
    return [ra, dec]


def RaDec2Healpix(ra=None, dec=None, nside=None, nest=False, cat=None):
    """
    Compute the HEALPix index for each RA/DEC pair. Numpy array is returned. 
  
    Parameters
    ----------
    cat (None/structured array)
        If not None, the structured data array (e.g. numpy recarray)
    ra (float array/str) 
        if `cat` is None, an array of the RA values. Otherwise, the column name for the RA column in `cat`.
    dec (float array/str)
        if `cat` is None, an array of the DEC values. Otherwise, the column name for the DEC column in `cat`.
    nside (int)
        HEALPix nside
    nest (bool)
        Whether or not to use nested format

    Returns
    -------
    index (int)
        Array of the HEALPix indexes for each RA/DEC pair.

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
        if ordering.upper()=='NESTED':
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


def ApplyMask(ra=None, dec=None, mask=None, nest=False, cat=None, nocut=False, val=1, cond='='):
    """
    Apply a HEALPix mask. 

    Parameters
    ----------
    mask (array/str) 
        Either a filename or the mask as an array (such as the first element returned by :func:`~suchyta_utils.hp.GetHPMap`).
        If you're going to be applying the same mask repeated times, it's faster to use the array so you don't have to read the file multiple times.
    cat (None/stuctured array)  
        If not None, the structured data array (e.g. numpy recarray)
    ra (float array/str)
        If `cat` is None, an array of the RA values. Otherwise, the column name for the RA column in `cat`.
    dec (float array/str)
        If `cat` is None, an array of the DEC values. Otherwise, the column name for the DEC column in `cat`.
    nest (bool)
        whether or not to use nested format. If `mask` is a file name, the `nest` keyword is irrelevant (because it is automatically determined from the header).
    nocut (bool)
        If True, do not apply the map to the data, and instead return the array of booleans if the entry should be kept or not.
    cond (str)
        How the cut is to be applied. Possible values: ``=``, ``<``, ``<=``, ``>``, ``>=``
    val (float)
        The value on the right side of `cond`.

    Returns
    -------
    rdata (bool array/structured array/[RA, DEC])
        | Three possibilities exist for how the output is returned:
        | 1) ``if nocut``, return a numpy array of booleans for which entries should be kept; i.e. the new catalog would be ``cat = cat[rdata]``
        | 2) ``if (cat is not None) and (not nocut)``, the portion of `cat` which survives the masking is returned
        | 3) ``if (cat is None) and (not nocut)``, return a 2D list of the [RA, DEC] that survive the masking.

    


    """

    if type(mask)==str:
        nest = _NestFromHeaderHP(mask, -1)
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

def GetHPMap(file):
    """
    Read a HEALPix map file in FITS (.fz) format,
    This function calls ``healpy.read_map(file, nest=nest)``, automatically populating `nest` by reading the image header.

    Parameters
    ----------
    file (str)
        Filename of the map

    Returns
    -------
    map (array)
        The HEALPix map. This amounts to an array, ordered by the HEALPix index.
    nest (bool)
        Whether or not the map is in nested format

    """
    nest = _NestFromHeaderHP(file, -1)
    map = _hp.read_map(file, nest=nest)
    return map, nest


def GetBorisMap(file):
    """
    Read one of Boris Leistedt's maps and convert it to a HEALPix map.

    Parameters
    ----------
    file (str)
        Filename of the map

    Returns
    -------
    map (array)
        The HEALPix map. This amounts to an array, ordered by the HEALPix index.
    nest (bool)
        Whether or not the map is in nested format

    """
    ext = -1
    nside = None
    nest = None

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
    Find the value of a map for each entry in the catalog.
  
    Parameters
    ----------
    map (structured array)
        The HEALPix map 
    cat (None/structured array)
        If not None, the structured data array (e.g. numpy recarray)
    ra (float array/str)
        If `cat` is None, an array of the RA values. Otherwise, the column name for the RA column in `cat`.
    dec (float array/str)
        if `cat` is None, an array of the DEC values. Otherwise, the column name for the DEC column in `cat`.
    nest (bool)
        Whether or not to use nested format.

    Returns
    -------
    vals (array)
        Array of the map values for each entry

    """

    if type(map)==str:
        nest = _NestFromHeaderHP(map, -1)
        mask = _hp.read_map(map, nest=nest)
    else:
        mask = map

    nside = _hp.npix2nside(mask.size)
    pix = RaDec2Healpix(ra=ra, dec=dec, nside=nside, nest=nest, cat=cat)
    return mask[pix]
