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
    hps = RaDec2Healpix(cat=cat, ra=ra, dec=dec, nside=nside)
    uhps = _np.unique(hps)
    area = len(uhps) * _hp.nside2pixarea(nside, degrees=True)
    return area

def Healpix2RaDec(pixels, nside=None, nest=False):
    _nsideExcept(nside)
    theta, phi = _hp.pix2ang(nside, pixels, nest=nest)
    dec = 90.0 - _np.degrees(theta)
    ra = _np.degrees(phi)
    return [ra, dec]


def RaDec2Healpix(ra=None, dec=None, nside=None, nest=False, cat=None):
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


def ApplyMask(ra=None, dec=None, mask=None, ext=None, nest=False, cat=None, nocut=False):
    if ext is None:
        ext = -1

    if type(mask)==str:
        nest = _NestFromHeaderHP(mask, ext)
        map = _hp.read_map(mask, nest=nest)
    else:
        map = mask

    nside = _hp.npix2nside(map.size)
    pix = RaDec2Healpix(cat=cat, ra=ra, dec=dec, nside=nside, nest=nest)
    use = (map[pix]==1)
    if nocut:
        return use
    elif cat is not None:
        return cat[use]
    else:
        return [ra[use], dec[use]]

def GetHPMap(mask):
    nest = _NestFromHeaderHP(mask, -1)
    map = _hp.read_map(mask, nest=nest)
    return map, nest


def GetBorisMap(file, ext=-1, nside=None, nest=None):
    nest, nside = _BFromHeader(file, ext, nside=nside, nest=nest)
    data = _pyfits.open(file)[ext].data
    pix = data['PIXEL']
    value = data['SIGNAL']
    map = _np.zeros(_hp.nside2npix(nside))
    map[:] = _hp.UNSEEN
    map[pix] = value
    return map, nest

def RaDec2MapValue(map=None, nest=False, cat=None, ra=None, dec=None):
    nside = _hp.npix2nside(map.size)
    pix = RaDec2Healpix(ra=ra, dec=dec, nside=nside, nest=nest, cat=cat)
    return map[pix]
