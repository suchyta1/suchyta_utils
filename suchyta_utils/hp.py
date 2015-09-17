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


def RaDec2Healpix(ra=None, dec=None, nside=None, nest=False, cat=None):
    _nsideExcept(nside)
    r, d = _CatOrArrays(cat, ra, dec)

    phi = _np.radians(r)
    theta = _np.radians(90.0 - d)
    hpInd = _hp.ang2pix(nside, theta, phi, nest=nest)
    return hpInd


def ApplyMask(ra=None, dec=None, mask=None, ext=None, nest=False, cat=None, nocut=False):
    if ext is None:
        ext = -1

    if type(mask)==str:
        try:
            ordering = _pyfits.read(mask)[ext]
            if ordering.upper()=='NEST':
                nest = True
            else:
                nest = False
        except:
            pass
        map = _hp.read_map(mask, nest=nest)
    else:
        map = mask

    nside = _hp.npix2nside(map.size)
    _nsideExcept(nside)
    r, d = _CatOrArrays(cat, ra, dec)
    pix = RaDec2Healpix(r, d, nside, nest=nest)
    use = (map[pix]==1)
    if nocut:
        return use
    elif cat is not None:
        return cat[use]
    else:
        return [r[use], d[use]]

