import matplotlib as _mpl
import matplotlib.pyplot as _plt
import os as _os
import suchyta_utils as _es
import healpy as _hp
import numpy as _np

from mpl_toolkits.basemap import Basemap as _Basemap
import matplotlib.cm as _cm


def Setup():
    dir = _os.path.dirname(_os.path.realpath(__file__))
    style = _ReadStyle(_os.path.join(dir,'custom-sytle.mpl'))
    style = _SetStyle(style)


def _ReadStyle(file):                                                                                     
    d = {}                                                                                               
    with open(file) as f:                                                                                
        for line in f:                                                                                   
            l = line.strip()                                                                             
            if (l=='') or (l[0]== '#') :                                                                 
                continue                                                                                 
                                                                                                         
            key, val = l.split(':')                                                                      
            k = key.strip()                                                                              
            v = val.strip()                                                                              
                                                                                                         
            if v[0] in ["'", '"']:                                                                       
                d[k] = v                                                                                 
            else:                                                                                        
                try:                                                                                     
                    d[k] = float(v)                                                                      
                except:                                                                                  
                    d[k] = v                                                                             
    return d                                                                                             

                                                                                                         
def _SetStyle(style):                                                                                     
    for k,v in style.iteritems():                                                                        
        _mpl.rcParams[k]=v 


def NTicks(ax, nxticks=None, nyticks=None):
    if nyticks is not None:
        ax.yaxis.set_major_locator( _mpl.ticker.MaxNLocator(nbins=(nyticks-1)) )
        ax.yaxis.set_major_locator( _mpl.ticker.MaxNLocator(nbins=(nyticks-1)) )
    if nxticks is not None:
        ax.xaxis.set_major_locator( _mpl.ticker.MaxNLocator(nbins=(nxticks-1)) )
        ax.xaxis.set_major_locator( _mpl.ticker.MaxNLocator(nbins=(nxticks-1)) )


def OffsetX(r, offset=0, log=False):
    if log:
        newr = _np.log10(r)
    else:
        newr = copy.copy(r)
    
    newr = newr + offset
    if log:
        newr = _np.power(10, newr)

    return newr


def LineSegment(ax=None, left=None, right=None, plotkwargs={}):
    if ax is None:
        fig, ax = _plt.subplots(1,1)

    ax.plot( [left[0], right[0]], [left[1], right[1]], **plotkwargs )
    return ax


def _getMapLocation(cat=None, ra=None, dec=None, nside=512, nest=False, map=None):
    ipix = _es.hp.RaDec2Healpix(cat=cat, ra=ra, dec=dec, nside=nside, nest=nest)
    upix = _np.unique(ipix)
    vals = map[upix]
    lon, lat = _es.hp.Healpix2RaDec(upix, nside=nside, nest=nest)
    return vals, lat, lon

def _getCountLocation(cat=None, ra=None, dec=None, nside=512, nest=False):
    ipix = _es.hp.RaDec2Healpix(cat=cat, ra=ra, dec=dec, nside=nside, nest=nest)
    bc = _np.bincount(ipix)
    pixels = _np.nonzero(bc)[0]
    bc = bc[bc>0] / _hp.nside2resol(nside, arcmin=True)**2 # in arcmin^-2
    lon, lat = _es.hp.Healpix2RaDec(pixels, nside=nside, nest=nest)
    return bc, lat, lon


def _Lon2RA(lon):
    h = lon / 360.0 * 24.0
    hh = int(h)
    m = _np.round( (h-hh)*60 )
    return "%d:%02dh" % (h,m)

def _lon2RA(lon):
    # reverse direction
    lon = 360 - lon
    hours = int(lon)/15
    minutes = int(float(lon - hours*15)/15 * 60)
    minutes = '{:>02}'.format(minutes)
    return "%d:%sh" % (hours, minutes)


def _BasePlot(ax=None, fig=None, nside=512, cat=None, ra=None, dec=None, nest=False, parallels=_np.arange(-75.,0.,5.), meridians=_np.arange(0.,360.,5.), dims=[20,20], center=[-75,-52.5], vmin=None, vmax=None, clabel='$n_s\ [\mathrm{arcmin}^{-2}]$', rafmt='d', raflip=True, xoffset=0, f=_getCountLocation, extrakwargs={} ):
    if fig is None:
        fig, ax = plt.subplots(1,1, figsize=(6.5*nside/512,6*nside/512))
    if ax is None:
        fig, ax = _plt.subplots(1,1)

    r = 6.371e6
    h = r * _np.radians(dims[1])
    w = r * _np.cos(_np.radians(center[1])) * _np.radians(dims[0])

    m = _Basemap(projection='aea', width=w, height=h, lat_0=center[1], lon_0=center[0], ax=ax)
    bc, lat, lon = f(cat=cat, ra=ra, dec=dec, nside=nside, nest=nest, **extrakwargs)
    x,y  = m(lon, lat)

    if vmin is None:
        vmin = _np.amin(bc)
    if vmax is None:
        vmax = _np.amax(bc)
    sc = m.scatter(x,y,c=bc, linewidths=0, s=9, marker='s', vmin=vmin, vmax=vmax, rasterized=True, cmap=_cm.YlOrRd, ax=ax)


    if rafmt=='h':
        m.drawmeridians(meridians,labels=[0,0,0,1], fmt=_Lon2RA, linewidth=0.5)
    elif rafmt=='d':
        m.drawmeridians(meridians,labels=[0,0,0,1], linewidth=0.5)

    if raflip:
        m.drawparallels(parallels,labels=[0,1,0,0], labelstyle="+/-", linewidth=0.5, xoffset=(r*_np.radians(xoffset)))
        ax.invert_xaxis()
    else:
        m.drawparallels(parallels,labels=[1,0,0,0], labelstyle="+/-", linewidth=0.5, xoffset=(r*_np.radians(xoffset)))


    cb = m.colorbar(sc,"right", size="3%", pad='0%')
    if clabel is not None:
        cb.set_label(clabel)
    #cb.set_ticklabels(_np.linspace(vmin, vmax, (vmax - vmin)/1 + 1, dtype='int32'))
    cb.solids.set_edgecolor("face")


def MapPlot(ax=None, fig=None, nside=512, cat=None, ra=None, dec=None, nest=False, parallels=_np.arange(-75.,0.,5.), meridians=_np.arange(0.,360.,5.), dims=[20,20], center=[-75,-52.5], vmin=None, vmax=None, clabel='$n_s\ [\mathrm{arcmin}^{-2}]$', rafmt='d', raflip=True, xoffset=0 ):

    _BasePlot(ax=ax, fig=fig, nside=nside, cat=cat, ra=ra, dec=dec, nest=nest, parallels=parallels, meridians=meridians, dims=dims, center=center, vmin=vmin, vmax=vmax, clabel=clabel, rafmt=rafmt, raflip=raflip, xoffset=xoffset, f=_getCountLocation, extrakwargs={} )


def MapValPlot(ax=None, fig=None, cat=None, ra=None, dec=None, nest=False, parallels=_np.arange(-75.,0.,5.), meridians=_np.arange(0.,360.,5.), dims=[20,20], center=[-75,-52.5], vmin=None, vmax=None, clabel='$n_s\ [\mathrm{arcmin}^{-2}]$', rafmt='d', raflip=True, xoffset=0, map=None ):
    nside = _hp.npix2nside(map.size)
    _BasePlot(ax=ax, fig=fig, nside=nside, cat=cat, ra=ra, dec=dec, nest=nest, parallels=parallels, meridians=meridians, dims=dims, center=center, vmin=vmin, vmax=vmax, clabel=clabel, rafmt=rafmt, raflip=raflip, xoffset=xoffset, f=_getMapLocation, extrakwargs={'map':map} )

