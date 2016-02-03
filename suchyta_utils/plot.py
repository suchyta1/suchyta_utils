""" 
(In my mind) :mod:`suchyta_utils.plot` makes plots looks nicer by default
and makes plotting maps in python more user-friendly.

Here's an example snippet from the git repository's `example.py <https://github.com/suchyta1/suchyta_utils/blob/master/examples/example.py>`_,
which shows some basic usage::

    import suchyta_utils as es
    import esutil
    import numpy as np
    import os
    import matplotlib.pyplot as plt

    # Setup automatic gridlines, among a few other changes to make the plots more readable
    es.plot.Setup()

    supdir = './supplementary'
    y1_data_file = os.path.join(supdir, 'balrog-small-y1.fits')
    y1_data = esutil.io.read(y1_data_file)

    # The plot submodule has a wrapper for making Basemap more user friendly to plot in projection
    map_fig, map_ax = plt.subplots(1,2, figsize=(16,4))
    es.plot.MapPlot(ax=map_ax[0], fig=map_fig, cat=y1_data, ra='alphawin_j2000_i', dec='deltawin_j2000_i', parallels=np.arange(-60.,0.,2.), meridians=np.arange(0.,360.,5.), dims=[20,5], clabel=r'$n_s\ [\mathrm{arcmin}^{-2}]$', vmin=10, vmax=30, center=[72,-57], nside=1024, xoffset=0.6)
    map_ax[0].set_title('Y1 Map')
"""

import matplotlib as _mpl
import matplotlib.pyplot as _plt
import os as _os
import suchyta_utils as _es
import healpy as _hp
import numpy as _np

from mpl_toolkits.basemap import Basemap as _Basemap
import matplotlib.cm as _cm


def Setup():
    """
    Setup some automatic styling of matplotlib plots. This started from
    a Supermongo-like style file from Matt Becker, that I probably tweaked to make things look nice to me.
    It loads `custom-style.mpl <https://github.com/suchyta1/suchyta_utils/blob/master/suchyta_utils/custom-sytle.mpl>`_, so look there to see all the settings. 
    The most notable thing that I wanted was drawing grid lines.

    .. note::
        I'll add a keyword to use any style file formatted like the default one.

    Parameters
    ----------
    None

    Returns
    -------
    None

    """

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
    """
    Set the maximum number of ticks along the axes. 
    This uses :mod:`matplotlib.ticker.MaxNLocator`.

    Parameters
    ----------
    ax (matplotlib axis object)
        The axis object containing the plot to format
    nxticks (int/None)
        The maximum number of ticks on the xaxis. None means keep whatever matplotlib does by default.
    nyticks (int/None)
        The maximum number of ticks on the yaxis. None means keep whatever matplotlib does by default.

    Returns
    -------
    None

    """

    if nyticks is not None:
        ax.yaxis.set_major_locator( _mpl.ticker.MaxNLocator(nbins=(nyticks-1)) )
        ax.yaxis.set_major_locator( _mpl.ticker.MaxNLocator(nbins=(nyticks-1)) )
    if nxticks is not None:
        ax.xaxis.set_major_locator( _mpl.ticker.MaxNLocator(nbins=(nxticks-1)) )
        ax.xaxis.set_major_locator( _mpl.ticker.MaxNLocator(nbins=(nxticks-1)) )


class Offset():

    def __init__(self, thing, by, log=False):
        if log:
            self.new = _np.log10(thing) + by
            self.new = _np.power(10, self.new)
        else:
            self.new = thing + by

    def __enter__(self):
        return self.new

    def __exit__(self, type, value, traceback):
        pass


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
    """
    Draw a line segement on a plot

    Parameters
    ----------
    ax (matplotlib axis object)
        The axis object containing the plot to format
    left (array)
        Starting coordinate for the line segment
    right (array)
        Ending coordinate for the line segment
    plotkwargs (dict)
        Dictionary of keyword-value pairs to give as `**kwargs` to ``ax.plot()`` for drawing the line segment


    Returns
    -------
    ax (matplotlib axis object)
        The modified axis object
        
    """
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


def _BasePlot(ax=None, fig=None, nside=512, cat=None, ra=None, dec=None, nest=False, parallels=_np.arange(-75.,0.,5.), meridians=_np.arange(0.,360.,5.), dims=[20,20], center=[75,-52.5], vmin=None, vmax=None, clabel='$n_s\ [\mathrm{arcmin}^{-2}]$', rafmt='d', raflip=True, xoffset=0, size=9, f=_getCountLocation, noplot=False, extrakwargs={} ):
    if fig is None:
        fig, ax = plt.subplots(1,1, figsize=(6.5*nside/512,6*nside/512))
    if ax is None:
        fig, ax = _plt.subplots(1,1)

    r = 6.371e6
    h = r * _np.radians(dims[1])
    w = r * _np.cos(_np.radians(center[1])) * _np.radians(dims[0])

    m = _Basemap(projection='aea', width=w, height=h, lat_0=center[1], lon_0=center[0], ax=ax)
    if rafmt=='h':
        m.drawmeridians(meridians,labels=[0,0,0,1], fmt=_Lon2RA, linewidth=0.5)
    elif rafmt=='d':
        m.drawmeridians(meridians,labels=[0,0,0,1], linewidth=0.5)

    if raflip:
        m.drawparallels(parallels,labels=[0,1,0,0], labelstyle="+/-", linewidth=0.5, xoffset=(r*_np.radians(xoffset)))
        ax.invert_xaxis()
    else:
        m.drawparallels(parallels,labels=[1,0,0,0], labelstyle="+/-", linewidth=0.5, xoffset=(r*_np.radians(xoffset)))


    if noplot:
        return m
    bc, lat, lon = f(cat=cat, ra=ra, dec=dec, nside=nside, nest=nest, **extrakwargs)
    x,y  = m(lon, lat)

    if vmin is None:
        vmin = _np.amin(bc)
    if vmax is None:
        vmax = _np.amax(bc)

    sc = m.scatter(x,y,c=bc, linewidths=0, s=size, marker='s', vmin=vmin, vmax=vmax, rasterized=True, cmap=_cm.YlOrRd, ax=ax)
    cb = m.colorbar(sc,"right", size="3%", pad='0%')
    if clabel is not None:
        cb.set_label(clabel)
    #cb.set_ticklabels(_np.linspace(vmin, vmax, (vmax - vmin)/1 + 1, dtype='int32'))
    cb.solids.set_edgecolor("face")


def MapPlot(ax=None, fig=None, nside=512, cat=None, ra=None, dec=None, nest=False, parallels=_np.arange(-75.,0.,5.), meridians=_np.arange(0.,360.,5.), dims=[20,20], center=[75,-52.5], vmin=None, vmax=None, clabel='$n_s\ [\mathrm{arcmin}^{-2}]$', rafmt='d', raflip=True, xoffset=0, size=9, noplot=False ):
    """
    Make a number density plot, showing the map in equal-area projection.
    This function uses :mod:`mpl_toolkits.basemap`, but I find that syntax terribly non-convenient and prefer the syntax here.
   
    .. note::
        Number densities to plot are computed in HEALPix cells, and each HEALPix cell will be drawn as a square,
        so setting a higher `nside` plots more points. The size of thise points is not currently adjustable,
        so things can end up looking weird for certain combinations of axis ratios, areas covered, and `nsides`.

    .. note::
        `cat`, `ra`, and `dec` behave the same way as they do in :mod:`suchyta_utils.hp`.

    Parameters
    ----------
    ax (matplotlib axis object/None)
        The axis object containing the plot to format. None will make a new axis.
    fig (matplotlib figure object)
        The figure containing the axis. None will make a new figure.
    nside (int)
        HEALPix nside to used for determining the number densities to plot
    cat (None/structured array)
        If not None, the structured data array (e.g. numpy recarray) whose number density is to be plotted
    ra (float array/str)
        If `cat` is None, an array of the RA values for the objects. Otherwise, the column name for the RA column in `cat`.
    dec (float array/str)
        if `cat` is None, an array of the DEC values for the objects. Otherwise, the column name for the DEC column in `cat`.
    nest (bool)
        Whether or not to use nested format for the plot. 
    parallels (array)
        Parallels to draw on the plot
    meridians (array)
        Meridians to draw on the plot
    dims (array)
        Approximate dimensions of the plot, in degrees. The first element is longitudonal width, and the second element is the latitudonal height.
    center (array)
        Approximate center of the plot, in degrees.The first element is longitudonal width, and the second element is the latitudonal height.
    vmin (float)
        Minimum number density (in armin^-2) for the color scale
    vmax (float)
        Maximum nunber density (in armin^-2) for the color scale
    clabel (str)
        Label for the color bar
    rafmt (str)
        How to label the meridians on the plot. Possible values are ``['d','h']``. ``'d'`` means degree format, and ``'h'`` means hour format
    raflip (bool)
        Whether or not to plot the larger RA values on the left
    xoffset (float)
        Offset for the labels of the paralles, in units of degrees on the plot. 

    Returns
    -------
    None

    """

    return _BasePlot(ax=ax, fig=fig, nside=nside, cat=cat, ra=ra, dec=dec, nest=nest, parallels=parallels, meridians=meridians, dims=dims, center=center, vmin=vmin, vmax=vmax, clabel=clabel, rafmt=rafmt, raflip=raflip, xoffset=xoffset, f=_getCountLocation, size=size, noplot=noplot, extrakwargs={} )


def MapValPlot(ax=None, fig=None, cat=None, ra=None, dec=None, nest=False, parallels=_np.arange(-75.,0.,5.), meridians=_np.arange(0.,360.,5.), dims=[20,20], center=[75,-52.5], vmin=None, vmax=None, clabel='$n_s\ [\mathrm{arcmin}^{-2}]$', rafmt='d', raflip=True, xoffset=0, size=9, map=None ):
    """
    Make a plot of a HEALPix map, showing the map in equal-area projection.
    This function uses :mod:`mpl_toolkits.basemap`, but I find that syntax terribly non-convenient and prefer the syntax here.
   
    .. note::
        This function plots the map in whatever resolution the map has.
        Things can end up looking weird for certain combinations of axis ratios, areas covered, and `nsides`.

    .. note::
        `cat`, `ra`, `dec`, and `map` behave the same way as they do in :mod:`suchyta_utils.hp`.

    Parameters
    ----------
    map (array)
        The HEALPix map. This amounts to an array, ordered by the HEALPix index.
    ax (matplotlib axis object/None)
        The axis object containing the plot to format. None will make a new axis.
    fig (matplotlib figure object)
        The figure containing the axis. None will make a new figure.
    cat (None/structured array)
        If not None, the structured data array (e.g. numpy recarray) whose number density is to be plotted
    ra (float array/str)
        If `cat` is None, an array of the RA values for the objects. Otherwise, the column name for the RA column in `cat`.
    dec (float array/str)
        if `cat` is None, an array of the DEC values for the objects. Otherwise, the column name for the DEC column in `cat`.
    nest (bool)
        Whether or not to use nested format for the plot. 
    parallels (array)
        Parallels to draw on the plot
    meridians (array)
        Meridians to draw on the plot
    dims (array)
        Approximate dimensions of the plot, in degrees. The first element is longitudonal width, and the second element is the latitudonal height.
    center (array)
        Approximate center of the plot, in degrees.The first element is longitudonal width, and the second element is the latitudonal height.
    vmin (float)
        Minimum number density (in armin^-2) for the color scale
    vmax (float)
        Maximum nunber density (in armin^-2) for the color scale
    clabel (str)
        Label for the color bar
    rafmt (str)
        How to label the meridians on the plot. Possible values are ``['d','h']``. ``'d'`` means degree format, and ``'h'`` means hour format
    raflip (bool)
        Whether or not to plot the larger RA values on the left
    xoffset (float)
        Offset for the labels of the paralles, in units of degrees on the plot. 

    Returns
    -------
    None

    """
    nside = _hp.npix2nside(map.size)
    return _BasePlot(ax=ax, fig=fig, nside=nside, cat=cat, ra=ra, dec=dec, nest=nest, parallels=parallels, meridians=meridians, dims=dims, center=center, vmin=vmin, vmax=vmax, clabel=clabel, rafmt=rafmt, raflip=raflip, xoffset=xoffset, f=_getMapLocation, size=9, extrakwargs={'map':map} )


