#!/usr/bin/env python

"""
This file is intended to show some of the most useful functionality.
It isn't comprehensive, but shows many of the relevant functions and the arguments you might want to tweak.

The supplmentary files are available here. (You'll need the usual DES credentials):
    http://www.physics.ohio-state.edu/~suchyta.1/suchyta_utils/supplementary/files/

"""

import suchyta_utils as es
import esutil
import numpy as np
import os
import matplotlib.pyplot as plt


if __name__ == "__main__":

    # Setup automatic gridlines, among a few other changes to make the plots more readable
    es.plot.Setup()

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

    # The plot submodule has a wrapper for making Basemap more user friendly to plot in projection
    map_fig, map_ax = plt.subplots(1,2, figsize=(16,4))
    es.plot.MapPlot(ax=map_ax[0], fig=map_fig, cat=y1_data, ra='alphawin_j2000_i', dec='deltawin_j2000_i', parallels=np.arange(-60.,0.,2.), meridians=np.arange(0.,360.,5.), dims=[20,5], clabel=r'$n_s\ [\mathrm{arcmin}^{-2}]$', vmin=10, vmax=30, center=[72,-57], nside=1024, xoffset=0.6)
    map_ax[0].set_title('Y1 Map')


    # SV files
    sv_data_file = os.path.join(supdir, 'balrog-small-sv.fits')
    sv_goodmask_file = os.path.join(supdir, 'sva1_gold_1.0.4_goodregions_04_equ_ring_4096.fits.gz')

    # Apply the SVA1 4% exclusion mask to the Balrog data. Do everything with file names.
    sv_data = esutil.io.read(sv_data_file)
    sv_ra, sv_dec = es.hp.ApplyMask(mask=sv_goodmask_file, ra=sv_data['alphawin_j2000_i'], dec=sv_data['deltawin_j2000_i'], val=1, cond='=')

    # The plot submodule has a wrapper for making Basemap more user friendly to plot in projection
    es.plot.MapPlot(ax=map_ax[1], fig=map_fig, cat=sv_data, ra='alphawin_j2000_i', dec='deltawin_j2000_i', parallels=np.arange(-60.,0.,2.), meridians=np.arange(0.,360.,5.), dims=[20,5], clabel=r'$n_s\ [\mathrm{arcmin}^{-2}]$', vmin=10, vmax=30, center=[72,-57], nside=1024, xoffset=0.6, rafmt='h')
    map_ax[1].set_title('SV Map')
    plt.tight_layout(w_pad=3, pad=3)



    # SLR correct the data
    # There are different SLR sets for Y1, that's what area='wide' is for. Could also have D04, D10, DFull
    y1_bands = ['g','r','i','z','y']
    y1_slr = es.slr.SLR(release='y1a1', area='wide', slrdir=supdir)
    for band in y1_bands:
        shift = y1_slr.GetMagShifts(band, y1_data['alphawin_j2000_i'], y1_data['deltawin_j2000_i'])
        y1_data['mag_auto_{0}'.format(band)] = y1_data['mag_auto_{0}'.format(band)] + shift

    # SV SLR didn't do Y-band
    sv_bands = ['g','r','i','z']
    sv_slr = es.slr.SLR(release='sva1', slrdir=supdir)
    for band in sv_bands:
        shift = sv_slr.GetMagShifts(band, sv_data['alphawin_j2000_i'], sv_data['deltawin_j2000_i'])
        sv_data['mag_auto_{0}'.format(band)] = sv_data['mag_auto_{0}'.format(band)] + shift


    # Find modest class for all objects 
    # See the caveat about Y1A1 version (http://www.physics.ohio-state.edu/~suchyta.1/suchyta_utils/doc/html/balrog.html#suchyta_utils.balrog.Modest)
    sv_modest = es.balrog.Modest(sv_data, release='sva1')
    y1_modest = es.balrog.Modest(y1_data, release='y1a1')

    sv_stars = sv_data[ sv_modest==2 ]
    y1_stars = y1_data[ y1_modest==2 ]


    # Plot the magnitude distributions
    rows = 2
    cols = int(np.ceil(float(len(y1_bands))/rows))
    fig, axarr = plt.subplots(rows, cols, figsize=(16,9), tight_layout=True)
    bins = np.arange(18, 26.5, 0.1)
    c = (bins[1:] + bins[:-1]) / 2.0
    axarr[-1][-1].axis('off')
    for i in range(len(y1_bands)):
        ii = i / cols
        jj = i % cols
        ax = axarr[ii][jj]

        h, b = np.histogram(y1_stars['mag_auto_{0}'.format(y1_bands[i])], bins=bins, density=True)
        ax.plot(c, np.log10(h), color='blue', label='Y1')
        h, b = np.histogram(sv_stars['mag_auto_{0}'.format(y1_bands[i])], bins=bins, density=True)
        ax.plot(c, np.log10(h), color='red', label='SV')
        ax.set_xlabel(r'\texttt{MAG\_AUTO\_%s}'%(y1_bands[i].upper()))
        ax.set_ylabel(r'$\log_{10}(p)$')

        # The number of ticks set is more like a maximum than an exact number
        es.plot.NTicks(ax, nxticks=6, nyticks=8)
        if i==0:
            ax.legend(loc='best')


    

    plt.show()
