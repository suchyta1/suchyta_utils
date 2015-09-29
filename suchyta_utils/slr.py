import os as _os
import imp as _imp
import numpy as _np


class SLR:

    def __init__(self, release='y1a1', area='wide', slrdir=None):
        self.slrdir = slrdir
        if self.slrdir is None:
            self.slrdir = _os.path.dirname(_os.path.realpath(__file__))

        self.slrfile = '%s_slr_shiftmap.py' %(release)
        _slr = _imp.load_source('_slr', _os.path.join(self.slrdir,self.slrfile))


        if release=='y1a1':
            self.slrfits = _os.path.join(slrdir, 'y1a1_%s_slr_wavg_zpshift2.fit'%(area))
            self.slrshift = _slr.SLRShift(self.slrfits, fill_periphery=True)
        elif reliase='sva1':
            self.slrfits = _os.path.join(slrdir, 'slr_zeropoint_shiftmap_v6_splice_cosmos_griz_EQUATORIAL_NSIDE_256_RING.fits'
            self.slrshift = _slr.SLRZeropintShift(self.slrfits, fill_periphery=True)



    def GetMagShifts(self, band, ra, dec):
        return self.slrshift.get_zeropoint_offset(band, ra, dec, interpolate=True)

    def GetFluxFactors(self, band, ra, dec):
        offsets = self.GetMagShifts(band, ra, dec)
        return _np.power(10.0, -offsets/2.5)

