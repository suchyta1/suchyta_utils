"""
The :mod:`suchyta_utils.slr` submodule is for applying the DES stellar locus regression (SLR) to the DES (or Balrog) data.
Both Y1A1 and SVA1 are supported.

.. note::
    DESers can download the SLR python scripts and FITS files from the 
    `supplementary directory <http://www.physics.ohio-state.edu/~suchyta.1/suchyta_utils/supplementary/files/>`_.
    ``y1a1_slr_shiftmap.py`` is the same file as available `Y1 Redmine page <https://cdcvs.fnal.gov/redmine/projects/des-y1/wiki/SLR_Adjustment_of_Y1A1>`_, 
    but ``sva1_slr_shiftmap.py`` is a modified version of what's available from the 
    `SVA1 wiki pages <https://cdcvs.fnal.gov/redmine/projects/descalibration/wiki/Stellar_Locus_Regression_for_SVA1_on_HEALPix_Grid>`_. 
    I changed some function names and returns to be consistent with the Y1 version. 
    In the future I should probably write my code to use the SVA1 file verbatim from Redmine, but that's not the case yet.

"""

import os as _os
import imp as _imp
import numpy as _np


class SLR:
    """
    Instantiate an SLR object. Both the FITS file and the python file for the release must both live in `slrdir`.
    One must also make sure to have the appropriate FITS file for Y1, between wide, D04, D10, and DFull.

    Parameters
    ----------
    release (str)
        'y1a1' or 'sva1' are valid
    area (str)
        Only relevant with ``release = y1a1``. Which dataset the SLR is for. Valid choices are ['wide', 'd04', 'd10', 'dfull']
    slrdir (str)
        Directory for the SLR FITS files and python scripts.

    Returns
    -------
    slr (object)
        An SLR object whose methods can be called to get corrections.

    """

    def __init__(self, release='y1a1', area='wide', slrdir=None):

        self.slrdir = slrdir
        if self.slrdir is None:
            self.slrdir = _os.path.dirname(_os.path.realpath(__file__))

        self.slrfile = '%s_slr_shiftmap.py' %(release)
        _slr = _imp.load_source('_slr', _os.path.join(self.slrdir,self.slrfile))


        if release=='y1a1':
            self.slrfits = _os.path.join(slrdir, 'y1a1_%s_slr_wavg_zpshift2.fit'%(area))
            self.slrshift = _slr.SLRShift(self.slrfits, fill_periphery=True)
        elif release=='sva1':
            self.slrfits = _os.path.join(slrdir, 'slr_zeropoint_shiftmap_v6_splice_cosmos_griz_EQUATORIAL_NSIDE_256_RING.fits')
            self.slrshift = _slr.SLRZeropointShiftmap(self.slrfits, fill_periphery=True)



    def GetMagShifts(self, band, ra, dec):
        """
        For each RA/DEC pair, find the SLR offset to be applied to the uncorrected magnitude measurements

        Parameters
        ----------
        band (str)
            DES pass band. Possible values are ['g','r','i','z','y'], but y-band was not computed for SV
        ra (float array)
            RA coordinates
        dec (float array)
            DEC coordinates

        Returns
        -------
        mag (float array)
            The shifts to be added to raw magnituded meausrements

        """
        return self.slrshift.get_zeropoint_offset(band, ra, dec, interpolate=True)


    def GetFluxFactors(self, band, ra, dec):
        """
        For each RA/DEC pair, find the SLR multiplicative to be applied to the flux measurements

        Parameters
        ----------
        band (str)
            DES pass band. Possible values are ['g','r','i','z','y'], but y-band was not computed for SV
        ra (float array)
            RA coordinates
        dec (float array)
            DEC coordinates

        Returns
        -------
        f (float array)
            The factors to multiply the raw flux meausrements by

        """
        offsets = self.GetMagShifts(band, ra, dec)
        return _np.power(10.0, -offsets/2.5)

