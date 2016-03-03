import fitsio
import healpy as hp
import numpy as np


class SLRShift:
    '''
    Class for applying SLR shifted zeropoints for Y1A1
    '''

    def _Print(self, msg, level='info'):
        if self.balrogprint is not None:
            self.balrog.SysInfoPrint(self.balrogprint, msg, level=level)
        else:
            print msg


    def __init__(self, zpfile, fill_periphery=True, balrogprint=None):
        '''
        Input: name of SLR map FITS file
        '''
        
        self.balrogprint = balrogprint
        if self.balrogprint is not None:
            import balrog
            self.balrog = balrog

        self.zpfile = zpfile

        zpshifts,hdr=fitsio.read(zpfile,ext=1,header=True)
        self.nside = hdr['NSIDE']
        self.zpshifts = zpshifts
        self.fill_periphery = fill_periphery

        # translate mags and indices...
        self.nmag = self.zpshifts['ZP_SHIFT'][0,:].size
        self.bands = np.zeros(self.nmag,dtype='S1')

        for i in xrange(self.nmag):
            magname = hdr['MAG%d' % (i)]
            parts=magname.split('_')
            # force to be upper case...
            self.bands[i] = parts[len(parts)-1].upper()  

        self._Print("Using zpshift file: %s" % (self.zpfile))
        self._Print("Bands: %s" % (str(self.bands)))
        self._Print("NSIDE = %d (~%.1f arcmin)" % (self.nside, hp.nside2resol(self.nside, True)))
        self._Print("Area = %.2f deg^2" % (hp.nside2pixarea(self.nside, degrees=True) * self.zpshifts.size))

        # now we want to build the healpix map...
        self.shiftmap = np.zeros((self.nmag,hp.nside2npix(self.nside)),dtype=np.float32) + hp.UNSEEN
        for i in xrange(self.nmag):
            # all 99s are un-fit
            gd,=np.where(self.zpshifts['ZP_SHIFT'][:,i] < 90.0)
            self.shiftmap[i,self.zpshifts['HPIX'][gd]] = self.zpshifts['ZP_SHIFT'][gd,i]
            
        if (self.fill_periphery):
            self._Print("Filling in periphery pixels...")
            self._fill_periphery()

            self._Print("Peripheral area:")
            for i in xrange(self.nmag):
                self._Print("  %s: %.2f deg^2" % (self.bands[i], self.peripheral_area[i]))
        else:
            self.peripheral_area = None

    def _fill_periphery(self):
        '''
        Fill in peripheral cells of shiftmap with average of nearest neighbors
        '''

        self.peripheral_area = np.zeros(self.nmag,dtype=np.float32)
        for i in xrange(self.nmag):
            all_neighbor_pix = np.unique(hp.get_all_neighbours(self.nside,np.nonzero(self.shiftmap[i,:] != hp.UNSEEN)[0]))

            # remove any negatives...
            gd,=np.where(all_neighbor_pix >= 0)
            all_neighbor_pix = all_neighbor_pix[gd]

            filled_pix, = np.nonzero(self.shiftmap[i,:] != hp.UNSEEN)
            periphery_pix = np.setdiff1d(all_neighbor_pix, filled_pix)
            
            shiftmap_filled = np.ma.masked_array(self.shiftmap[i,:], self.shiftmap[i,:] == hp.UNSEEN)
            self.shiftmap[i,periphery_pix] = np.array(np.mean(shiftmap_filled[hp.get_all_neighbours(self.nside,periphery_pix)],axis=0),dtype=shiftmap_filled.dtype)

            self.peripheral_area[i] = hp.nside2pixarea(self.nside, degrees=True) * periphery_pix.size
            
    def get_zeropoint_offset(self, band, ra, dec, interpolate=True):
        '''
        Inputs: band {g,r,i,z,Y}, arrays for ra (deg), dec (deg)
        Return: array of zeropoint offsets to be added
        Note: returns 99 for bad fits outside region
        '''

        # get band index
        bind,=np.where(band.upper() == self.bands)
        if (bind.size == 0):
            self._Print("Error: could not find band %s in list of bands %s" % (band.upper(), str(self.bands)), level='warning')
            return np.zeros(ra.size,dtype=np.float32) + 99.0
        bind=bind[0]
        
        # get theta/phi
        theta = np.radians(90.0-dec)
        phi = np.radians(ra)

        # interpolate zeropoints (if necessary)
        if (interpolate):
            # use masked array here...
            m=np.ma.array(self.shiftmap[bind,:],mask=(self.shiftmap[bind,:] == hp.UNSEEN))
            zp_shift = np.ma.filled(hp.get_interp_val(m,theta,phi),hp.UNSEEN)
        else :
            pix = hp.ang2pix(self.nside, theta, phi)
            zp_shift = self.shiftmap[bind,pix]

        bd,=np.where(np.abs(zp_shift) > 10.0)
        if (bd.size > 0):
            zp_shift[bd] = 99.0

        return zp_shift

    def add_zeropoint_offset(self, band, ra, dec, mag, interpolate=True):
        '''
        Inputs: band {g,r,i,z,Y}, arrays for ra (deg), dec (deg)
        Return: array of zeropoint offsets to be added
        Note: returns 99 for bad fits outside region
        '''

        zp_shift = self.get_zeropoint_offset(band, ra, dec, interpolate=interpolate)

        return mag + zp_shift
