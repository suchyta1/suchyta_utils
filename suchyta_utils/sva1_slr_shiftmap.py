import pyfits
import healpy
import numpy

'''
try:
    import pylab
    pylab.ion()
except:
    print 'Could not import pylab, no plotting functions will be available.'
'''

############################################################

class SLRZeropointShiftmap:
    '''
    Class for applying stellar locus regression (SLR) shifted zeropoints
    '''

    def __init__(self, zeropoint_shiftmap_file, fill_periphery=True):
        '''
        Input: name of SLR zeropoint map FITS file
        '''

        self.zeropoint_shiftmap = pyfits.open(zeropoint_shiftmap_file)[1]
        self.nside = healpy.npix2nside(len(self.zeropoint_shiftmap.data))
        self.fill_periphery = fill_periphery
        
        #print 'SLR zeropoint shiftmap = %s'%(zeropoint_shiftmap_file)
        #print 'Filters = %s'%(', '.join(self.zeropoint_shiftmap.data.names))
        #print 'NSIDE = %i (~%.1f arcmin)'%(self.nside, healpy.nside2resol(self.nside, True))
        #print 'Area = %.2f deg**2'%(healpy.nside2pixarea(self.nside, degrees=True) * numpy.sum(self.zeropoint_shiftmap.data.field('n_des') > 0))
        
        if self.fill_periphery:
            #print 'Filling in periphery pixels...'
            for filter in self.zeropoint_shiftmap.data.names:
                if 'n_' not in filter:
                    self._fillPeriphery(filter)
                    
            #print 'Peripheral Area = %.2f deg**2'%(self.peripheral_area)
        else:
            self.peripheral_area = 0.

    def _fillPeriphery(self, filter):
        '''
        Fill in peripheral cells of shiftmap with average of nearest neighbors.
        '''
        all_neighbor_pix = numpy.unique(healpy.get_all_neighbours(self.nside,
                                                                  numpy.nonzero(self.zeropoint_shiftmap.data.field(filter) != healpy.UNSEEN)[0]))
        filled_pix = numpy.nonzero(self.zeropoint_shiftmap.data.field(filter) != healpy.UNSEEN)[0]
        periphery_pix = numpy.setdiff1d(all_neighbor_pix, filled_pix)
        shiftmap_filled = numpy.ma.masked_array(self.zeropoint_shiftmap.data.field(filter),
                                                self.zeropoint_shiftmap.data.field(filter) == healpy.UNSEEN)
        self.zeropoint_shiftmap.data.field(filter)[periphery_pix] = numpy.array(numpy.mean(shiftmap_filled[healpy.get_all_neighbours(self.nside,
                                                                                                                                     periphery_pix)],
                                                                                           axis=0), dtype=shiftmap_filled.dtype)
        self.peripheral_area = healpy.nside2pixarea(self.nside, degrees=True) * len(periphery_pix)
        
    def plot(self, filter, ra=None, dec=None):
        '''
        Plot zoomed zeropoint shiftmap for chosen filter centered on (ra, dec) given in degrees,
        or simply Mollweide all-sky map by default.
        '''
        if ra is not None and dec is not None:
            healpy.gnomview(self.zeropoint_shiftmap.data.field(filter), rot=(ra, dec, 0))
        else:
            healpy.mollview(self.zeropoint_shiftmap.data.field(filter))

    def addZeropoint(self, filter, ra, dec, mag, interpolate=True, plot=False):
        '''
        Inputs: filter name {g, r, i, z}, arrays for ra (deg), dec (deg), and magnitude
        Return: array of corrected magnitudes, array quality flag
        Note: uses median fitted zeropoint shift for objects having (rad, dec) outside fitted or periphery pixels
        '''
        pix = healpy.ang2pix(self.nside, numpy.radians(90. - dec), numpy.radians(ra))

        # Interpolate zeropoints
        if interpolate:
            zeropoint_shift = healpy.get_interp_val(self.zeropoint_shiftmap.data.field(filter),
                                                    numpy.radians(90. - dec), numpy.radians(ra))
        else:
            zeropoint_shift = self.zeropoint_shiftmap.data.field(filter)[pix]

        # Quality flag
        quality_flag = self.zeropoint_shiftmap.data.field('n_des')[pix]
        quality_flag[quality_flag > 0] = 0
        quality_flag[quality_flag < 0] = 2
        if self.fill_periphery:
            quality_flag[numpy.logical_and(quality_flag == 2, numpy.fabs(zeropoint_shift) <= 10.)] = 1

        # Use median zeropoint shift for objects outside interpolated region 
        median_fitted_zeropoint_shift = numpy.median(zeropoint_shift[numpy.fabs(zeropoint_shift) <= 10.])
        zeropoint_shift[numpy.fabs(zeropoint_shift) > 10.] = median_fitted_zeropoint_shift
        
        return mag + zeropoint_shift, quality_flag

    def GetZeropoint(self, filter, ra, dec, mag, interpolate=True, plot=False):
        '''
        Inputs: filter name {g, r, i, z}, arrays for ra (deg), dec (deg), and magnitude
        Return: array of corrected magnitudes, array quality flag
        Note: uses median fitted zeropoint shift for objects having (rad, dec) outside fitted or periphery pixels
        '''
        pix = healpy.ang2pix(self.nside, numpy.radians(90. - dec), numpy.radians(ra))

        # Interpolate zeropoints
        if interpolate:
            zeropoint_shift = healpy.get_interp_val(self.zeropoint_shiftmap.data.field(filter),
                                                    numpy.radians(90. - dec), numpy.radians(ra))
        else:
            zeropoint_shift = self.zeropoint_shiftmap.data.field(filter)[pix]

        # Quality flag
        quality_flag = self.zeropoint_shiftmap.data.field('n_des')[pix]
        quality_flag[quality_flag > 0] = 0
        quality_flag[quality_flag < 0] = 2
        if self.fill_periphery:
            quality_flag[numpy.logical_and(quality_flag == 2, numpy.fabs(zeropoint_shift) <= 10.)] = 1

        # Use median zeropoint shift for objects outside interpolated region 
        median_fitted_zeropoint_shift = numpy.median(zeropoint_shift[numpy.fabs(zeropoint_shift) <= 10.])
        zeropoint_shift[numpy.fabs(zeropoint_shift) > 10.] = median_fitted_zeropoint_shift
        
        return zeropoint_shift, quality_flag


    def get_zeropoint_offset(self, filter, ra, dec, interpolate=True):
        '''
        Inputs: filter name {g, r, i, z}, arrays for ra (deg), dec (deg), and magnitude
        Return: array of corrected magnitudes, array quality flag
        Note: uses median fitted zeropoint shift for objects having (rad, dec) outside fitted or periphery pixels
        '''
        pix = healpy.ang2pix(self.nside, numpy.radians(90. - dec), numpy.radians(ra))

        # Interpolate zeropoints
        if interpolate:
            zeropoint_shift = healpy.get_interp_val(self.zeropoint_shiftmap.data.field(filter),
                                                    numpy.radians(90. - dec), numpy.radians(ra))
        else:
            zeropoint_shift = self.zeropoint_shiftmap.data.field(filter)[pix]

        # Use median zeropoint shift for objects outside interpolated region 
        median_fitted_zeropoint_shift = numpy.median(zeropoint_shift[numpy.fabs(zeropoint_shift) <= 10.])
        zeropoint_shift[numpy.fabs(zeropoint_shift) > 10.] = median_fitted_zeropoint_shift
        
        return zeropoint_shift

############################################################


'''
# Fill in peripheral cells of shiftmap with average of nearest neighbors
neighbor = numpy.take(self.zeropoint_shiftmap.data.field(filter),
healpy.get_all_neighbours(self.nside, range(0, healpy.nside2npix(self.nside))))
neighbor_mask = numpy.ma.masked_array(neighbor, neighbor == healpy.UNSEEN)

mean_neighbor = neighbor_mask.mean(axis=0)
indices = numpy.nonzero(numpy.logical_and(mean_neighbor > 0,
self.zeropoint_shiftmap.data.field(filter) == healpy.UNSEEN))[0]

neighbor_shiftmap = self.zeropoint_shiftmap.data.field(filter)        
neighbor_shiftmap[indices] = numpy.array(mean_neighbor[indices].compressed().tolist(), dtype=neighbor_shiftmap.dtype)
'''
# Plot
#if plot:
#    ra_center, dec_center = numpy.median(ra), numpy.median(dec)
#    healpy.gnomview(self.zeropoint_shiftmap.data.field(filter), rot=(ra_center, dec_center, 0))
#    healpy.gnomview(neighbor_shiftmap, rot=(ra_center, dec_center, 0))
