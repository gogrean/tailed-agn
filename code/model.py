import numpy as np
import astropy.io.fits as fits
from astropy.cosmology import FlatLambdaCDM

from polygon import Poly

class Model(object):
    '''
    Model of the galaxy distribution expected in a cluster.
    '''
    def __init__(self, gal_file, ext=0):
        '''
        Instantiate image and coordinate boundaries used to calculate
        the model.
        '''
        self.gal_file = gal_file
        hdulist = fits.open(self.gal_file)
        self.img = hdulist[ext].data
        self.hdr = hdulist[ext].header
        # might change this if we decide to limit the coordinates
        # 1 is subtracted from everything since count starts at 0 in Python
        self.xmin, self.ymin = 0, 0
        self.xmax = hdulist[ext].header['NAXIS2'] - 1
        self.ymax = hdulist[ext].header['NAXIS1'] - 1

    def histogram(self, polygon, z, bin_edges=range(0,3000,50)):
        '''
        Make histogram of the expected distribution of galaxies with respect
        to distance (in kpc) from the radio relics. The radio relics are
        defined by a polygon object. The redshift is needed to convert
        distances from arcsec to kpc.
        '''
        y_rel, x_rel = polygon.to_pixels(self.gal_file)
        cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
        scale = cosmo.kpc_proper_per_arcmin(z).value
        pix2arcsec = self.hdr['CD1_1'] * 60. * scale
        dist = []
        for x in range(self.xmin, self.xmax + 1):
            for y in range(self.ymin, self.ymax + 1):
                d = np.min(np.sqrt((x-x_rel)**2 + (y-y_rel)**2)) * pix2arcsec
                # converting galaxy density from per arcmin**2 to per kpc**2
                # that should be okay when adding up histograms???
                # REINOUT?
                dist.extend([d] * int(np.round(self.img[x, y] * 1e5 / scale**2)))
        # returns hist, edges
        return np.histogram(dist, bins=bin_edges)
