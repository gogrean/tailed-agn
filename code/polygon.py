import astropy.io.fits as fits
import astropy.wcs as wcs

class Poly(object):
    '''
    Polygon region that traces the 3-sigma contour level around a radio relic.

    The object is created from a DS9 region file saved in J2000 coordinates (degrees).
    '''
    def __init__(self, reg_file):
        '''Initialze RA and DEC coordinates by parsing them from the DS9 region file.'''
        with open(reg_file, 'r') as f:
            ds9_poly_reg = f.readlines()[-1]
        coords_start = ds9_poly_reg.find("(") + 1
        coords_end = ds9_poly_reg.find(")")
        coords_def = ds9_poly_reg[coords_start:coords_end]
        coords = coords_def.split(",")
        self.ra = [float(x) for x in coords[::2]]
        self.dec = [float(x) for x in coords[1::2]]

    @staticmethod
    def get_wcs(img_file, ext=0):
        hdulist = fits.open(img_file)
        w = wcs.WCS(hdulist[ext].header, hdulist)
        hdulist.close()
        return w

    def to_pixels(self, img, ext=0):
        '''
        Convert polygon coordinates from degrees to pixels.

        The function requires the filename of an image whose header is used to create
        the WCS object needed for the coordinate transformation. By default, the 0
        extension of the file is used for the header. A different extension can be
        provided by defining the ext attribute.
        '''
        self.wcs = self.get_wcs(img, ext=ext)
        return self.wcs.all_world2pix(self.ra, self.dec, 1)

    def to_degrees(self, x, y, img, ext=0):
        '''
        Convert polygon coordinates from pixels to degrees.

        The function converts two lists of x and y pixel coordinates to degrees using
        a given WCS object. By default, the 0 extension of the file is used for the
        WCS. A different extension can be provided by defining the ext attribute.
        '''
        self.wcs = self.get_wcs(img, ext=ext)
        return self.wcs.all_pix2world(x, y, 1)
