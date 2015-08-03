"""
=======================
MockGFA
=======================

This produced mock ETCGuiderPack objects using
galsim to produce fake stars

"""

from __future__ import division
from __future__ import print_function
from astropy import units as u
import numpy as np
import numpy.random as npran
import galsim
from Telescope import Telescope
from MockConditions import MockConditions
from ETCGuiderPack import ETCGuiderPack

class MockGFA(object):
    def __init__(self, Nstars=10, Nskies=10, star_rmags=18.0, Npix=30, etime=0.5*u.s, telescope=Telescope(), 
                 cond=MockConditions(), Tccdnoise=0.1, Tccdoffset=0.3,
                 filter_midpt=645*u.nm, filter_bandwidth=0.16*u.um, filter_transmission=0.9,
                 QE=0.9,pixel_size=15*u.um, rms_read_noise=20.0, poisson=True):
        """
        Constructor with the following options:
        Nstars = number of postage stamps with stars to mock
        Nskies = number of postage stamps with sky to mock
        star_rmags=r band magnitudes of stars. Either a number (all stars the same)
                   or a list.
        Nside = side of the postage stamp image
        etime = exposure time
        telescope = telescope model (for throughput, etc.)
        cond = conditions object
        Tccdnoise = rms noise in K of CCD temperature
        Tccdoffset = offset in K of CCD temperature
        filter_midp = GFA CCD filter mid point
        filter_bandwidth = GFA CCD filter bandwidth 
        filter_transmission = 1-fraction of light lost in the filter
        QE = quantum efficiency
        pixels_size = GFA pixel size
        rms_read_noise = readout noise in electrons
        poisson = if true, make real poisson noise realization, otherwise use gaussian of same variance
        Note: dark current is determined from CCD temperature as given by the conditions object

        """
        self.Nstars=Nstars
        self.Nskies=Nskies
        self.Npix=Npix
        self.etime=etime
        self.cond=cond
        self.telescope=telescope
        self.Tccdnoise=Tccdnoise
        self.Tccdoffset=Tccdoffset
        self.filter_midpt = filter_midpt
        self.throughput = QE*filter_transmission*filter_bandwidth
        self.pixel_size = pixel_size
        self.rms_read_noise = rms_read_noise
        self.poisson = poisson
        try:
            if (len(star_rmags)==Nstars):
                self.star_rmags=star_rmags
            else:
                print("List of magnitudes different from # of stars")
                stop()
        except TypeError:
            self.star_rmags=[star_rmags]*Nstars
        self.EffArea = telescope.effective_area*telescope.throughput*self.throughput
        

    def MockOnePostage (self, rmag,atmospheric_transmission,FWHM5400, sky, dark_current):
        do_star= not (rmag==None)
        if do_star:
            src_photons = (6.25e10/(u.m**2*u.s*u.um))*10**(-0.4*rmag)*self.etime
            src_electrons = src_photons*self.EffArea*atmospheric_transmission
            assert src_electrons == src_electrons.value,(
                'Source electrons not dimensionless: %r' % src_electrons)
        pixel_size = self.pixel_size/self.telescope.plate_scale
        sky_electrons_per_pixel = sky*self.etime*self.EffArea*pixel_size**2
        assert sky_electrons_per_pixel == sky_electrons_per_pixel.value,(
                'Skye electrons not dimensionless: %r' % src_electrons)
        
        # Calculate the mean image centered in a ROI, with a random sub-pixel offset.
        pixel_size_arcsec = pixel_size.to(u.arcsec).value
        ROI = galsim.Image(self.Npix,self.Npix,scale = pixel_size_arcsec)
        if do_star:
            dx0,dy0 = pixel_size_arcsec*np.random.uniform(low = -1,high = +1,size = 2)
            PSF = src_electrons.value*self.telescope.get_psf(self.filter_midpt,FWHM5400).shift(dx = dx0,dy = dy0)
            PSF.drawImage(image = ROI)
        #Add sky noise and dark current        
        extra_el=sky_electrons_per_pixel + dark_current*self.etime

        ROI +=extra_el
        if (self.poisson):
            ROI.addNoise(galsim.PoissonNoise(sky_level=float(extra_el)))
        else:
            ROI.array+= npran.normal(0.0,sqrt(ROI.array))
            
        ROI.addNoise(galsim.GaussianNoise(sigma=self.rms_read_noise))
        return ROI
        
        

    def Yield (self, t):
        """
        Returns a GFAGuiderPack, based on options
        that were set up above
        """
        atmospheric_transmission=self.cond.atmospheric_transmission(t)
        FWHM5400=self.cond.FWHM5400(t)
        sky=self.cond.sky(t)
        
        Tccd=self.cond.Tccd(t)/(u.K)
        dark_current = 19335.*Tccd**3*np.exp(-6400./Tccd)/u.s # per pixel per second
        stars=[self.MockOnePostage(rmag,atmospheric_transmission,FWHM5400, sky, dark_current) for rmag in self.star_rmags]
        skies=[self.MockOnePostage(None,atmospheric_transmission,FWHM5400, sky, dark_current) for c in range(self.Nskies)]
        Tccd_measured=Tccd+npran.normal(self.Tccdoffset,self.Tccdnoise)

        return ETCGuiderPack(stars, skies, Tccd)

