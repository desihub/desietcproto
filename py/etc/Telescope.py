"""
=======================
Telescope
=======================

Telescope Object copied from workbook

"""
from __future__ import division
from __future__ import print_function
import numpy as np
from astropy import units as u
import galsim

class Telescope(object):
    """
    Represents a telescope.
    """
    def __init__(self,diameter=3.80*u.m,obscuration_area_fraction=0.25,throughput=0.95*0.77,plate_scale=67.40*u.um/u.arcsec):
        self.diameter = diameter
        self.obscuration_area_fraction = obscuration_area_fraction
        self.throughput = throughput
        self.plate_scale = plate_scale
        self.effective_area = np.pi*diameter**2/4.*(1-obscuration_area_fraction)
    def get_optical_psf(self,wavelength):
        #Convert dimensionless lam/D to arcsec units.
        lam_over_diam_arcsec = ((wavelength/self.diameter)*u.rad).to(u.arcsec)
        # Airy requires floats as inputs, not numpy scalars.
        return galsim.Airy(lam_over_diam=float(lam_over_diam_arcsec.value),
            obscuration=float(np.sqrt(self.obscuration_area_fraction)))
    def get_atmospheric_psf(self,wavelength,fwhm5400):
        wlen_ratio = (wavelength/(5400*u.Angstrom)).si
        assert wlen_ratio == wlen_ratio.value,'wavelength has invalid units.'
        fwhm = fwhm5400.to(u.arcsec).value*wlen_ratio**(-0.2)
        # Kolmogorov requires floats as inputs, not numpy scalars.
        return galsim.Kolmogorov(fwhm=float(fwhm))
    def get_psf(self,wavelength,fwhm5400,rms_jitter=0.1*u.arcsec):
        components = [ self.get_atmospheric_psf(wavelength,fwhm5400),self.get_optical_psf(wavelength) ]
        # Include a Gaussian pointing jitter, if requested.
        if rms_jitter is not None:
            components.append(galsim.Gaussian(sigma = rms_jitter.to(u.arcsec).value))
        return galsim.Convolve(components)
    def plot_fwhm(self,wlen1,wlen2,fwhm5400,nwlen=50):
        wlen = np.linspace(wlen1,wlen2,nwlen)
        optical = np.empty_like(wlen)
        atmospheric = np.empty_like(wlen)
        for i in range(nwlen):
            print(wlen[i])
            optical[i] = self.get_optical_psf(wlen[i]).getFWHM()
            atmospheric[i] = self.get_atmospheric_psf(wlen[i]).getFWHM()
        plt.plot(wlen.value,optical,'r-',label='Optical')
        plt.plot(wlen.value,atmospheric,'b-',label='Atmospheric')
        plt.legend()
        plt.xlabel('Wavelength (%s)',wlen.unit)
        plt.ylabel('FWHM (arcsec)')


