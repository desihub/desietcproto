"""
=======================
MockConditions
=======================

This produces mock changes in brightness, CCD temperature, etc.

"""

from __future__ import division
from __future__ import print_function
from astropy import units as u


class MockConditions(object):
    def __init__(self, sky=546/(u.m**2*u.s*u.um*u.arcsec**2), 
                 Tccd=293.*u.K, 
                 atmospheric_transmission=0.76,
                 FWHM5400=1*u.arcsec,
                 dsky_dt=0.0, dTccd_dt=0.0, datmospheric_transmission_dt=0.0, dFWHM5400_dt=0.0):
        """
        Mocks variable condition inside the dome/sky. Really a wrapper object
        to allow study of quick changes in conditions isolated in a simple object.
        sky = sky brightness
        Tccd = CCD temperature
        atmospheric_transmission = amount of light absorbed
        FWHM5400 = size of PSF at 5400A
        *_dt = their time derivatives
        
        """
        self.sky_=sky
        self.Tccd_=Tccd
        self.atmo_=atmospheric_transmission
        self.psf_=FWHM5400

        self.dskydt=dsky_dt
        self.dTccddt=dTccd_dt
        self.datmodt=datmospheric_transmission_dt
        self.dpsfdt=dFWHM5400_dt

    def Tccd(self,t):
        return self.Tccd_+self.dTccddt*t

    def sky(self,t):
        return self.sky_+self.dskydt*t

    def atmospheric_transmission(self, t):
        return self.atmo_+self.datmodt*t

    def FWHM5400(self, t):
        return self.psf_+self.dpsfdt*t
        

        
