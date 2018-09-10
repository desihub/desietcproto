"""Signal and background calibration models.
"""
from __future__ import print_function, division

import numpy as np


class SampleHold(object):

    def __init__(self, num_initial=5):
        self.samples = []
        self.num_initial = num_initial

    def add(self, value):
        if len(self.samples) < self.num_initial:
            self.samples.append(value)
            if len(self.samples) == self.num_initial:
                self.hold = np.mean(self.samples)
            return value
        else:
            return self.hold


class SignalCalib(object):

    def __init__(self, a=2.04760, b=-1.18590, c=0.21231, rGFA=1.0,
                 airmass_exponent=-1.25, dust_coef=3.303):
        self.a = a
        self.b = b
        self.c = c
        self.rGFA = rGFA
        self.flux_sh = SampleHold()
        self.airmass_exponent = airmass_exponent
        self.dust_coef = dust_coef

    def rate(self, seeing_fwhm_arcsec, flux):
        """Uncalibrated signal rate.
        """
        f_seeing = self.a + self.b * seeing_fwhm_arcsec + self.c * seeing_fwhm_arcsec ** 2
        flux0 = self.flux_sh.add(flux)
        f_flux = 1 + self.rGFA * (flux - flux0) / flux0
        return f_seeing * f_flux

    def alpha(self, airmass, Ebv, alpha0=6.9):
        """Signal rate calibration.
        """
        value = (alpha0 *
            airmass ** self.airmass_exponent *
            10 ** (-2 * self.dust_coef * Ebv / 2.5))
        return value, 0.1 * value


class BackgroundCalib(object):

    def __init__(self, rSC=1.5):
        self.rSC = rSC
        self.flux_sh = SampleHold()

    def rate(self, flux):
        """Uncalibrated sky background rate.
        """
        flux0 = self.flux_sh.add(flux)
        return 1 + self.rSC * (flux - flux0) / flux0

    def beta(self, beta0=2400.):
        """Sky background rate calibration.

        TODO: add calculation of moon factor.
        """
        value = beta0
        return value, 0.1 * value
