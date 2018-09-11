"""Signal and background calibration models.
"""
from __future__ import print_function, division

import numpy as np


class SampleHold(object):

    def __init__(self, num_initial=5):
        self.samples = []
        self.num_initial = num_initial

    def reset(self):
        self.samples = []

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
        self.flux_sh = SampleHold(num_initial=5)
        self.airmass_exponent = airmass_exponent
        self.dust_coef = dust_coef

    def reset(self):
        self.flux_sh.reset()

    def rate(self, fwhm, dfwhm, flux, dflux):
        """Uncalibrated signal rate.
        """
        f_seeing = self.a + self.b * fwhm + self.c * fwhm ** 2
        flux0 = self.flux_sh.add(flux)
        f_flux = 1 + self.rGFA * (flux - flux0) / flux0
        y = f_seeing * f_flux
        dy_dfwhm = (self.b + 2 * self.c * fwhm) * f_flux
        dy_dflux = self.rGFA / flux0 * f_seeing
        dy = np.sqrt((dy_dfwhm * dfwhm) ** 2 + (dy_dflux * dflux) ** 2)
        return y, dy

    def alpha(self, airmass, Ebv, alpha0=22.):
        """Signal rate calibration.
        """
        value = (alpha0 *
            airmass ** self.airmass_exponent *
            10 ** (-2 * self.dust_coef * Ebv / 2.5))
        return value, 0.2 * value


class BackgroundCalib(object):

    def __init__(self, rSC=1.5, Bread=40000.):
        self.rSC = rSC
        self.Bread = Bread
        self.sky_sh = SampleHold(num_initial=1)

    def reset(self):
        self.sky_sh.reset()

    def rate(self, sky, dsky):
        """Uncalibrated sky background rate.
        """
        sky0 = self.sky_sh.add(sky)
        y = 1 + self.rSC * (sky - sky0) / sky0
        dy = self.rSC / sky0 * dsky
        return y, dy

    def beta(self, airmass, beta0=2400.):
        """Sky background rate calibration.

        TODO: add calculation of moon factor.
        """
        value = beta0 * airmass
        return value, 0.2 * value
