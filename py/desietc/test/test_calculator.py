from __future__ import print_function, division

import unittest

from ..calculator import *


class TestConfig(unittest.TestCase):

    def test_basic(self):
        """Test that we can create a calculator and add signal, bg"""
        calc = Calculator(
            1., 0.1, 1., 0.1,
            1.0, 0.5, 2000.,
            1.0, 0.5, 2000.,
            0., 10.)
        calc.update_signal(1000., 1., 0.1)
        calc.update_background(900., 0.5, 0.2)

    def test_initial(self):
        """Verify initial forecast, before any rate updates"""
        snr_goal = 10.
        dsig0 = 0.1
        dbg0 = 0.2
        t0, tcorr = 1e6, 1e3
        for alpha in (0.5, 2.0):
            for beta in (0.5, 2.0):
                for sig0 in (0.9, 1.1):
                    for bg0 in (0.9, 1.1):
                        calc = Calculator(
                            alpha, 0, beta, 0, sig0, dsig0, tcorr,
                            bg0, dbg0, tcorr, t0, snr_goal, seed=123)
                        assert not calc.will_timeout()
                        tpred = snr_goal ** 2 * (
                            alpha * sig0 + beta * bg0) / (alpha * sig0) ** 2
                        assert np.allclose(tpred, calc.get_remaining(t0))
