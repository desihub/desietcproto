from __future__ import print_function, division

import unittest

from ..simulate import *


class TestSimulate(unittest.TestCase):

    def test_reproducible(self):
        """Verify that random numbers are reproducible."""
        dt_pred = np.linspace(100., 200., 50)
        interval = 12.
        rel_accuracy = 0.1
        initial = 1.5
        slope = 1e-3
        gen1 = np.random.RandomState(123)
        dt1, rate1, error1, truth1 = simulate_measurements(
            dt_pred, interval, rel_accuracy, initial, slope, gen1)
        gen2 = np.random.RandomState(123)
        dt2, rate2, error2, truth2 = simulate_measurements(
            dt_pred, interval, rel_accuracy, initial, slope, gen2)
        assert np.array_equal(dt1, dt2)
        assert np.array_equal(rate1, rate2)
        assert np.array_equal(error1, error2)
        assert np.array_equal(truth1, truth2)

    def test_different(self):
        """Verify that different seeds give different measurements."""
        dt_pred = np.linspace(100., 200., 50)
        interval = 12.
        rel_accuracy = 0.1
        initial = 1.5
        slope = 1e-3
        gen1 = np.random.RandomState(123)
        dt1, rate1, error1, truth1 = simulate_measurements(
            dt_pred, interval, rel_accuracy, initial, slope, gen1)
        gen2 = None
        dt2, rate2, error2, truth2 = simulate_measurements(
            dt_pred, interval, rel_accuracy, initial, slope, gen2)
        assert np.array_equal(dt1, dt2)
        assert not np.array_equal(rate1, rate2)
        assert np.array_equal(error1, error2)
        assert np.array_equal(truth1, truth2)

    def test_simulate_exposure_basic(self):
        """Verify basic operation of an exposure simulation."""
        gen = np.random.RandomState(123)
        (calc, tnow, snr_actual, snr_true,
         telapsed, tremaining, snr_range) = simulate_exposure(
            snr_goal=7, texp_goal=1500, rho=1, s0err=0.1, b0err=0.05,
            alpha_err=0.1, beta_err=0.1, gen=gen)
