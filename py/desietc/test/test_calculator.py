from __future__ import print_function, division

import unittest

from ..calculator import *


class TestCalculator(unittest.TestCase):

    def test_basic(self):
        """Test that we can create a calculator and add signal, bg"""
        calc = Calculator(
            1., 0.1, 1., 0.1,
            1.0, 0.5, 2000.,
            1.0, 0.5, 2000.,
            0., 10.)
        calc.update_signal(1000., 1., 0.1)
        calc.update_background(900., 0.5, 0.2)

    def test_initial_forecast(self):
        """Verify initial forecast"""
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

    def test_initial_samples(self):
        """Verify mean and RMS of initial samples"""
        alpha, beta = 0.5, 1.2
        sig0, bg0 = 0.9, 1.1
        dsig0, dbg0 = 0.10, 0.15
        t0, snr_goal = 1e6, 10.
        dbeta = 0.15
        dtmax = 4000.
        for dalpha in (0., 0.10):
            for tcorr in (1e-2 * dtmax, 1e2 * dtmax):
                calc = Calculator(
                    alpha, dalpha, beta, dbeta, sig0, dsig0, tcorr,
                    bg0, dbg0, tcorr, t0, snr_goal, seed=123)
                # Predict mean and stddev of initial samples.
                Spred = alpha * sig0
                dSpred = Spred * np.sqrt(
                    (dalpha / alpha) ** 2 + (dsig0 / sig0) ** 2)
                Bpred = beta * bg0
                dBpred = Bpred * np.sqrt(
                    (dbeta / beta) ** 2 + (dbg0 / bg0) ** 2)
                # Compare with numerical results at each time step.
                S_samples, B_samples, _ = calc.get_samples(
                    calc.dt_pred, nsamples=10000)
                Smean = np.mean(S_samples, axis=0)
                Sstd = np.std(S_samples, axis=0)
                Bmean = np.mean(B_samples, axis=0)
                Bstd = np.std(B_samples, axis=0)
                # Check the level of agreement.
                assert np.allclose(Smean, Spred, rtol=0.01)
                assert np.allclose(Bmean, Bpred, rtol=0.01)
                assert np.allclose(dSpred, Sstd, rtol=0.03)
                assert np.allclose(dBpred, Bstd, rtol=0.03)

    def test_short_correlation_time(self):
        alpha, beta = 0.5, 1.2
        dalpha, dbeta = 0.1, 0.05
        sig0, bg0 = 0.9, 1.1
        dsig0, dbg0 = 0.25, 0.25
        t0, snr_goal = 1e6, 10.
        dtmax = 4000.
        # Calculate predictions for S, B before any updates.
        Spred0 = alpha * sig0
        dSpred0 = Spred0 * np.sqrt((dalpha / alpha) ** 2 + (dsig0 / sig0) ** 2)
        Bpred0 = beta * bg0
        dBpred0 = Bpred0 * np.sqrt((dbeta / beta) ** 2 + (dbg0 / bg0) ** 2)
        for tcorr in (1e-2 * dtmax, 1e2 * dtmax):
            calc = Calculator(
                alpha, dalpha, beta, dbeta, sig0, dsig0, tcorr,
                bg0, dbg0, tcorr, t0, snr_goal, seed=123)
            # Update signal and background rate estimates.
            idx = len(calc.dt_pred) // 2
            tupdate = t0 + calc.dt_pred[idx]
            sig, dsig = 2.0, 0.02
            bg, dbg = 0.5, 0.01
            # Calculate the predictions for S, B after the update.
            Spred = alpha * sig
            dSpred = Spred * np.sqrt((dalpha / alpha) ** 2 + (dsig / sig) ** 2)
            Bpred = beta * bg
            dBpred = Bpred * np.sqrt((dbeta / beta) ** 2 + (dbg / bg) ** 2)
            calc.update_signal(tupdate, sig, dsig)
            calc.update_background(tupdate, bg, dbg)
            # Calulate the predicted mean, std as a function of time.
            S_samples, B_samples, _ = calc.get_samples(
                calc.dt_pred, nsamples=10000)
            Smean = np.mean(S_samples, axis=0)
            Sstd = np.std(S_samples, axis=0)
            Bmean = np.mean(B_samples, axis=0)
            Bstd = np.std(B_samples, axis=0)
            # Verify that predictions match the update at the update time.
            assert np.allclose(Smean[idx], Spred, rtol=0.01)
            assert np.allclose(Sstd[idx], dSpred, rtol=0.03)
            assert np.allclose(Bmean[idx], Bpred, rtol=0.01)
            assert np.allclose(Bstd[idx], dBpred, rtol=0.03)
            for idx0 in (0, len(calc.dt_pred) - 1):
                if tcorr < dtmax:
                    # Verify that predictions match the prior far from the
                    # update time when using a short correlation time.
                    assert np.allclose(Smean[idx0], Spred0, rtol=0.01)
                    assert np.allclose(Sstd[idx0], dSpred0, rtol=0.03)
                    assert np.allclose(Bmean[idx0], Bpred0, rtol=0.01)
                    assert np.allclose(Bstd[idx0], dBpred0, rtol=0.03)
                else:
                    # Verify that the predictions match the update far from the
                    # update time when using a long correlation time.
                    assert np.allclose(Smean[idx0], Spred, rtol=0.01)
                    assert np.allclose(Sstd[idx0], dSpred, rtol=0.03)
                    assert np.allclose(Bmean[idx0], Bpred, rtol=0.01)
                    assert np.allclose(Bstd[idx0], dBpred, rtol=0.03)
