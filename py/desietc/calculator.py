"""High-level calculator for online exposure-time forecasts.
"""
from __future__ import print_function, division

import numpy as np


class Calculator(object):
    """Create a new calculator object for a single spectrograph exposure.

    Model the time-evolution of independent signal and background rate
    estimators, and use models to forecast the currently acheived signal-to-
    noise ratio (SNR) and the remaining exposure time.

    Register updates to the estimated signal and background rates using
    :meth:`update_signal` and :meth:`update_background`, respectively. The
    model is initialized with an initial guess at these rates which is
    gradually replaced by actual rate updates.  The signal and background
    rate estimates are assumed to be statistically independent.

    The target SNR is calculated as :math:`S/\sqrt(S+B)` where
    :math:`S = \\alpha s` and :math:`B = \\beta b` are calibrated versions of
    the estimates :math:`s` and :math:`b` passed to :meth:`update_signal` and
    :meth:`update_background`.  The dimensions of :math:`a`, :math:`b`,
    :math:`\\alpha`, and :math:`\\beta` are arbitrary, but the combinations
    :math:`S = \\alpha s` and :math:`B = \\beta b` must be dimensionless and
    scaled appropriately to have counting statistics.

    The calibration constants :math:`\\alpha`, and :math:`\\beta` can have
    associated errors, :math:`\\sigma_{\\alpha}` and :math:`\\sigma_{\\beta}`
    that are propagated to SNR estimates. In general, the calibration depends
    on the program (DARK, GRAY, BRIGHT) and can be refined by comparing this
    object's predictions with the pipeline outputs from previous
    spetrograph exposures.

    Parameters
    ----------
    alpha : float
        Constant conversion factor from raw signal rate to fiducial target
        signal rate. Must be > 0.
    dalpha : float
        Error on alpha. Must be > 0.
    beta : float
        Constant conversion factor from raw signal rate to fiducial target
        signal rate. Must be > 0.
    dbeta : float
        Error on beta. Must be > 0.
    sig0 : float
        Prior on raw (uncalibrated) signal rate. Must be > 0.
    dsig0 : float
        One sigma error on ``sig0``.  Must be > 0.
    tsig0 : float
        Correlation time for changes in signal rate. Acts as a prior on
        how rapidly the signal rate changes. Must be > 0.
    bg0 : float
        Prior on raw (uncalibrated) background rate. Must be > 0.
    dbg0 : float
        One sigma error on ``bg0``.  Must be > 0.
    tsig0 : float
        Correlation time for changes in background rate. Acts as a prior on
        how rapidly the background rate changes. Must be > 0.
    t0 : float
        Timestamp for when shutter was opened for the current exposure,
        in units of seconds with an arbitrary origin.
    snr_goal : float
        Value of calibrated SNR that we ideally want to achieve. Must be > 0.
    dtmax : float
        Maximum allowed total exposure time in seconds.  Must be > 0.
    npred : int
        Number of equally spaced times spanning [0, dtmax] where predictions
        are calculated internally.
    """
    def __init__(self, alpha, dalpha, beta, dbeta,
                 sig0, dsig0, tsig0, bg0, dbg0, tbg0,
                 t0, snr_goal, dtmax=5000., npred=1000):
        # Remember the exposure start time and SNR goal.
        self.t0 = t0
        assert snr_goal > 0
        self.snr_goal = snr_goal
        # Remember calibration constants.
        assert alpha > 0 and dalpha >= 0, 'Invalid alpha, dalpha'
        assert beta > 0 and dbeta >= 0, 'Invalid beta, dbeta'
        self.alpha = alpha
        self.beta = beta
        self.dalpha2frac = (dalpha / alpha) ** 2
        self.dbeta2frac = (dbeta / beta) ** 2
        # Remember priors.
        assert sig0 > 0 and dsig0 > 0, 'Invalid sig0, dsig0'
        assert bg0 > 0 and dbg0 > 0, 'Invalid bg0, dbg0'
        assert tsig0 > 0 and tbg0 > 0, 'Invalid tsig0, tbg0'
        self.sig0 = sig0
        self.dsig0 = dsig0
        self.bg0 = bg0
        self.dbg0 = dbg0
        self.tsig0 = tsig0
        self.tbg0 = tbg0
        # Initialize uncalibrated signal and background rate estimates.
        self.sig = []
        self.dsig = []
        self.dtsig = []
        self.bg = []
        self.dbg = []
        self.dtbg = []
        # Initialize times relative to t0 where model is predicted.
        self.dt_pred = np.linspace(0., dtmax, npred)
        # Initialize models and our prediction.
        self.sig_pred, self.dsig_pred, self.sig_model = self._update_model(
            [], [], [], self.sig0, self.dsig0, self.tsig0)
        self.bg_pred, self.dbg_pred, self.bg_model = self._update_model(
            [], [], [], self.bg0, self.dbg0, self.tbg0)
        self._update_snr()

    def _update_model(self, dt, rate, drate, rate0, drate0, smoothness):
        """Update rate model"""
        return None, None, None

    def _update_snr(self):
        """Update values of snr_now and t_remain using _update_model() results.
        """
        pass
