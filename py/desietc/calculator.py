"""High-level calculator for online exposure-time forecasts.
"""
from __future__ import print_function, division

import numpy as np

import sklearn.gaussian_process


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
    tcorr_sig : float
        Correlation time for changes in signal rate. Acts as a prior on
        how rapidly the signal rate changes. Must be > 0.
    bg0 : float
        Prior on raw (uncalibrated) background rate. Must be > 0.
    dbg0 : float
        One sigma error on ``bg0``.  Must be > 0.
    tcorr_bg : float
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
                 sig0, dsig0, tcorr_sig, bg0, dbg0, tcorr_bg,
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
        assert tcorr_sig > 0 and tcorr_bg > 0, 'Invalid tcorr_sig, tcorr_bg'
        self.sig0 = sig0
        self.dsig0 = dsig0
        self.bg0 = bg0
        self.dbg0 = dbg0
        self.tcorr_sig = tcorr_sig
        self.tcorr_bg = tcorr_bg
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
            [], [], [], self.sig0, self.dsig0, self.tcorr_sig)
        self.bg_pred, self.dbg_pred, self.bg_model = self._update_model(
            [], [], [], self.bg0, self.dbg0, self.tcorr_bg)
        self._update_snr()

    def update_signal(self, tsig, sig, dsig):
        """Update the signal estimate.

        Can be called with timestamps ``tsig`` in any order.

        Parameters
        ----------
        tsig : float
            Timestamp in seconds for this signal rate estimate, using the same
            (arbitrary) origin as the ``t0`` passed to the constructor.
            The elapsed time ``tsig - t0`` must be >= 0.
        sig : float
            Estimated raw (uncalibrated) signal rate at ``tsig``. Must be > 0.
        dsig : float
            Estimated 1-sigma error on the value of ``sig``. Must be > 0.
        """
        assert tsig >= self.t0, 'Expected tsig >= t0'
        assert sig > 0 and dsig > 0, 'Invalid sig, dsig'
        self.sig.append(sig)
        self.dsig.append(dsig)
        self.dtsig.append(tsig - self.t0)
        # Update signal rate model.
        self.sig_pred, self.dsig_pred, self.sig_model = self._update_model(
            self.dtsig, self.sig, self.dsig, self.sig0, self.dsig0,
            self.tcorr_sig)
        # Update estimates of calibrated SNR.
        self._update_snr()

    def update_background(self, tbg, bg, dbg):
        """Update the background estimate.

        Can be called with timestamps ``tbg`` in any order.

        Parameters
        ----------
        tbg : float
            Timestamp in seconds for this background rate estimate, using the
            same (arbitrary) origin as the ``t0`` passed to the constructor.
            The elapsed time ``tbg - t0`` must be >= 0.
        bg : float
            Estimated raw (uncalibrated) background rate at ``tbg``.
            Must be > 0.
        dbg : float
            Estimated 1-sigma error on the value of ``bg``. Must be > 0.
        """
        assert tbg >= self.t0, 'Expected tbg >= t0'
        assert bg > 0 and dbg > 0, 'Invalid bg, dbg'
        self.bg.append(bg)
        self.dbg.append(dbg)
        self.dtbg.append(tbg - self.t0)
        # Update background rate model.
        self.bg_pred, self.dbg_pred, self.bg_model = self._update_model(
            self.dtbg, self.bg, self.dbg, self.bg0, self.dbg0, self.tcorr_bg)
        # Update estimates of calibrated SNR.
        self._update_snr()

    def _update_model(self, dt, rate, drate, rate0, drate0, tcorr):
        """Internal method to update a rate model.

        This method is used by both :meth:`update_signal` and
        :meth:`update_background` to update their respective models.

        Call :meth:`update_snr` after calling this method to calculate the
        update SNR.

        Parameters
        ----------
        dt : array
            1D array of N elapsed times in seconds since ``self.t0``, not
            necessarily in increasing order. An empty array (N=0) is allowed.
        rate : array
            1D array of N raw (uncalibrated) rate estimates corresponding to
            each ``dt`` value. Must be > 0.
        drate : array
            1D array of N 1-sigma errors on each ``rate`` value.  Must be > 0.
        rate0 : float
            Prior on the raw (uncalibrated) rate that is combined with the
            rate estimates to obtain a posterior estimate of the rate. Must
            be > 0.
        drate0 : float
            1-sigma error on the prior value ``rate0``. Must be > 0.
        tcorr : float
            Correlation time in seconds for changes in the raw rate. Acts as a
            prior on how rapidly the rate changes. Must be > 0.

        Returns
        -------
        tuple
            Tuple (pred, dpred, gp) where pred and dpred are arrays of predicted
            rates and their 1-sigma errors for each elapsed time in
            ``self.dt_steps``, and ``gp`` is the updated Gaussian process model.
        """
        dt = np.asarray(dt)
        rate = np.asarray(rate)
        drate = np.asarray(drate)
        assert dt.shape == rate.shape == drate.shape
        # Fit a Gaussian process regression model.
        kernel = (
            sklearn.gaussian_process.kernels.ConstantKernel(drate0 ** 2) *
            sklearn.gaussian_process.kernels.RBF(tcorr))
        gp = sklearn.gaussian_process.GaussianProcessRegressor(
            # See https://github.com/scikit-learn/scikit-learn/pull/10026
            kernel=kernel, optimizer=None, alpha=drate ** 2)
        if len(dt):
            gp.fit(dt.reshape(-1, 1), rate - rate0)
        # Evaluate the model.
        pred, dpred = gp.predict(self.dt_pred.reshape(-1, 1), return_std=True)
        pred += rate0
        return pred, dpred, gp

    def _update_snr(self):
        """Internal method to update values of snr_now and t_remain.
        """
        pass
