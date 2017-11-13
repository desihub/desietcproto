"""High-level calculator for online exposure-time forecasts.
"""
from __future__ import print_function, division

import numpy as np

import scipy.interpolate

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

    Use :meth:`will_timeout`, :meth:`get_snr_now` and :meth:`get_remaining` to
    query the latest models.

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
    seed : int or None
        Random number seed to use for reproducible random sampling of the
        signal and background rate models.
    """
    def __init__(self, alpha, dalpha, beta, dbeta,
                 sig0, dsig0, tcorr_sig, bg0, dbg0, tcorr_bg, t0, snr_goal,
                 dtmax=4000., npred=401, seed=None):
        self.t0 = t0
        assert snr_goal > 0
        self.snr_goal = snr_goal
        # Remember calibration constants.
        assert alpha > 0 and dalpha >= 0, 'Invalid alpha, dalpha'
        assert beta > 0 and dbeta >= 0, 'Invalid beta, dbeta'
        self.alpha = alpha
        self.beta = beta
        self.dalpha = dalpha
        self.dbeta = dbeta
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
        # Initialize random numbers.
        self.gen = np.random.RandomState(seed)
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

    def will_timeout(self):
        """Is this exposure expected to time out before reaching its goal SNR?

        Based on all signal and background rate updates recorded so far.

        Returns
        -------
        bool
            True if the current exposure is not expected to achieve ``snr_goal``
            with an exposure duration less than ``dtmax``.
        """
        return self.dt_goal >= self.dt_pred[-1]

    def get_remaining(self, tnow):
        """Forecast the remaining exposure time in seconds.

        Based on all signal and background rate updates recorded so far.

        Parameters
        ----------
        tnow : float
            Current time used for forecasting.  Must be >= t0.

        Returns
        -------
        float
            Remaining time in seconds, which might be < 0 if the goal SNR
            has already been achieved.
        """
        assert tnow >= self.t0, 'Expected tnow >= t0'
        return self.dt_goal - (tnow - self.t0)

    def get_snr_now(self, tnow, CL=0.6827, nsamples=1000):
        """Estimate the current integrated SNR.

        Based on all signal and background rate updates recorded so far.

        This calculation is relatively expensive since it requires generating
        random samples from the signal and background rate models.

        Parameters
        ----------
        tnow : float
            Current time used for forecasting.  Must be between t0, t0+dtmax.
        CL : float
            Confidence level to use for the returned range.  Must be 0-1.
        nsamples : int
            Number of random samples of the signal and background rate models
            to generate. A larger number gives more accurate ranges but takes
            longer to run. Must be > 0.

        Returns
        -------
        tuple
            Tuple (lo, hi) containing the true integrated SNR with posterior
            probability CL.
        """
        assert tnow >= self.t0, 'Expected tnow >= t0'
        assert tnow <= self.t0 + self.dt_pred[-1], 'Expected tnow <= t0 + dtmax'
        assert 0 < CL < 1, 'Invalid CL'
        assert nsamples > 0, 'Expected nsamples > 0'
        # Generate relative timestamps covering [0, tnow - t0].
        last = np.argmax(self.t0 + self.dt_pred >= tnow)
        assert self.t0 + self.dt_pred[last] >= tnow
        assert last == 0 or (self.t0 + self.dt_pred[last - 1] < tnow)
        dt = self.dt_pred[:last + 1].copy()
        dt[last] = tnow - self.t0
        # Generate random realizations of SNR at tnow.
        _, _, snr_now = self.get_samples(dt, nsamples)
        assert snr_now.shape == (nsamples, last + 1)
        # Estimate percentiles that cover a central fraction CL.
        lo = 50 * (1.0 - CL)
        cuts = (lo, 100 - lo)
        assert np.allclose(cuts[1] - cuts[0], 100 * CL)
        return np.percentile(snr_now[:, -1], cuts)

    def _update_model(self, dt, rate, drate, rate0, drate0, tcorr):
        """Internal method to update a rate model.

        This method is used by both :meth:`update_signal` and
        :meth:`update_background` to update their respective models.

        Call :meth:`update_snr` after calling this method to calculate the
        updated SNR.

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
            ``self.dt_pred``, and ``gp`` is the updated Gaussian process model.
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

    def _eval_snr(self, t, S, B):
        """Evaluate the SNR for calibrated rates S and B at increasing times t.

        Automatically broadcasts over input arrays whose last index corresponds
        to time.

        Parameters
        ----------
        t : array
            1D array of N increasing times in seconds.
        S : array
            Array with shape (...,N) of calibrated signal rates in counts per
            second at each time.
        B : array
            Array with shape (...,N) of calibrated background rates  in counts
            per second at each time.

        Returns
        -------
        array
            Array with shape (...,N) of cummulative signal-to-noise ratios at
            each time calculated as :math:`S / np.sqrt(S + B)`.  Any time when
            :math:`S+B \le 0` will have SNR = 0.
        """
        assert len(t.shape) == 1, 't must be 1D array'
        # Use trapezoid rule to integrate cummulatively over each interval.
        tstep = np.diff(t)
        assert np.all(tstep > 0)
        assert S.shape == B.shape, 'S, B must have same shape'
        Scum = np.cumsum(0.5 * tstep * (S[...,:-1] + S[...,1:]), axis=-1)
        Bcum = np.cumsum(0.5 * tstep * (B[...,:-1] + B[...,1:]), axis=-1)
        snr = np.zeros_like(S)
        nonzero = Scum + Bcum > 0
        snr[...,1:][...,nonzero] = Scum[nonzero] / np.sqrt(
            Scum[nonzero] + Bcum[nonzero])
        return snr

    def get_samples(self, dt, nsamples):
        """Generate random samples of our signal and background rate models.

        Parameters
        ----------
        dt : array of floats
            Array of times where models should be sampled. Times are in
            seconds relative to the exposure start. Must all be
            between 0 and dtmax.
        nsamples : int
            Number of independent samples to generate. Must be > 0.

        Returns
        -------
        tuple
            Tuple (S, B, snr) of calibrated signal and background rate samples,
            and the corresponding SNR values. Each is an array of floats with
            dimensions (nsamples, len(dt)).
        """
        dt = np.asarray(dt)
        assert (np.all(dt >= 0) and
                np.all(dt <= self.dt_pred[-1])), 'dt not in range [0, dtmax]'
        assert nsamples > 0, 'Expected nsamples > 0'
        # Generate random realizations of uncalibrated signal, bg rates
        # from the latest Gaussian process rate models.
        X = dt.reshape(-1, 1)
        S_samples = (
            self.sig0 + self.sig_model.sample_y(X, nsamples)).T
        B_samples = (
            self.bg0 + self.bg_model.sample_y(X, nsamples)).T

        # Apply calibration with random errors.
        if self.dalpha > 0:
            S_samples *= self.gen.normal(
                loc=self.alpha, scale=self.dalpha, size=(nsamples, 1))
        else:
            S_samples *= self.alpha
        if self.dbeta > 0:
            B_samples *= self.gen.normal(
                loc=self.beta, scale=self.dbeta, size=(nsamples, 1))
        else:
            B_samples *= self.beta

        # Calcuate SNR for each (S,B) pair.
        snr_samples = self._eval_snr(dt, S_samples, B_samples)

        return S_samples, B_samples, snr_samples

    def _update_snr(self):
        """Internal method to update SNR forecast.

        Uses the the most recent updates to the signal and background models
        to forecast the calibrated SNR out to ``dtmax``.  This forecast is
        then used to estimate the current SNR and the time remaining until
        ``snr_goal`` is achieved.

        Sets the flag ``attr:timeout`` if the updated model indicates that
        this exposure will not complete before ``dtmax`` is reached.
        """
        # Calculate nominal calibrated signal and background rates.
        S = self.alpha * np.asarray(self.sig_pred)
        B = self.beta * np.asarray(self.bg_pred)

        # Evaluate the corresponding nominal SNR model.
        snr = self._eval_snr(self.dt_pred, S, B)
        assert np.all(np.diff(snr) >= 0), 'nominal SNR is not increasing'

        # Estimate when the nominal SNR model hits the SNR goal using
        # linear interpolation.
        # If the result equals self.dt_pred[-1], this indicates that
        # the exposure is not expected to reach its SNR goal within the
        # maximum allowed exposure time.
        interpolator = scipy.interpolate.interp1d(
            snr, self.dt_pred, kind='cubic', assume_sorted=True,
            bounds_error=False, fill_value=self.dt_pred[-1])
        self.dt_goal = interpolator(self.snr_goal)
