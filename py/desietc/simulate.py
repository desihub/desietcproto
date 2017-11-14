"""Plot utilities for the exposure-time calculator package
"""
from __future__ import print_function, division

import numpy as np

import desietc.calculator


def simulate_measurements(dt_pred, interval, rel_accuracy,
                          initial=1., slope=0., gen=None):
    """Simulate a time series of signal or background rate measurements.

    Parameters
    ----------
    dt_pred : array of floats
        Increasing times in seconds when the true values are tabulated.
        Measurements are simulated up to dt_pred[-1].
    interval : float
        Time between measurements in seconds. The first measurement occurs
        after interval seconds, not at zero seconds.
    rel_accuracy : float
        Relative accuracy of each measurement. Gaussian noise with this
        fractional RMS will be simulated.
    initial : float
        Mean initial measurement value. Actual initial value will have
        measurement noise added.
    slope : float
        Rate of change per second of mean measurement value.
    gen : numpy.random.RandomState or None
        Use the specified random number generator for reproducible results.

    Returns
    -------
    tuple
        Tuple (dt, rate, error, truth) of simulated measurements rate +/- error
        at equally spaced times dt, with corresponding true (noise-free) values.
    """
    if gen is None:
        gen = np.random.RandomState()
    # Tabulate true values.
    truth = initial + slope * dt_pred
    # Tabulate measured values.
    n = int(np.floor(dt_pred[-1] / float(interval)))
    dt = interval * np.arange(1, n + 1)
    rate = initial + slope * dt
    error = rel_accuracy * rate
    rate += gen.normal(scale=error)
    return dt, rate, error, truth


def simulate_exposure(snr_goal, texp_goal, rho, s0err, b0err, alpha_err,
                      beta_err, sprior=1., bprior=1., priorerr=0.5,
                      estimate_snr=False, sig_period=120., bg_period=60.,
                      etc_period=60., gen=None):
    """Simulate the ETC during a single exposure.

    Assume that the true signal and background rates are constant during the
    exposure (for now).

    Parameters
    ----------
    snr_goal : float
        Goal SNR value to achieve.
    texp_goal : float
        Time when snr_goal should be achieved, in seconds.
    rho : float
        Ratio of the initial calibrated signal and background rates.
    s0err : float
        Fractional RMS error in the uncalibrated signal rate estimates.
        Averages down as more rate estimates are accumulated.
    b0err : float
        Fractional RMS error in the uncalibrated background rate estimates.
        Averages down as more rate estimates are accumulated.
    alpha_err : float
        Fractional RMS error of the signal calibration constant.
        Does not average down as rate estimates are accumulated.
    beta_err : float
        Fractional RMS error of the signal calibration constant.
        Does not average down as rate estimates are accumulated.
    sprior : float
        Ratio of prior on signal rate and actual initial rate. Must be > 0.
    bprior : float
        Ratio of prior on background rate and actual initial rate. Must be > 0.
    priorerr : float
        Fractional error on the signal and background rate priors. Must be > 0.
    estimate_snr : bool
        Estimate SNR after each simulation step when True. Simulations are
        faster but do not return snr_range values when False.
    sig_period : float
        Period for signal rate measurements in seconds.
    bg_period : float
        Period for background rate measurements in seconds.
    etc_period : float
        Period for exposure-time calculator forecast updates in seconds.
        Should be the lowest-common multiple of sig_period and bg_period
        to avoid any discretization effects.
    gen : numpy.random.RandomState or None
        Use the specified random number generator for reproducible
        random samples. Ignored unless estimate_snr is True.

    Returns
    -------
    tuple
        Tuple (calc, tnow, snr_actual, snr_true, telapsed, tremaining,
        snr_range).  The snr_range array is empty when called with
        estimate_snr False.
    """
    if gen is None:
        gen = np.random.RandomState()

    # Calculate the constant calibrated signal and background rates to
    # simulate in order to reach snr_goal at texp_goal given rho.
    B0 = snr_goal ** 2 / texp_goal * (1 + rho) / rho ** 2
    S0 = rho * B0

    # Fix calibration constants to 1 w/o loss of generality.
    alpha, beta = 1., 1.
    dalpha, dbeta = alpha_err * alpha, beta_err * beta

    # Calculate uncalibrated constant rates and errors.
    s0, b0 = S0 / alpha, B0 / beta
    ds0, db0 = s0err * s0, b0err * b0

    t0, dtmax = 1e6, 4000.

    # Use rate priors that differ from the true rates with 25% rms.
    s0prior = s0 * sprior
    b0prior = b0 * bprior
    # Assign 50% error to the priors.
    ds0prior = priorerr * s0prior
    db0prior = priorerr * b0prior
    # Fix correlation times of 500s, 1500s.
    stcorr, btcorr = 500., 1500.

    # Initialize an ETC for this exposure.
    calc = desietc.calculator.Calculator(
        alpha, dalpha, beta, dbeta,
        s0prior, ds0prior, stcorr, b0prior, db0prior, btcorr,
        t0=t0, snr_goal=snr_goal, dtmax=dtmax)
    assert not calc.will_timeout()

    # Generate signal rate updates.
    dtsig, sig, dsig, sigtrue = simulate_measurements(
        calc.dt_pred, sig_period, s0err, s0, 0., gen)

    # Generate background rate updates.
    dtbg, bg, dbg, bgtrue = simulate_measurements(
        calc.dt_pred, bg_period, b0err, b0, 0., gen)

    # Calculate the true snr evolution.
    snr_true = calc._eval_snr(calc.dt_pred, alpha * sigtrue, beta * bgtrue)

    # Record the initial forecast (based only the priors).
    telapsed = [0.]
    tremaining = [calc.get_remaining(t0)]
    if estimate_snr:
        snr_range = [calc.get_snr_now(t0, gen=gen)]
    else:
        snr_range = []

    # Loop over simulated time steps.
    tnow = t0
    isig, ibg = 0, 0
    while calc.get_remaining(tnow) > 0.01:
        # Update the simulation time.
        tnow += min(etc_period, calc.get_remaining(tnow))
        # Add any signal or background updates during this step.
        while t0 + dtbg[ibg] <= tnow:
            calc.update_background(t0 + dtbg[ibg], bg[ibg], dbg[ibg])
            ibg += 1
        while t0 + dtsig[isig] <= tnow:
            calc.update_signal(t0 + dtsig[isig], sig[isig], dsig[isig])
            isig += 1
        # Record the new forecast.
        telapsed.append(tnow - t0)
        tremaining.append(calc.get_remaining(tnow))
        if estimate_snr:
            snr_range.append(calc.get_snr_now(tnow, gen=gen))

    # Calculate the true SNR when the exposure ends.
    snr_actual = np.interp(tnow - t0, calc.dt_pred, snr_true)

    # Convert python arrays to numpy.
    telapsed = np.array(telapsed)
    tremaining = np.array(tremaining)
    snr_range = np.array(snr_range)

    return calc, tnow, snr_actual, snr_true, telapsed, tremaining, snr_range
