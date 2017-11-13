"""Plot utilities for the exposure-time calculator package
"""
from __future__ import print_function, division

import numpy as np


def simulate_measurements(dt_pred, interval=60, rel_accuracy=0.04,
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
