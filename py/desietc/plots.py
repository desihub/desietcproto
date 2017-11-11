"""Plot utilities for the exposure-time calculator package
"""
from __future__ import print_function, division

import numpy as np


def plot_calculator(calc, tnow, nsigma=1.0, nsamples=0, save=None):
    """Plot the current state of a calculator.

    The matplotlib and basemap packages must be installed to use this function.

    Parameters
    ----------
    nsamples : int
        Number of random samples to plot on each graph.  No samples are
        plotted when zero. Must be >= 0.
    save : string or None
        Name of file where plot should be saved.  Format is inferred from
        the extension.

   Returns
    -------
    tuple
        Tuple (figure, axes) returned by ``plt.subplots()``.
    """
    import matplotlib.pyplot as plt

    assert nsamples >= 0, 'Expected nsamples >= 0'

    dt_now = tnow - calc.t0
    dt = calc.dt_pred

    # Get the estimated time remaining.
    remaining = calc.get_remaining(tnow)
    dt_goal = dt_now + remaining

    # Get the 68% CL estimated range of SNR now.
    snr_lo, snr_hi = calc.get_snr_now(tnow)

    # Calculate the mean predicted S, B, SNR.
    S = calc.alpha * calc.sig_pred
    B = calc.beta * calc.bg_pred
    snr = calc._eval_snr(dt, S, B)

    # Calculate the predicted errors with and w/o calibration uncertainties.
    dS = calc.alpha * calc.dsig_pred
    dB = calc.beta * calc.dbg_pred
    dSc = np.sqrt((calc.dalpha / calc.alpha * S) ** 2 + dS ** 2)
    dBc = np.sqrt((calc.dbeta / calc.beta * B) ** 2 + dB ** 2)

    # Prepare arrays for filled plots.
    dt_fill = np.concatenate([dt, dt[::-1]])
    S_fill = np.concatenate([S - nsigma * dS, (S + nsigma * dS)[::-1]])
    B_fill = np.concatenate([B - nsigma * dB, (B + nsigma * dB)[::-1]])
    Sc_fill = np.concatenate([S - nsigma * dSc, (S + nsigma * dSc)[::-1]])
    Bc_fill = np.concatenate([B - nsigma * dBc, (B + nsigma * dBc)[::-1]])

    # Initialize the plot to fill an 8.5" x 11" page.
    fig, ax = plt.subplots(3, 1, figsize=(8.5, 11), sharex=True)

    # Plot signal rate data and model.
    for S_sample in calc.S_samples[:nsamples]:
        ax[0].plot(dt, S_sample, 'g--', lw=0.5)
    ax[0].fill(dt_fill, Sc_fill, alpha=.15, fc='g', ec='None')
    ax[0].fill(dt_fill, S_fill, alpha=.25, fc='g', ec='None')
    ax[0].plot(dt, S, 'g-', alpha=0.5, lw=2)
    if calc.dtsig:
        ax[0].errorbar(calc.dtsig, calc.alpha * np.asarray(calc.sig),
                       calc.alpha * np.asarray(calc.dsig),
                       fmt='k.', ms=10, lw=1)
    ax[0].set_ylabel('Signal Rate $S$')
    ax[0].axvline(dt_goal, ls='-', c='r')
    ax[0].axvline(dt_now, ls='--', c='r')
    ax[0].grid()

    # Plot background rate data and model.
    for B_sample in calc.B_samples[:nsamples]:
        ax[1].plot(dt, B_sample, 'b--', lw=0.5)
    ax[1].fill(dt_fill, Bc_fill, alpha=.15, fc='b', ec='None')
    ax[1].fill(dt_fill, B_fill, alpha=.25, fc='b', ec='None')
    ax[1].plot(dt, B, 'b-')
    if calc.dtbg:
        ax[1].errorbar(calc.dtbg, calc.beta * np.asarray(calc.bg),
                       calc.beta * np.asarray(calc.dbg),
                       fmt='k.', ms=10, lw=1)
    ax[1].set_ylabel('Background Rate $B$')
    ax[1].axvline(dt_goal, ls='-', c='r')
    ax[1].axvline(dt_now, ls='--', c='r')
    ax[1].grid()

    # Plot SNR predictions.
    for snr_sample in calc.snr_samples[:nsamples]:
        ax[2].plot(dt, snr_sample, 'k--', lw=0.5)
    ax[2].plot(dt, snr, '-', c='gray', lw=2)
    ax[2].errorbar(dt_now, 0.5 * (snr_lo + snr_hi),
                   yerr=0.5 * (snr_hi - snr_lo), fmt='k', lw=6, ms=0,
                   capsize=8, capthick=1)
    ax[2].axvline(dt_goal, ls='-', c='r')
    ax[2].axvline(dt_now, ls='--', c='r')
    ax[2].axhline(calc.snr_goal, ls=':', c='r')
    ax[2].set_ylim(0, 1.2 * calc.snr_goal)
    ax[2].set_xlim(dt[0], dt[-1])
    ax[2].set_xlabel('Exposure Duration [s]')
    ax[2].set_ylabel('Signal-to-Noise Ratio $S/\sqrt{S+B}$')
    ax[2].grid()

    # Add labels.
    def tfmt(s):
        m = int(np.floor(s / 60.))
        s -= 60 * m
        return '{:d}:{:04.1f}'.format(m, s)

    label = 'Elapsed: {}'.format(tfmt(dt_now))
    xy = 0.97, 0.20
    ax[2].annotate(label, xy, xy, 'axes fraction', 'axes fraction',
                   horizontalalignment='right', fontsize=24)
    label = 'Remaining: {}'.format(tfmt(remaining))
    xy = 0.97, 0.05
    ax[2].annotate(label, xy, xy, 'axes fraction', 'axes fraction',
                   horizontalalignment='right', fontsize=24)

    plt.subplots_adjust(hspace=0.05, left=0.05, right=)
    if save:
        plt.savefig(save)
    return fig, ax
