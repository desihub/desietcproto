from __future__ import print_function, division

import unittest

from ..calib import *


class TestCalib(unittest.TestCase):

    def test_sample_hold(self):
        sh = SampleHold()
        xvec = np.arange(10)
        yvec = []
        for x in xvec:
            yvec.append(sh.add(x))
        assert np.array_equal(yvec, [0, 1, 2, 3, 4, 2, 2, 2, 2, 2])

    def test_signal_calib(self):
        scalib = SignalCalib()
        assert np.allclose(scalib.rate(1.1, 0.5), 1.)
        alpha, _ = scalib.alpha(1.0, 0.0, 1.23)
        assert np.allclose(alpha, 1.23)

    def test_background_calib(self):
        bcalib = BackgroundCalib()
        assert np.allclose(bcalib.rate(0.5), 1.)
        beta, _ = bcalib.beta(1.23)
        assert np.allclose(beta, 1.23)
