from __future__ import print_function, division

import unittest

import desietc.calculator

from ..plots import *

import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!


class TestPlots(unittest.TestCase):

    def test_calculator_plot(self):
        """Test that we can create a calculator and plot its state"""
        calc = desietc.calculator.Calculator(
            1., 0.1, 1., 0.1,
            1.0, 0.5, 2000.,
            1.0, 0.5, 2000.,
            0., 10.)
        calc.update_signal(1000., 1., 0.1)
        calc.update_background(900., 0.5, 0.2)
        plot_calculator(calc, tnow=1500.)
