from __future__ import print_function, division

import unittest
import os
import tempfile
import shutil

import desietc.calculator

from ..plots import *

import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!


class TestPlots(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Create a temporary directory.
        cls.tmpdir = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        # Remove the directory after the test.
        shutil.rmtree(cls.tmpdir)

    def test_calculator_plot(self):
        """Test that we can create a calculator and plot its state"""
        calc = desietc.calculator.Calculator(
            1., 0.1, 1., 0.1,
            1.0, 0.5, 2000.,
            1.0, 0.5, 2000.,
            0., 10.)
        calc.update_signal(1000., 1., 0.1)
        calc.update_background(900., 0.5, 0.2)
        # Save plot to temporary directory.
        save = os.path.join(TestPlots.tmpdir, 'test.png')
        plot_calculator(calc, tnow=1500., nsamples=10, save=save)
