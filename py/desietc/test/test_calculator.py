from __future__ import print_function, division

import unittest

from ..calculator import *


class TestConfig(unittest.TestCase):

    def test_ctor(self):
        calc = Calculator(
            1., 0.1, 1., 0.1,
            1.0, 0.5, 2000.,
            1.0, 0.5, 2000.,
            0., 10.)
