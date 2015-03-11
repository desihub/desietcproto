"""
=======================
MockGuiderPackProcessor
=======================

This is like GuiderPackProcessor, only that it returns fake values

"""

from __future__ import division
from __future__ import print_function

class GuiderPackProcessor:
    def __init__(self, sig0=10, dsigdt=0.1, noise0=5; dnoisedt=-0.1):
        """ Mock data. We have a very good observation starting with sig0 and noise0
        and changing by dsigdt and dnoisedt in each step """
        self.sig=sig0
        self.noise=noise0
        self.dsigdt=dsigdt
        self.dnoisedt=dnoisedt

    def Estimate (gpack):
        """ Takes a GuiderPack and returns an 
        estimate of curent dSigDt and DNoiseDt  in a tuple"""
        toret=(self.sig, self.noise)
        self.sig+=self.dsigdt
        self.noise+=self.dnoisedt
        return toret
