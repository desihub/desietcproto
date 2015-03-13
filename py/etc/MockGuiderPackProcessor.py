"""
=======================
MockGuiderPackProcessor
=======================

This is like GuiderPackProcessor, only that it returns fake values

"""

from __future__ import division
from __future__ import print_function

class GuiderPackProcessor:
    def __init__(self, sig0=5, dsigdt=0.1, noiseVar0=1, dnoiseVardt=+0.05):
        """ Mock data. We have a very good observation starting with sig0 and noiseVar0
        and changing by dsigdt and dnoiseVardt in each step """
        self.sig=sig0
        self.noiseVar=noiseVar0
        self.dsigdt=dsigdt
        self.dnoiseVardt=dnoiseVardt

    def Estimate (self,gpack):
        """ Takes a GuiderPack and returns an 
        estimate of curent dSigDt and DNoiseVarDt  in a tuple"""
        toret=(self.sig, self.noiseVar)
        self.sig+=self.dsigdt
        self.noiseVar+=self.dnoiseVardt
        return toret
