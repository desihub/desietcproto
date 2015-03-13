"""
=====================
MockTargetsProcessor
=====================

Like TargetsProcessor, but mock
"""

from __future__ import division
from __future__ import print_function
import math

class TargetsProcessor(object):
    def __init__ (self):
        pass

    def Completion(self,Sig,NoiseVar,Time):
        ## lets have something trivial
        TNoiseVar=NoiseVar+0.1*Time
        SNR=Sig/math.sqrt(TNoiseVar)
        ## Lets our completion be 0 at SNR=20 and 1 at SNR=30
        if (SNR<20):
            return 0.0
        elif (SNR>30):
            return 1.0
        else:
            return (SNR-20.)/10.
