""""
=========
ETCResuts
=========

This is the object containing the results of current analysis to be pushed
to OCS

"""

import scipy
from scipy.interpolate import interp1d
from __future__ import division
from __future__ import print_function

class ETCResults:
    def __init__ (self, timestamp, completion, recordedTimes,
                  recordedCompletion, extrapolatedTimes, extrapolatedCompletion):
        self.timestamp=timestamp
        self.completion=completion
        self.times=recodedTimes
        self.completion=recordedCompletion
        self.etimes=extrapolatedTimes
        self.ecompletion=extrapolatedCompletion


    def timeToCompletion(self, target=1.0):
        """ 
        Tries to estimate when are we going to be complete, returns None, if it can't do that
        """
        ### are we already done by any chance?
        if (self.completion>=target):
            return self.timestamp
            
        ## Now, let's hope we can actually do it
        if (target>self.self.ecompletion[-1]):
            return None
            
        if (target<self.self.ecompletion[0]):
            # This doesn't make any sense, either we are complete,
            # we can predictit, or it is too far in the future to predict
            stop()

        ## we interpolated etimes as a function of completion and estimates time
        ## when we finish
        tend=intedp1d(self.ecompletion,self.etimes)(target)

        ## in principle, we could estimate the error on this quantity (and also
        ## the etc module could estimate the errors, but we know nobody is going to bother
        ## taking that into account

        return tend
