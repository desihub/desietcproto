""""
=========
ETCResuts
=========

This is the object containing the results of current analysis to be pushed
to OCS

"""

from __future__ import division
from __future__ import print_function
import scipy
import scipy.interpolate  as intp

class ETCResults:
    def __init__ (self, timestamp, completion, ETCSignal, ETCNoise,
                  recordedTimes, recordedCompletion,
                  extrapolatedTimes, extrapolatedCompletion,
                  extrapolatedCompletionP, extrapolatedCompletionM):
        self.timestamp=timestamp
        self.ETCSignal=ETCSignal
        self.ETCNoise=ETCNoise
        self.completion=completion
        self.times=recordedTimes
        self.recordedCompletion=recordedCompletion
        ## extrapolated times
        self.etimes=extrapolatedTimes
        ## extrapolated completions and P/M errors 
        self.ecompletion=extrapolatedCompletion
        self.ecompletionP=extrapolatedCompletionP
        self.ecompletionM=extrapolatedCompletionM



    def timeToCompletion(self, target=1.0):
        """ 
        Tries to estimate when are we going to be complete, returns None, if it can't do that
        The return value is a tripled (tm, t, tp), where we have around 68% confidence, we will
        finish between tm and tp with best guess at t
        """
        ### are we already done by any chance?
        if (self.completion>=target):
            return None
            
        ## Now, let's hope we can actually do it
        if (target>self.ecompletion[-1]):
            return float("inf"), float("inf"), float("inf")
            
        if (target<self.ecompletion[0]):
            # This doesn't make any sense, either we are complete,
            # we can predictit, or it is too far in the future to predict
            stop()



        ## we interpolated etimes as a function of completion and estimates time
        ## This is pointless, since we know one has completion = 1
        #tend=intp.interp1d(self.ecompletion,self.etimes)(target)
        #tendP=intp.interp1d(self.ecompletionP,self.etimes)(target)
        #tendM=intp.interp1d(self.ecompletionM,self.etimes)(target)
        tendP,tend, tendM=None,None,None
        for t, cm,c,cp in zip(self.etimes, self.ecompletionM, 
                              self.ecompletion, self.ecompletionP):
            #print (tendM,cm,t)
            if not tendM and cm>=1.0:
                tendM=t
            if not tend and c>=1.0:
                tend=t
            if not tendP and cp>=1.0:
                tendP=t

        return tendP, tend, tendM

    def printStatus(self):
        print ("Time=",self.timestamp, "ETC SIG/NOISE=",
               self.ETCSignal, self.ETCNoise, " Completion=",
               self.completion,end="") 
        try:
            t1,t2,t3=self.timeToCompletion()
            print (" Expected completion:",
               t1,t2,t3)
        except:
            print ("")

    def complete(self):
        return (self.completion>=1.0)
