"""
===
etc
===

The DESI Exposure Time Calculator

Main interface to other parts.

"""

### These three modules are our window to the outside world
## ETC Targest contains information about sources we are observing,
## magnitudes, etc.
import ETCTargets as ETCt
## GuiderPack is a packaing of a single packet of information coming
## from GFA
import ETCGuiderPack as ETCgp
## ETC results contains the results of our processing. This is the
## object that can tell you the completion(t) for the past as well
## as extrapolates for the future. You can also ask it to estimate time to
## completion.
import ETCResults import ETCres

### These following two modules are doing the grunt of analysis
## GuiderPackProcessor implements getting from Guider Images
## to our estimate of disgnal/dt and dnoise/dt
import GuiderPackProcessor as gpp
## Targets Processor impements getting from Targets and a
## tripled of (signal,noise,time) to completion.
import TargetsProcessor as tp

import scipy
from __future__ import division
from __future__ import print_function

class ETC:
    """ This is the main Exposure Time Calculator module.
    An instance of ETC eats GuiderPack objects,
    processes them and calculate relevant quanties, such as
    current level of completeness and time to end.
    """

    verbose=True
    
    def __init__ (self):
        self.integrating=False
        ## These three lists holds the result so far of our
        ## integrated signal and noise
        ## The units of Sig and Noise are ETC unites, i.e. how much a synthethic 
        ## fiber with the same bandpass as guider would accumulate
        self.times=[]
        self.dSigDt=[]
        self.dNoiseDt=[]

        # GuiderPackProcessor is an auxiliary object processing GuiderPacks and returing
        # SNR estimates
        self.gpp = gpp.GuiderPackProcessor()

    def StartIntegration(self, timestamp, targets):
        """ This function starts the integration 
        We need to supply current timestamp and targets.
        targets is of type ETCTargets and contains the list of objects and all relevant
        information to 

        """
        self.integrating=True
        ## We will later replace dSigDt with real values after we have the first estimate
        self.times=[timestamp]
        self.dSigDt=[0.0]
        self.dNoiseDt[0.0]
        self.completion=[0.0]
        ## Targets are used to initialize TargetsProcess
        self.tp=tp.TargetsProcessor(targets)
        
    def StopIntegration(self,timestamp):
        """ This is the opposite of start integration. We call this at the end"""
        self.integrating=False
        ## We use the latest estimate of SNR
        self.times.append(timestamp)
        self.dSigDt.append(dSigDt[-1])
        self.dNoiseDt.append(dNoiseDt[-1])
        ## update our results
        self.UpdateResults(timestamp)

    def ProcessGFAPacket(self, timestamp, GPack):
        """ This processes the GFA object and updates our current state.
            GPack is an object of type GuiderPack
        """
        if not self.integrating:
            return
        self.times.append(timestamp)
        ## Get the current derivatives of noise and signal from GPProcessor
        dSigDt,dNoiseDt =  self.gpp.Estimate(GPack)
        self.times.append(timestamp)
        self.dSigDt.append(dSigDt)
        self.dNoiseDt.append(dNoiseDt)
        if (len(self.times)==2):
            # This is our first estimate. Let's extrapolate our current
            # quantities to time zero
            self.dSigDt[0]=self.dSigDt[1]
            self.dNoiseDt=self.dNoiseDt[1]
        # Update results of where we are at the moment 
        self.UpdateResults(timestampe)
        
        self.PushResults()
            
    def UpdateResults(self,timestamp):
        """ This integrates the signa and noise in our units 
            and produces a results object
        """
        ### whatever the approraite routiner here
        Sig=scipy.integrate(self.times, self.dSigDt)
        Noise=scipy.integrate(self.times, self.dNoiseDt)
        timepassed=self.times[-1]-self.times[0]
        
        ## what is the completion at the moment? TargetProcessor can
        ## convert our sig,noise, timepassed into "completion" number

        completion=self.tp.Completion(Sig,Noise,timepassed)
        self.completion.append(completion)
        ## hopefully all three should now have the same size
        assert(len(self.completion)==len(self.dSigDt)) #etc.

        # if not complete, predict
        if (completion<1.0):
            ### Let's calculate predicted completion
            ### take average timestamp
            dt=time/len(self.times)
            ## we use latest estimate
            dSig=self.dSigDt[-1]*dt
            dNoise=self.dNoiseDt[-1]*dt
            ctime=self.times[-1]

            ## extrapolate quantitties
            etimes=[]
            ecompletion=[]

            while True:
                ctime+=dt
                cSig+=dSig
                cNoise+=dNoise
                cCompletion=self.tp.Completion(cSig,cNoise,ctime-self.times[0])
                etimes.append(ctime)
                ecompletion.append(cCompletion)
                ## stop when either complete or you've been predicting for one hour ahead
                if (cCompletion>1.0) or (len(etimes)>3600):
                    break
        else:
        # we are complete, just have name
            etimes=None
            ecompletion=None

        self.results=ETCres.ETCResults(timestamp, completion, self.times,self.completion, self.etimes, self.ecompletion)

    def PushResults(self):
        if (hasattr(self,"results")):
            ## push self.results
            pass
        
