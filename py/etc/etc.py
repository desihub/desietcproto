"""
===
etc
===

The DESI Exposure Time Calculator

Main interface to other parts.

"""

from __future__ import division
from __future__ import print_function
import scipy,math
import scipy.integrate

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
import ETCResults as ETCres

### These following two modules are doing the grunt of analysis
## GuiderPackProcessor implements getting from Guider Images
## to our estimate of disgnal/dt and dnoisevar/dt
import GuiderPackProcessor as gpp
## Targets Processor impements getting from Targets and a
## tripled of (signal,noiseVar,time) to completion.
import TargetsProcessor as tp


class ETC:
    """ This is the main Exposure Time Calculator module.
    An instance of ETC eats GuiderPack objects,
    processes them and calculate relevant quanties, such as
    current level of completeness and time to end.
    """

    verbose=True
    
    def __init__ (self, gpp=None):
        self.integrating=False
        ## These three lists holds the result so far of our
        ## integrated signal and noiseVar
        ## The units of Sig and NoiseVar are ETC unites, i.e. how much a synthethic 
        ## fiber with the same bandpass as guider would accumulate
        self.times=[]
        self.dSigDt=[]
        self.dNoiseVarDt=[]

        # GuiderPackProcessor is an auxiliary object processing GuiderPacks and returing
        # SNR estimates
        ## we allow the caller to supply its wn GPP (for testing purposes, for example):
        if (gpp):
            self.gpp=gpp
        else:
            self.gpp = gpp.GuiderPackProcessor()


    def StartIntegration(self, timestamp, targets, tp=None):
        """ This function starts the integration 
        We need to supply current timestamp and targets.
        targets is of type ETCTargets and contains the list of objects and all relevant
        information to 

        """
        self.integrating=True
        ## We will later replace dSigDt with real values after we have the first estimate
        self.times=[timestamp]
        self.dSigDt=[0.0]
        self.dNoiseVarDt=[0.0]
        self.completion=[0.0]
        ## Targets are used to initialize TargetsProcess
        if tp:
            self.tp=tp
        else:
            self.tp=tp.TargetsProcessor(targets)
        
    def StopIntegration(self,timestamp):
        """ This is the opposite of start integration. We call this at the end"""
        self.integrating=False
        ## We use the latest estimate of SNR
        self.times.append(timestamp)
        self.dSigDt.append(self.dSigDt[-1])
        self.dNoiseVarDt.append(self.dNoiseVarDt[-1])
        ## update our results
        self.UpdateResults(timestamp)

    def ProcessGFAPacket(self, timestamp, GPack):
        """ This processes the GFA object and updates our current state.
            GPack is an object of type GuiderPack
        """
        if not self.integrating:
            return
        ## Get the current derivatives of noiseVar and signal from GPProcessor
        dSigDt,dNoiseVarDt =  self.gpp.Estimate(GPack)
        self.times.append(timestamp)
        self.dSigDt.append(dSigDt)
        self.dNoiseVarDt.append(dNoiseVarDt)
        if (len(self.times)==2):
            # This is our first estimate. Let's extrapolate our current
            # quantities to time zero
            self.dSigDt[0]=self.dSigDt[1]
            self.dNoiseVarDt[0]=self.dNoiseVarDt[1]
        # Update results of where we are at the moment 
        self.UpdateResults(timestamp)
        
        self.PushResults()
            
    def UpdateResults(self,timestamp):
        """ This integrates the signa and noiseVar in our units 
            and produces a results object
        """
        ### whatever the approraite routiner here
        Sig=scipy.integrate.simps(self.dSigDt,self.times)
        NoiseVar=scipy.integrate.simps(self.dNoiseVarDt,self.times)
        timepassed=self.times[-1]-self.times[0]
        
        ## what is the completion at the moment? TargetProcessor can
        ## convert our sig,noiseVar, timepassed into "completion" number

        completion=self.tp.Completion(Sig,NoiseVar,timepassed)
        self.completion.append(completion)
        ## hopefully all three should now have the same size
        assert(len(self.completion)==len(self.dSigDt)) #etc.

        # if not complete, predict
        if (completion<1.0):
            ### Let's calculate predicted completion
            ### take average timestamp
            dt=timepassed/(len(self.times)-1)/10.0
            ## we use latest estimate for actual values
            dSig=self.dSigDt[-1]*dt
            dNoiseVar=self.dNoiseVarDt[-1]*dt
            ctime=self.times[-1]
            ## but variance for uncertainty
            dSigErr=math.sqrt(scipy.array(self.dSigDt).var())
            dNoiseVarErr=math.sqrt(scipy.array(self.dNoiseVarDt).var())

            ## extrapolate quantitties
            etimes=[]
            ecompletion=[]
            ecompletionP=[]
            ecompletionM=[]

            #forecasted signal
            cSig=Sig
            cNoiseVar=NoiseVar

            #forecasted signal+uncert. Let's do 1/sqrt(2) +sig -noiseVar 
            #i.e. 1 sigma fast
            dSigP = dSig + (dSigErr)/math.sqrt(2.0)*dt
            dNoiseVarP = dNoiseVar - (dNoiseVarErr)/math.sqrt(2.0)*dt
            cSigP=Sig
            cNoiseVarP=NoiseVar
            #forecasted signal-uncert
            dSigM = dSig - (dSigErr)/math.sqrt(2.0)*dt
            dNoiseVarM = dNoiseVar + (dNoiseVarErr)/math.sqrt(2.0)*dt
            cSigM=Sig
            cNoiseVarM=NoiseVar


            while True:
                ctime+=dt
                cSig+=dSig
                cNoiseVar+=dNoiseVar
                cCompletion=self.tp.Completion(cSig,cNoiseVar,ctime-self.times[0])
                cSigP+=dSigP
                cNoiseVarP+=dNoiseVarP

                cCompletionP=self.tp.Completion(cSigP,cNoiseVarP,ctime-self.times[0])
                cSigM+=dSigM
                cNoiseVarM+=dNoiseVarM
                cCompletionM=self.tp.Completion(cSigM,cNoiseVarM,ctime-self.times[0])

                #print (ctime, cSig,cNoiseVar, cCompletion, cCompletionP, cCompletionM)

                etimes.append(ctime)
                ecompletion.append(cCompletion)
                ecompletionP.append(cCompletionP)
                ecompletionM.append(cCompletionM)

                ## stop when either complete or you've been predicting for one hour ahead
                if (cCompletionM>=1.0) or (len(etimes)>3600):
                    break
        else:
        # we are complete, just have name
            etimes=None
            ecompletion=None
            ecompletionP=None
            ecompletionM=None

        self.results=ETCres.ETCResults(timestamp, completion, Sig,
                                       NoiseVar, self.times,
                                       self.completion, etimes,
                                       ecompletion, ecompletionP,
                                       ecompletionM )

    def PushResults(self):
        if (hasattr(self,"results")):
            ## push self.results
            pass
        
