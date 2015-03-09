"""
================
TargetsProcessor
================

Given targets, the main job of this piece of shit is to 
convert the magnitudes and object types in the ETCTargets
into a single number called "Completion", between 0 and 1
"""

from __future__ import division
from __future__ import print_function


class TargetsProcessor:
    def __init__ (self, targets):
        self.targets=targets

    def Completion(self,Sig,Noise,Time):
        """ This is where the devil resides. Given our self.targets
        and ETC normalized values of Sig, Noise and Time, return 
        completion, which is a number between 0.0 and 1.0.
        """
        pass
        return 0.0

    
