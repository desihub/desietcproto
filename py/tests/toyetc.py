#!/usr/bin/env python
import sys
sys.path.append("../etc")

import etc
import MockGuiderPackProcessor
import MockTargetsProcessor

## Start our etc
## Our instance
test=etc.ETC(gpp=MockGuiderPackProcessor.GuiderPackProcessor())
## Start integration

time=0
test.StartIntegration(time, None, tp=MockTargetsProcessor.TargetsProcessor())
while True:
    time+=1
    ## send a fake packet
    gfapack=None
    test.ProcessGFAPacket(time,gfapack)
    ## pull results out
    r=test.results
    ## print where we are (including projection)
    r.printStatus()
    ## if we are done, break
    if (r.complete()):
        break
test.StopIntegration()
print "Done"
