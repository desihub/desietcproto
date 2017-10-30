#!/usr/bin/env python
import sys
sys.path.append("../etc")
sys.path.append("etc")

import etc
import MockGuiderPackProcessor
import MockTargetsProcessor

## Start our etc
## Our instance
test=etc.ETC(gpp=MockGuiderPackProcessor.GuiderPackProcessor())
## Start integration

time=0.0
test.StartIntegration(time, None, tp=MockTargetsProcessor.TargetsProcessor())
while True:
    time+=1.0
    ## send a fake packet
    gfapack=None
    test.ProcessGFAPacket(time,gfapack)
    ## pull results out and print status (including proj)
    test.results.printStatus()
    ## if we are done, break
    if (test.results.complete()):
        break
time+=0.1
test.StopIntegration(time)
test.results.printStatus()
print "Done"
