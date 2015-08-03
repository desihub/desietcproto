#!/usr/bin/env python
import sys
sys.path.append("../etc")
sys.path.append("etc")

from MockGFA import MockGFA
import pylab

gfagp=MockGFA(Nstars=4, Nskies=4,star_rmags=17).Yield(0.0)
for i in range(4):
    pylab.subplot(2,4,1+i)
    pylab.imshow(gfagp.stars[i].array)
    pylab.subplot(2,4,5+i)
    pylab.imshow(gfagp.skies[i].array)

pylab.show()

