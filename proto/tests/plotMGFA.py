#!/usr/bin/env python
import sys
import numpy as np

sys.path.append("../etc")
sys.path.append("etc")

from MockGFA import MockGFA
import pylab

rmags=np.linspace(18,15,8)
gfagp=MockGFA(Nstars=8, Nskies=8,star_rmags=rmags).Yield(0.0)

mn=min([x.array.min() for x in gfagp.stars+gfagp.skies])
mx=max([x.array.max() for x in gfagp.stars+gfagp.skies])

pylab.figure(figsize=(12,12))
for i in range(8):
    pylab.subplot(4,4,1+i)
    pylab.imshow(gfagp.stars[i].array,vmin=mn,vmax=mx,interpolation='nearest', cmap=pylab.cm.binary)
    pylab.xlabel('r=%g'%rmags[i])
    pylab.subplot(4,4,9+i)
    pylab.imshow(gfagp.skies[i].array,vmin=mn,vmax=mx,interpolation='nearest',cmap=pylab.cm.binary)
    pylab.xlabel('sky')
pylab.show()

