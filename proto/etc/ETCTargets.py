"""
==============
ETCTargets
=============

A class containing relevant information that allows us to estimate time to completion

"""
from __future__ import division
from __future__ import print_function
import astropy.io.fits as pyfits

class ETCTargets(object):

    def __init__(self, fibermapfile):
        """ Initializes ETCTargets from a fibermap file. 
        If fibermapfile is a string -> open fits file
        otherwise assume it is already a fits file handle"""
        if type(fibermapfile)==type("string"):
            fibermapfile=pyfits.open(fibermapfile)
        self.objects=[]
        for line in fibermapfile[1].data:
            objtype=line['OBJTYPE']
            if objtype in ['LRG','ELG','QSO']:
                self.objects.append((objtype,line['MAG']))
            else:
                if not (objtype in ['SKY','STD']):
                    print("Warning: ETCTargets does not know about objtype=",objtype)
        print("ETCTargets initialized with ",len(self.objects), " targets")
