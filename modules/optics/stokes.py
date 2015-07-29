import math
import cmath
import wavelength as wl
from opticalray import Ray
from matplotlib.pyplot import polar,show
from jones import *
#
#                 Python classes to implement Stokes and Muller
#                 polarsiation classes.



class StokesVector(list,Ray):
    #
    #            Constuctor for a Stokes Vector
    #
    def __init__(self,s0_or_sv,s1=None,s2=None,s3=None,wavelength = wl.Green):
        Ray.__init__(self,wavelength) # Set inderlying class
        self = [0.0]*4                # The 4 paraeters
        if isinstance(s0_or_sv,list) :
            for i in s0_or_sv:
                self[i] = s0_or_sv[i]
        else:
            self[0] = float(s0_or_sv)
            self[1] = float(s1)
            self[2] = float(s2)
            self[3] = float(s3)
            

    #        String to return Stokes vector
    def __str__(self):
        return "StokesVector : [" + str(self[0]) + " , " + str(self[1]) +\
            " , " + str(self[2]) + " , " + str(self[3]) + "] w: " +\
        str(self.wavelength)

