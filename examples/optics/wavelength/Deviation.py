"""
Programme to calcualte and plot the angle of minumum deviations 
for a prism of specified material against wavelength
"""

import optics.wavelength as w
import matplotlib.pyplot as plt
import math
import tio
import numpy as np


def deviation(pa,n):
    """
    Deviation angle function
    
    :param pa: prism angle in degrees
    :param n: refratcive index
    """
    a = math.radians(pa)     # Pa in radians
    sa = n*math.sin(a/2)
    s = math.asin(sa)
    return math.degrees(2*s - a)



def main():

    #        Get prism angle an material
    prismAngle = tio.getFloat("Prism angle in degrees",60,10,90)
    index = w.MaterialIndex()
    
    wData = np.linspace(w.BlueLimit,w.RedLimit,100)    # wavelength data
    devData = np.empty(len(wData))                     # deviation angle

    #          Fill up the deviation array
    for i,wave in enumerate(wData):
        devData[i] = deviation(prismAngle,index.getValue(wave))

    #          Do the plotting
    plt.plot(wData,devData,"g")
    plt.title("Angle of Minium Deviation for {0:s}".format(index.title))
    plt.xlabel("Wavelength")
    plt.ylabel("Angle in degrees")
               
    plt.show()
    
        
main()
