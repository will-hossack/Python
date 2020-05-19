"""
       Set of utils used in the optics package
"""
import numpy as np
from scipy.special import j1
import matplotlib.pyplot as plt



def sinc(x):
    """
    A simple sinc (witout the stray pi sin Scipy)
    """
    if x == 0.0 :
        return 1.0
    else:
        return np.sin(x)/x

def jinc(x):
    """
    The jinc() normalsied to 1.0
    """
    if x.any() == 0.0:
        return 1.0
    else:
        return 2*j1(x)/x


def main():
    #s = sinc(0.0)
    #print(str(s))
    #j = jinc(0.0)
    #print(str(j))
    xData = np.linspace(0.0,10,300)
    
    plt.plot(xData,sinc(xData))
    #plt.plot(xData,jinc(xData))
    plt.show()
    
    
main()