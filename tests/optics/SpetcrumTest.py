"""
   Test for Spectrum
"""

import optics.wavelength as w
import sys
import matplotlib.pyplot as plt

def main() :

    s = w.TriColourSpectrum(0.5,0.8,0.2)
    
    print(repr(s))
    s.draw()
    plt.legend()
    plt.title(repr(s))
    plt.show()
    




main()
