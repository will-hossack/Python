"""
      Test of spectrums
"""

import optics.wavelength as w
import matplotlib.pyplot as plt

def main():

    s = w.PhotopicSpectrum(10)

    s.draw()
    plt.legend(loc="lower left")
    plt.show()

main()
