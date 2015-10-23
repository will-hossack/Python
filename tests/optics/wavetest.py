"""
      Test for the wavelength.py module
"""

import optics.wavelength as w
import matplotlib.pyplot as plt

def main():

    print("Default wavelength is : " + str(w.Default))

    print("New call to default is : " + str(w.Default))
    w.setDefaultWavelength(0.66)
    print("New Default wavelength is : " + str(w.Default))

    n = w.SimpleCauchyIndex(562568)

    print(str(n))

    print("nd is : " + str(n.getNd()) + " and Vd is : " + str(n.getVd()))
    n.draw("g")
    plt.show()

main()
