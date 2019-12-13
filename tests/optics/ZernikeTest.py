"""      Test for Seidel wavefront
"""
from optics.wavefront import ZernikeWaveFront
import tio as t
import matplotlib.pyplot as plt

def main():

    s = ZernikeWaveFront(1.0,1.0,12.0,3.0,4.0,-6.0,4.0)
    t.tprint(s)

    s.draw(xtilt=3.0)
    plt.show()

    im = s.plotPSF(log=False)
    plt.show()

main()
