"""      Test for Seidel wavefront
"""
from optics.wavefront import SeidelWaveFront
import tio as t
import matplotlib.pyplot as plt

def main():
    coef =[0.0,-6.0,0.0,2.0,0.0,0.0]

    s = SeidelWaveFront(1.0,0.5,1.0,coef)
    t.tprint(repr(s))

    s.plotImage(xtilt=1.0)
    plt.show()

    im = s.plotPSF(log=False)
    plt.show()

main()
