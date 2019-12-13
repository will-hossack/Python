"""      Test for Seidel wavefront
"""
from optics.wavefront import SeidelWaveFront
import tio as t
import matplotlib.pyplot as plt

def main():
    coef =[12.0,0.0,0.0,8.0,6.0,5.0]

    s = SeidelWaveFront(1.0,1.0,coef)
    t.tprint(s)

    s.draw(xtilt=3.0)
    plt.show()

    im = s.plotPSF(log=False)
    plt.show()

main()
