

import tio as t
import optics.wavefront as f
import matplotlib.pyplot as plt

def main():
    #a = [0.0,1.0,-2.0,-1.0,0.5,0.7]
    #a = [4.0,0.0,0.0,0.0,0.0,0.0]
    a = [4.0,0.0,0.0,0.0]

    wf = f.SeidelWaveFront(a)
    t.tprint(repr(wf))

    psf = f.ScalarPSF()
    psf.addWaveFront(wf)
    psf.draw()
    plt.show()
main()
