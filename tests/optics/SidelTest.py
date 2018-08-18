

import tio as t
import optics.wavefront as f
import matplotlib.pyplot as plt

def main():
    sa = [0.0,1.0,-2.0,-1.0,0.5,0.7]
    #sa = [4.0,0.0,0.0,0.0,0.0,0.0]
    zc = [1.0,0.0,0.0,2.0]
    pc = [0.0,0.0,0.0,4.0,0.0,4.0]

    #sid =f.ZernikeWaveFront(zc)
    #sid = f.SeidelWaveFront(sa,theta=0.7,radius=4.0)
    sid = f.PolynomialWaveFront(pc,radius = 3.0)
    t.tprint(repr(sid))

    f.Interferometer().addWaveFront(sid,ytilt=4.0)
    plt.show()

main()
