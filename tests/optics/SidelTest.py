

import tio as t
import optics.wavefront as f
import matplotlib.pyplot as plt

def main():
    #sa = [0.0,1.0,-2.0,-1.0,0.5,0.7]
    #sa = [4.0,0.0,0.0,0.0,0.0,0.0]
    zc = [1.0,0.0,0.0,2.0]

    sid =f.Zernike(zc)
    #sid = f.Seidel(sa,theta=0.7,radius=4.0)
    t.tprint(repr(sid))

    f.Interferometer().addWaveFront(sid,xtilt=4.0)
    plt.show()

main()
