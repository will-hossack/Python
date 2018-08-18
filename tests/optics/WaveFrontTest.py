

import tio as t
import optics.wavefront as f
import matplotlib.pyplot as plt

def main():


    wf = f.WaveFront().readFromFile()
    
    t.tprint(repr(wf))

    f.Interferometer().addWaveFront(wf,ytilt=4.0)
    plt.show()

main()
