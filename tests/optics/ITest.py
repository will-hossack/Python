import optics.wavefront as wf
import matplotlib.pyplot as plt
import math


def main():
    wave = wf.WaveFront().fromFile()
    print(str(wave))
    inter = wf.Interferometer(wave,xtilt=0.0,ytilt=0.0)

    inter.draw()
    plt.show()

main()
