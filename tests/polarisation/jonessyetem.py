import optics.jones as j
import matplotlib.pyplot as plt
import tio
import math

def main():

    beam = j.LinearPolarisedBeam("h")
    print(repr(beam))

    
    fp = j.LinearPolariser("h")
    mp = j.QuarterWavePlate()
    bp = j.LinearPolariser("v")
    sys = j.JonesMatrixSystem(fp,mp,bp)


    sys.polarPlot(beam,1)


    plt.show()

main()
