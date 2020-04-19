import optics.jones as j
import matplotlib.pyplot as plt
from optics.wavelength import Red,Green,Blue

def main():

    bbeam = j.LinearPolarisedBeam("h",wavelength = Blue)
    gbeam = j.LinearPolarisedBeam("h",wavelength = Green)
    rbeam = j.LinearPolarisedBeam("h",wavelength = Red)

    beam = [bbeam,gbeam,rbeam]
    
    fp = j.LinearPolariser("h")
    mp = j.QuarterWavePlate(order=0)
    bp = j.LinearPolariser("v")
    sys = j.JonesMatrixSystem(fp,mp,bp)


    sys.polarPlot(beam,1)

    plt.show()

main()
