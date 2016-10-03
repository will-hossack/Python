import matplotlib.pyplot as plt
import optics.wavelength as w
import tio as t

def main():
    
    sp = w.PlanckSpectrum()
    sp.minWavelength = 0.15
    tr = [10000.0,8000.0,6000.0,5000.0,4000.0,3000.0]
    for t in tr:
        sp.setTemperature(t)
        sp.draw()
    plt.show()
 

main()
