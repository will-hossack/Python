import matplotlib.pyplot as plt
import optics.wavelength as w
import tio as t

def main():
    peak = t.getFloat("Peak")
    width = t.getFloat("Width")

    sp = w.GaussianSpectrum(peak,width,5.0)
    sp.draw()
    plt.show()
 

main()
