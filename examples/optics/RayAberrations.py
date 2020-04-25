"""
    Example program to plot the Ray Aberrations for a lens
Author: Will Hossack, The Unievsrity of Edinburgh
"""
import optics.lens as len
import optics.wavelength as wl
import matplotlib.pyplot as plt
import optics.wavefront as a
import math
import tio

def main():

    # Get lens and other info

    lens = len.DataBaseLens()
    angle = math.radians(tio.getFloat("Angle in degrees"))
    wave = tio.getFloat("Wavelength of plot",wl.Default)

    wa = a.WaveFrontAnalysis(lens)
    wa.drawAberrationPlot(angle,wave=wave,legend = "lower left")

    
    plt.show()

main()
            
    

