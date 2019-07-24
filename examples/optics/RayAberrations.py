"""
    Example program to plot the Ray Aberrations for a lens
Author: Will Hossack, The Unievsrity of Edinburgh
"""
import optics.lens as len
import optics.wavelength as wl
import matplotlib.pyplot as plt
import optics.analysis as an
import math
import tio

def main():

    # Get lens and other info

    lens = len.DataBaseLens()
    angle = math.radians(tio.getFloat("Angle in degrees"))
    wave = tio.getFloat("Wavelength of plot",wl.Default)

    # plot = an.aberrationPlot(lens,angle,wave)
    abb = an.AberrationPlot(lens)
    abb.draw(angle)

    plt.xlim(-1.0,1.0)
    plt.legend(loc="lower right",fontsize="small")
    plt.grid()
    plt.show()

main()
            
    

