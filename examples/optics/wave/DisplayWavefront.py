"""
Exaample program to read in a wavefront from a file, display the contents 
and plot display the wavefront an an interferogram with optional x/y tilts
"""

import optics.wavefront as wf
import matplotlib.pyplot as plt
import tio as t


def main():
    #       Read the wavefront in from a .wf file and display it contents
    wave = wf.WaveFront().fromFile()
    t.tprint(repr(wave))
    
    #       Get the tilts
    xt = t.getFloat("Xtilt",0.0)
    yt = t.getFloat("Ytilt",0.0)
    
    #       Make the interferometed
    inter = wf.Interferometer(wave)
    inter.setTilt(xt,yt)             # Set the tilts
    inter.draw()                     # Display
    plt.show()

main()
