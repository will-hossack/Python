"""
   Read in a matrix from a from a DataBase matrix lens, scale focal
   Length and, and display planes in a diagram with equal axis scale.

"""

import optics.matrix as m
import matplotlib.pyplot as plt
import tio as t


def main() :

    file = t.openFile("Matrix",defaulttype="matrix")
    pg = m.DataBaseMatrix(file)    # Read in the file
    t.tprint(pg.getInfo())

    fl = t.getFloat("Set focal Length to", pg.backFocalLength())
    pg.setFocalLength(fl)
    t.tprint(pg.getInfo())


    fig,ax = plt.subplots()        # Get the axis
    pg.draw()                      # Draw
    ax.set_aspect(1.0)             # set unit aspect
    plt.show()
    




main()
