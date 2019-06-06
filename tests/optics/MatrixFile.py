"""
   Test for materials 
"""

import optics.matrix as m
import matplotlib.pyplot as plt
import tio as t


def main() :

    file = t.openFile("Matrix",defaulttype="matrix")
    pg = m.DataBaseMatrix(file)

    t.tprint(pg.getInfo())

    pg.draw()
    plt.show()
    




main()
