import matplotlib.pyplot as plt
import optics.analysis as a
import numpy as np

def main():

    image = a.ColourImage().addTestGrid(12,6,intensity=[1.0,1.0,0.0])
    image.draw()
    plt.show()

main()
