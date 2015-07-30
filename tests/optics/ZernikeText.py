import numpy as np
import matplotlib.pyplot as plt
import optics.zernike as zern
import tio

def main():
    
    z = zern.ZernikeExpansion(6.0,0.566,[1.0,0.0,0.0,1.0,-2.0,3.0])
    z.draw()
    plt.show()

main()
    
