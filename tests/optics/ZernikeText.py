import numpy as np
import matplotlib.pyplot as plt
import optics.zernike as z
import tio as t

def main():


    n = t.getInt("n value")
    r = np.linspace(-1.0,1.0,100)

    for m in range(n%2,n + 1,2):
        plt.plot(r,z.radial(n,m,r))
    plt.show()

main()
    
