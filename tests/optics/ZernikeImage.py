from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import matplotlib.pyplot as plt
#from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import optics.zernike as z
import tio as t

def main():



    exp = z.ZernikeExpansion(1.0,0.55,[1.0,0.0,-0.5,-1.0])
    im = exp.getImage(32)

    x = np.linspace(-1,1,32)
    y = np.linspace(-1,1,32)
    X,Y = np.meshgrid(x,y)

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    surf = ax.plot_surface(X, Y, im)   # cmap=cm.coolwarm,linewidth=0, antialiased=False)
    
    plt.show()

main()
    
