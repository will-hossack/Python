import optics.zernike as zernike
import numpy as np
import matplotlib.pyplot as plt
import math
import cmath
import tio as t


def main():

    
    ze = zernike.ZernikeExpansion(1.0,0.0,0.0,0.0,4.0)


    ze.plotOTF()
    plt.show()

    shiftData = np.linspace(0.0,1.0,otfData.size)
    plt.plot(shiftData,otfData)
    plt.xlim(0.0,1.0)
    plt.grid()
    plt.xlabel("Normalised spatial frequency")
    plt.ylabel("OFT")
    plt.title("Plot of OTF")
    #plt.imshow(psf,cmap=plt.cm.gray,extent=(-1.0,1.0,-1.0,1.0))
    plt.show()

    
    im = ze.getImage()
    
    r = np.cos(im)
    i = np.sin(im)
    z = r + 1j*i
    zc = np.conj(z)

    plt.imshow(r,cmap=plt.cm.gray,extent=(-1.0,1.0,-1.0,1.0))
    plt.show()

    horizontal = False

    xsize,ysize = im.shape

    if horizontal:
        shiftSize = xsize
        fullRange = range(0,ysize)
    else:
        shiftSize = ysize
        fullRange = range(0,xsize)
        
    shiftData = np.linspace(0.0,1.0,shiftSize)
    otfData = np.zeros(shiftSize)

    for shift in range(0,shiftSize):
        otf = 0.0
        shiftRange = range(shift,shiftSize)
        for i in shiftRange:
            for j in fullRange:
                ish = i - shift
                if horizontal:
                    zr = z[i,j]
                    zl = zc[ish,j]
                else:
                    zr = z[j,i]
                    zl = zc[j,ish]
                if not (cmath.isnan(zr) or cmath.isnan(zl)) :
                        otf += (zr * zl).real
        otfData[shift] = otf

    max = otfData[0]
    otfData /= max
    #psf = ze.getPSF(log=False)
    #otf = np.fft.fft2(psf)
    #otf = np.fft.fftshift(otf)
    #otf = abs(otf)
    #line = psf[128]
    
    plt.plot(shiftData,otfData)
    plt.xlim(0.0,1.0)
    plt.grid()
    plt.xlabel("Normalised spatial frequency")
    plt.ylabel("OFT")
    plt.title("Plot of OTF")
    #plt.imshow(psf,cmap=plt.cm.gray,extent=(-1.0,1.0,-1.0,1.0))
    plt.show()

main()
