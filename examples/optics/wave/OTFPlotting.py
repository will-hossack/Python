"""
      Explore the OFT plotting

"""
import scipy.integrate as integrate
import optics.wavefront as wf
import math
import tio as t
import matplotlib.pyplot as plt
import numpy as np

wave = None


def fntwo(x,y,a,b):
    
    z = wave.getComplexValue(x,y)
    zs = wave.getComplexValue(x + a,y + b)
    val = (z*zs.conjugate()).real
    if math.isnan(val):
        return 0.0
    else:
        return 1.0
    
def fnone(x,s):
    z = wave.getComplexValue(x,0.0)
    zs = wave.getComplexValue(x + s,0.0)
    al = (z*zs.conjugate()).real
    if math.isnan(val):
        val = 0.0
    return float(val)

def shear(x,y,r,a,b):
    rsqr = r*r
    if x*x + y*y > rsqr:
        return 0.0
    
    xp = x + a
    yp = y + b
    if xp*xp + yp*yp > rsqr:
        return 0.0
    
    #z = wave.getComplexValue(x,y)
    #zs = wave.getComplexValue(xp, yp)
    #    val = (z*zs.conjugate()).real
    #val = z.real*zs.real + z.imag*zs.imag
    
    p = wave.getValue(x,y)
    ps = wave.getValue(xp,yp)
    val = math.cos(p - ps)
    
    if math.isnan(val):
        return 0.0
        print("nan Trap neeed")
    else:
        return val

def main():
    global wave
    wave = wf.WaveFront().fromFile()
    print(repr(wave))
    
    wave.plotOTF(100,"b")
    
    """sData = np.linspace(0,2.0*wave.radius,128)
    oData = np.zeros(sData.size)
    
    xyRange = np.linspace(-wave.radius,wave.radius,50)
    
    for i,s in enumerate(sData):
        otf = 0.0
        for y in xyRange:
            for x in xyRange:
                otf += shear(x,y,wave.radius,s,0.0) 


        #           area = integrate.dblquad(shear,-radius,radius,-radius,radius,args = (radius,x,y))
        oData[i] = otf
    
    oData /= oData[0]
        
        
        
        
        
        #otf = integrate.quad(fnone,-wave.radius,wave.radius,args= (a))
    
        #otf = integrate.dblquad(fntwo,-wave.radius,wave.radius,lambda x: -wave.radius,lambda x: wave.radius,args = (x,y))
    
        #print("Val : " + str(otf))
    
        #oData[i] = otf[0]
        
    
    plt.plot(sData,oData)
    """
    plt.show()
    
        
    
main()
