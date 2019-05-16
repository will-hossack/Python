"""
     Test use of invverse square calcualtions
"""

import tio as t
import vector as v
import matplotlib.pyplot as pt

def main():


    xdata = []
    ydata = []
    
    a = t.getVector3d("Central",[5,1,0])
    t.tprint(a)

    for i in range(0,100):
        x = i/10.0
        xdata.append(x)
        b = v.Vector3d(x,0,0)
        f = a.inverseSquare(b,-1.0)
        t.tprint(f),
        ydata.append(abs(f))

    pt.plot(xdata,ydata)
    pt.show()
main()
