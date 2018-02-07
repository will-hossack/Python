"""
    Simple appliaction to perform a Focul Knife Edge test of
    a specified lens .

Author: Will Hossack, The University of Edinburgh
"""
import optics.lens as len
import optics.analysis as anal
import matplotlib.pyplot as plt
import vector
import math
import tio

def main():

    lens = len.DataBaseLens()
    angle = tio.getFloat("Angle in degrees")
    opt = tio.getBool("At optimal focus",True)
    output = anal.knifeEdgeTest(lens,math.radians(angle),0.0,optimal=opt)
    
    output.draw()
    plt.show()

main()
    
