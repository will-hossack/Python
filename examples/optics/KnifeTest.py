"""
    Simple appliaction to perform a Focul Knife Edge test of
    a specified lens .

Author: Will Hossack, The University of Edinburgh
"""
import optics.lens as len
import optics.analysis as anal
import matplotlib.pyplot as plt
import math
import tio

def main():

    #         Get the lens from database
    lens = len.DataBaseLens()
    angle = math.radians(tio.getFloat("Angle in degrees"))
    opt = tio.getBool("At optimal focus",True)
    kt = anal.KnifeEdgeTest(lens,angle)
    kt.setKnife(0.0,math.radians(90))
    kt.getImage(opt).draw()
    plt.show()

main()
    
