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
    angle = math.radians(tio.getFloat("Ray Angle in degrees"))
    opt = tio.getInt("Focal Option",1)
    kangle = math.radians(tio.getFloat("Knife angle in degrees",0.0))
    kt = anal.KnifeTest(lens,angle,opt)
    knife = tio.getFloat("Knife",0.0)
    kt.setKnife(knife,kangle)
    kt.getImage().draw()
    plt.show()

main()
    
