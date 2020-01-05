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
    lens = len.DataBaseLens("$LENS/Linos140Doublet")
    angle = math.radians(tio.getFloat("Ray Angle in degrees"))
    opt = tio.getInt("Focal Option",1)
    kangle = math.radians(tio.getFloat("Knife angle in degrees"))
    kt = anal.KnifeTest(lens,angle)
    kt.setWire(True)
    knife = tio.getFloat("Knife")
    kt.setKnife(knife,kangle,shift=-.02)
    kt.getImage(opt).draw()
    plt.show()

main()
    
