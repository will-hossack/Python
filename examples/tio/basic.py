"""
   basic use of the tio input/output module

Author: Will Hossack, The Univesrity of Edinburgh
"""

import tio as t
import math

def main():
    #     Get a float with default and range checking
    f = t.getFloat("Give a float",3.0,0.0,5.0)
    t.tprint("Float is : ", f)

    #     Get a logical within a if statemnet
    if t.getBool("Logical",True):
        print("True")
    else:
        print("False")

    #     Get an integer with no checking
    i = t.getInt("Integer Value",min=5,max=10)
    t.tprint("Integer Value is : ", i)

    c = t.getComplex("Complex Number",3+4j,100.0)
    t.tprint("Complex is : ", c)


main()
