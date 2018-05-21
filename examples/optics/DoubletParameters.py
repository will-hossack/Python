"""
View an investigate double parameters.

Author: Will Hossack, The University of Edinburgh
"""
import optics.lens as len
import tio as t
import matplotlib.pyplot as plt

def main():
    #
    #      Create default doublet
    #
    lens = len.Doublet()

    t.tprint("Power is : ",1.0/lens.focalLength())
    t.tprint("Petzal sum is : ",lens.petzvalSum())

    first = lens.paraxialMatrix(first=0,last=1)
    second = lens.paraxialMatrix(first=1,last=2)

    t.tprint("First matrix ", first)
    t.tprint("Second matrix", second)

    fp = first.backPower()
    sp = second.backPower()
    t.tprint("First power is : ",fp)
    t.tprint("Second power is : ",sp)
    t.tprint("Total power is ", fp - sp)

main()
