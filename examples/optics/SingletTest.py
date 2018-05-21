"""
   Test the methids to read in a Singlet
"""
import optics.lens as l
import tio as t

def main():

    lens = l.Singlet()
    t.tprint("Lens is : ",lens)

    nl = l.Singlet().setFromString("point: (5.0,-5.0,20) index: F4 focal: 50 fno: 4 bend planoconvex")
    t.tprint("Updated Lens is : ",nl)

main()
