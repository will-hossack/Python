""""
Program to calulate the resolution limit of a prism and compare it to the
Na Double
"""


import optics.wavelength as w
import tio
import matplotlib.pyplot as plt

def main():

    res = w.Sodium_D/(w.Sodium_D2 - w.Sodium_D1)
    tio.tprint("Resolution target for Na Doublet is : ",res)

    index = w.MaterialIndex()
    tio.tprint(repr(index))
    tio.tprint("Nd index : ",index.getNd()," Abbe No: ",index.getVd())

    
    dn = index.getDerivative(w.Sodium_D)
    tio.tprint("dn / d :",dn)
    d = tio.getFloat("d in mm")
    d *= 1000
    pr = d*dn

    tio.tprint("Prism resolution  : ",pr)

    if abs(pr) > abs(res):
        tio.tprint("Doublet resolved")
    else:
        tio.tprint("Doublet not resolved")
    
main()
