""" 
Program to read in wavelength , index values from CSV file and 
fit first order Sellmier 

This used the tio module for input/output
"""

import csvfile as csv
from optics.wavelength import Sellmeier,Sodium_D
import tio
import matplotlib.pyplot as plt

def main():

    #      Read in csv filename
    file = tio.getFilename("File","csv")   # Open File
    wave,ref = csv.readCSV(file)           # Read to two np arrays

    #    Create a Sellmier index and it it to the read in np arrays
    index = Sellmeier().fitIndex(wave,ref)
       
    #     Print out results and Nd / Vd values                   
    tio.tprint("Index : " , repr(index))
    tio.tprint("Nd : ",index.getNd()," Vd : ",index.getVd())
    
    #     Get the digital rerivative at the Sodium Doublet
    dn = index.getDerivative(Sodium_D)
    tio.tprint("Derivative is : ", dn)
    
    #     Plot out the data valuse a s x
    plt.plot(wave,ref,"x")
    index.draw()       # Defaul draw of the index from 0.35 -> 0.7
    plt.title("Plot of Sellmeier index")
    plt.xlabel("Wavelength in microns")
    plt.ylabel("Index")
    plt.show()

    

    

    

    
main()

