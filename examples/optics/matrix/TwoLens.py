"""
Example code to take a system of two thick singlet lenses the first on at 50 mm
and plot the  focal length of as the second lens is moved from 60 to 100 mm.

"""
import optics.matrix as mat
import matplotlib.pyplot as plt
import numpy as np

def main():

    fp = 50.0            # Poistion of front lens
    front = mat.ParaxialThickLens(fp,0.025,1.61,6.0,-0.035,5.0)
    print(repr(front))   # Info about front lens
    print("Focal length of front : " + str(front.backFocalLength()))
    
    bp = 60.0            # Start position of back lens
    back = mat.ParaxialThickLens(bp,-0.020,1.57,2.0,0.025,5.0)
    print(repr(back))    # Info about back lens
    print("Focal length of back : " + str(back.backFocalLength()))

    
    zpos = np.linspace(bp,100.0,30)     # np array of positions of back lens
    focal = np.zeros(zpos.size)         # np array to hold focal lengths
    for i,z in enumerate(zpos):         # Loop across zpos
        back.setInputPlane(z)           # Set postion of back lens
        system = front*back             # form sytsem
        focal[i] = system.backFocalLength() # calcuate back focal length

    plt.plot(zpos,focal)                # Do the plot with matplotlib
    plt.title("Focal Length Plot")
    plt.xlabel("Back lens position")
    plt.ylabel("Focal Length")
    plt.ylim(bottom = 0.0)
    plt.show()

main()
