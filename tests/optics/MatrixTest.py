
import optics.matrix as mat
import matplotlib.pyplot as plt

def main():

    
    l = mat.ThickLensMatrix(0.02,1.5,10.0,-0.02)
    print(repr(l))

    nl = mat.ParaxialMatrix(l)
    print(repr(nl))
    

    #cavity = mat.CavityMatrix(-0.01,40,0.015)
    #print("Cavity : " + repr(cavity) + " with trace " + str(cavity.trace()))

    pg = mat.ParaxialGroup(100,l,10)

    print(repr(pg))


    planes = pg.planePair(-2.0)
    print(str(planes))

    cp = pg.cardinalPoints()
    print(str(cp))
    
    pt = pg.draw()
    plt.show(pt)

main()
