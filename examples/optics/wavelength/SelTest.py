import optics.wavelength as w
import tio
import matplotlib.pyplot as plt

def main():



    l0 = tio.getFloat("Lambda_0",0.08)
    beta = tio.getFloat("Beta",1.25)
    index = w.Sellmeier(beta,l0)

    nd = index.getNd()
    vd = index.getVd()
    tio.tprint("Nd index : ",nd," Abbe No: ",vd)

        
    index.draw()
    plt.show()
    
main()
