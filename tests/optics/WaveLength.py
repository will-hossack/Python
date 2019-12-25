import optics.wavelength as w
import tio as t

def main():

    wave = w.Default
    t.tprint("Default is : ",wave)

    w.setDefaultWavelength(0.65)
    wave = w.Default
    t.tprint("Updated default is : ",wave)
    

    design = w.Design
    t.tprint("Design is : ",design)

    w.setDesignWavelength(0.75)
    design = w.Design
    t.tprint("Updated design is : ",design)
    
main()
    
