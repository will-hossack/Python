import sys
import optics.wavelength as w

def main():


    wave = w.Default
    design = w.Design
    for l in sys.argv:
        if l.startswith("-w"):
            tokens = l.split("=")
            wave = float(tokens[1])
            w.setDefaultWavelength(wave)
        if l.startswith("-d"):
            tokens = l.split("=")
            design = float(tokens[1])
            w.setDesignWavelength(design)
        print(l)

    print("Wave is " + str(w.Default) + " and Design is " + str(w.Design))
main()
