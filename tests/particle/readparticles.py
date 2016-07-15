import tio as t
import vector as v
import particle as p

def main():
    file = t.openFile("Particle files",defaulttype = "data")

    system = p.ParticleSystem().readFile(file)

    t.tprint(system)


main()
