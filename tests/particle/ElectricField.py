import tio as t
import particle as p
import vector as v


def main():
    
    ps = p.ParticleSystem().readFile(t.openFile("File",defaulttype="part"))
    t.tprint(ps)

    while True:
        pos = t.getVector3d("Position")
        ef = ps.getElectrostaticPotential(pos)
        t.tprint("Field at : ",pos," is : ",ef)

main()
