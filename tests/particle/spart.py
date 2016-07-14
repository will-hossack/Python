import vector as v
import particle as p
import tio as t

def main():
    
    s = p.ParticleSystem("system")

    for i in range(10):
        pos = v.Vector3d(i,-i,2*i)
        vel = v.Vector3d(-i,i,2*i)
        part = p.Particle(pos,vel,mass = 2)
        s.append(part)


    t.tprint("Total Kinetic Engery is : " + str(s.getKineticEnergy()))
    t.tprint("Total Linear Momentum is : " + str(s.getLinearMomentum()))
    t.tprint("Centre of mass is : " + str(s.getCentreOfMass()))

main()
