import tio as t
import vector as c
import particle as p

def main():
    #
    ps = p.ParticleSystem().readFile()
    
    t.tprint("Initial total Linear Momentum is : ",ps.getLinearMomentum())
    t.tprint("Initial total Angular Momnetum is : ",ps.getAngularMomentum())
    t.tprint("Initial Total Energy is : ",ps.getKineticEnergy())


    if t.getBool("Elastic",True):
        ps[0].elasticCollision(ps[1])
    else:
        ps[0].inelasticCollision(ps[1])

    t.tprint("Left particle: ",ps[0])
    t.tprint("Right particle: ",ps[1])

    t.tprint("Final total Linear Momentum is : ",ps.getLinearMomentum())
    t.tprint("Final totalAngular Momnetum is : ",ps.getAngularMomentum())
    t.tprint("Final Total Energy is : ",ps.getKineticEnergy())


main()
