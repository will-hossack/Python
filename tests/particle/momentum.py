import particle as p
import vector as v
import tio as t

def main():
    vel = t.getVector3d("Velocity :")
    pos = t.getVector3d("Position :",[0,0,0])
    # org = t.getVector3d("Origin :",[0,0,0])

    left = p.Particle(pos,vel,mass=5,title="Left")
    
    right = p.Particle(pos,-3.*vel,mass=10,title="Right")
    t.tprint(left)
    t.tprint(right)


    mom = left.getLinearMomentum() + right.getLinearMomentum()
    k = left.getKineticEnergy() + right.getKineticEnergy()

    t.tprint("Inital linear momentum is : ",mom)
    t.tprint("Initial Energy is : ",k)

    left.elasticCollision(right)

    t.tprint(left)
    t.tprint(right)

    mom = left.getLinearMomentum() + right.getLinearMomentum()
    k = left.getKineticEnergy() + right.getKineticEnergy()

    t.tprint("Final linear momentum is : ",mom)
    t.tprint("Final Energy is : ",k)

    

main()
