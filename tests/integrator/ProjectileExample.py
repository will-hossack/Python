"""
Demo programme for the Integrator Projectile example
Plots of a trajectory of a specified particle under
damping launched from (0,0,0) with a angle in range 5 -> 85 degrees.            
"""
from integrator import *
from integratorexamples import *
from vector import *
from particle import *
import tio
import matplotlib.pyplot as plt
#
#             
def main():

    speed = tio.getFloat("Initial speed of particle",100.0,0.0)
    damp = tio.getFloat("Damping",3.0,0.0)

    mass = 5.0
    radius = 0.02


    for a in range(5,85,5):
        sp = Vector3d(0,0,0)                          # Start position
        sv = Vector3d().setPolarDegrees(speed,a,90) # Start velocity
        particle = Particle(sp,sv,mass,0.0,radius)    # Make particle

        #          Make the equations
        eqn = Projectile(particle,damp)

        #          Make the integrator
        integrator = ImprovedEuler(eqn)
        integrator.run(speed/1000)       # Run with step size of 
        print("Angle: {0:d} Range : {1:8.3f} Final speed : {2:8.3f} Energy : {3:8.3f}:".\
              format(a,abs(particle.position),abs(particle.velocity),particle.getKineticEnergy()))

        #          Get plot information
        plt.plot(eqn.x,eqn.z,label="Angle: {0:3.1f}".format(a))
    
    #                  Make plot nice
    plt.title("Projectile Motion with Drag of {0:4.2f}".format(damp))
    plt.legend(loc="upper right",fontsize="small")
    plt.xlabel("Range")
    plt.ylabel("Height")
    plt.ylim(0.0)
    plt.show()
    
   
main()
