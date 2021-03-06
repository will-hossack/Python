from integrator import *
from integratorexamples import *
from vector import *
from particle import *
import matplotlib.pyplot as plt

#
#                  Example of forced damped SHM

def main():
    
    sp = Vector3d(-10,0,0)             # Start position
    sv = Vector3d(-5,0,0)               # start velocity
    mass = 2.0
    radius = 0.2
    particle = Particle(sp,sv,mass,0.0,radius)  # particle
    origin = Vector3d(0,0,0)         # Origin

    naturalFreq = float(input("Natrual Requency : "))
    spring = particle.mass * naturalFreq**2
    print("Spring Constant is : " + str(spring))

    damp = float(input("Damping : "))
    
    drivedir = Vector3d(1,0,0)
    driveamp = 0.1                    # drive amp
    
    omega = float(input("Drive freqency : "))
    simulationtime = 50.0*2.0*math.pi/naturalFreq

    eqn = ForcedDampedSHM(particle,origin,spring,damp,drivedir,driveamp,omega)

    step = simulationtime/10000
    integrator = RungeKuttaFour(eqn)
    integrator.run(step,10000)

    plt.plot(eqn.t,eqn.x,'r')
    plt.show()

main()

    
    
