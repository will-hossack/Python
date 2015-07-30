"""
Two example classes for the integrator class to implement 
a) projectile motion via the Projectile Class
b) force - dampled SHM with the ForceDampledSHM Class

Author: Will Hossack, The University of Edinburgh
"""
from integrator import *
from vector import *
from particle import *
#
#
class Projectile(Equations):
    """
    Particle projectile example with linear drag
    """
    
    #           
    def __init__(self,particle,damp):
        """
        Define constructor
        param particle Particle in initial state.
        param damp float damping (or drag) factor
        """
        self.particle = particle
        #
        #            Point of time, position and vectovity
        self.point = Point(self.particle.time,self.particle.position,\
                           self.particle.velocity)
        self.damp = damp
        #            List to return x and z locations
        self.x = []
        self.z = []

    #
    #
    def derivative(self,pt):
        """
        Method to get the derivative"
        """
        part = self.particle.copy()    # Take copy of local particle
        part.time = pt[0]              # Update time to current
        part.position = pt[1]          # Update position to current
        part.velocity = pt[2]          # Update velocity to current
        #
        #         Calculate force or particle
        force = part.getLinearDrag(self.damp) + part.getGravity()
        #
        #          Get acceleration (note this is Vector3d)
        acc = force/part.mass
        #
        #          Build the Derivative
        return Derivative(1.0,part.velocity,acc)

    #                
    #
    def monitor(self):
        """
        Method to record the x/z posutions of the particle for plotting
        """
        self.x.append(self.particle.position.x)
        self.z.append(self.particle.position.z)

    #                
    def terminate(self):
        """
         Method to terminate when z < 0
        """
        return self.particle.position.z < 0

#
#           
#
class ForcedDampedSHM(Equations):
    """
    Class for Forced damped SHM of a particle
    """
    #
    #
    def __init__(self,particle,origin,spring,damp,drivedir,driveamp,omega):
        """
        Define constrcutor, all needed
        param Particle the particle in inital state.
        param origin Vector3d where particle is attached by spring
        param spring float spring constant
        param damp float damping factor
        param drivedir Vector3d, directon of drive
        param driveamp float amplitude of drive
        param omega float angular frequency of drive
        """
        self.particle = particle
        self.origin = origin
        self.spring = spring
        self.damp = damp
        self.drivedir = drivedir
        self.driveamp = driveamp
        self.omega = omega
        #             Form point of time, position and velocity of particle
        self.point = Point(self.particle.time,self.particle.position,\
                           self.particle.velocity)
        self.t = []                  # data to record  (t and x)
        self.x = []


    def derivative(self,pt):
        """
        Method to for the derivative
        """
        part = self.particle.copy()    # Make copy of particle
        part.time = pt[0]              # Update time
        part.position = pt[1]          # Update poistion
        part.velocity = pt[2]          # Update velocity
        
        #
        #     Form force in two parts 
        dv = part.getLinearForce(self.origin,self.spring,self.damp)
        dv += part.getHarmonicDrivingForce(self.drivedir,self.driveamp,\
                                           self.omega)
        #     Form acceleration
        acc = dv/part.mass
        #
        #     Form deriavtive of step, veloicity, acceleration
        return Derivative(1.0,part.velocity,acc)


    def monitor(self):
        """
        Monitor time and x position or graph
        """
        self.t.append(self.particle.time.t)            # Get t position
        self.x.append(self.particle.position.x)        # Get x position


#
#
class Cyclotron(Equations):
    """
    Class for Cyclotron example for changed particle in constant B feild
    """
    
    def __init__(self,particle,bfield):
        """
        Constructor with:
        param partice Particle in initial state
        param bfeild Vector3d the constand B-field
        """
        self.particle = particle
        self.bfield = bfield
        #
        #        Form point of time, position and velocity
        self.point = Point(self.particle.time,self.particle.position,\
                           self.particle.velocity)
        self.x = []
        self.y = []
    
    #
    #
    def derivative(self,pt):
        """
        Form the derivatives.
        """
        part = self.particle.copy()     # Make copy of particle
        part.time = pt[0]               # Update time
        part.position = pt[1]           # Update position
        part.velocity = pt[2]           # Update velocity

        #
        #       Get magnetic force of particle
        f = part.getMagneticForce(self.bfield)
        #
        #       Apply accelerating e-field (bit of a hack)
        if abs(part.position.x) < 5:
            if part.velocity.x > 0:
                f += Vector3d(50,0,0)
            else:
                f += Vector3d(-50,0,0)
        #
        #       Form acceleration
        acc = f/part.mass
        #
        #       Derivative is timestep, velocity, acceleration
        return Derivative(1.0,part.velocity,acc)

    def monitor(self):
        """
        Monitor x/y posiition and print kinetice enery to terminal
        """
        self.x.append(self.particle.position.x)          # Get x position
        self.y.append(self.particle.position.y)          # Get y position
        print("Energy is : {0:6.3e}".format(self.particle.getKineticEnergy()))
