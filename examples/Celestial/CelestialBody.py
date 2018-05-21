"""
Classes for Celelstial body similation.
"""
import vector as v
import tio as t
import math

G = 6.674e-11 

class CelestialBody(object):
    """
    The CelestialBody class to hold a general cestial body
    """
    def __init__(self,m = 1.0 ,pos = [0,0,0], vel = [1,0,0],name = "Planet"):
        """ Constrcutor"""

        self.mass = float(m)
        self.position = v.Vector3d(pos)
        self.velocity = v.Vector3d(vel)
        self.name = name

    def __str__(self):
        """ to str method to monitor progress """
        
        return "m = {0:8.4e} p = {1:s} v = {2:s} n = {3:s}".\
            format(self.mass,str(self.position),str(self.velocity),self.name)

    def kineticEnergy(self):
        """   Kinetic energy of body """
        return 0.5*self.mass*self.velocity.absSquare()

    def potentialEnergy(self,body):
        """    Potential energy between self and another body or list of bodies
        """
        if isinstance(body,list):
            pe = 0.0
            for b in body:
                if b != self:
                    pe += self.potentialEnergy(b)
                return pe
        else:
            return -G*self.mass*body.mass/self.position.distance(body.position)

    def gravitationalForce(self,body):
        """     Force due to gravity on self due to other body or list of bodies
        """
        if isinstance(body,list):
            f = v.Vector3d()    # initialise force
            for b in body:    # 
                if b != self:
                    f += self.gravitationalForce(b)
                    
            return f
        else:                 # single body
            r = body.position - self.position
            return G*self.mass*body.mass/r.absCube()*r

        

class CelestialSystem(list):
    """      List of Clestical Objects for the simulation"""

    def __init__(self,*args):
        """ Constructor to append CelestialBodies
        """
        list.__init__(self)
        for body in args:
            self.append(body)
        self.t = 0.0

    def __str__(self):
        """    Print out as string """
        s = ""
        for body in self:
            s += str(body) + "\n"
        return s

    def kineticEnergy(self):
        """  Total kinetic energy is system
        """
        ke = 0.0
        for b in self:
            ke += b.kineticEnergy()
        return ke

    
    def accelerations(self):
        """    Get the acceletaions as a list """
        acc = []
        for b in self:
            a = b.gravitationalForce(self)/b.mass
            acc.append(a)
        return acc

    def initialise(self):
        """    Initialise the system bu setting the current acceleations
        """
        self.previousaccel = self.accelerations()

    def update(self, dt):
        """     Update timestep dt
        """
        acc = self.accelerations()
        #    Update position
        for i in range(len(self)):
            self[i].position += self[i].velocity*dt + (4.0*acc[i] - self.previousaccel[i])*dt*dt/6.0
        #     Get new acceleration
        newacc = self.accelerations()
        #    Update velovity
        for i in range(len(self)):
            self[i].velocity += (2.0*newacc[i] + 5.0*acc[i] - self.previousaccel[i])*dt/6.0

        #    Set preious accel
        self.previousaccel = acc
        self.t += dt
        return self.t
    
def main():
    sun = CelestialBody(1.989e30,[0,0,0],[0,0,0],"Sun")
    rad = 149.6e9
    vel = math.sqrt(G*sun.mass/rad)
    earth = CelestialBody(5.972e24,[rad,0,0],[0,vel,0],"Earth")

    solar = CelestialSystem(sun,earth)
    t.tprint(solar)

    t.tprint("Kinetic energy of earth is ", solar[1].kineticEnergy())
    t.tprint("Potental energy due to sun is ",solar[1].potentialEnergy(solar[0]))
    t.tprint("Gravitational force on earth is ",solar[1].gravitationalForce(solar))

    solar.initialise()
    while(solar.update(100.0) < 1000.):
        t.tprint(solar)

    
main()
