"""
Set of classes to implement Paticle dynamics for meachanics and
Electrostatic problms.

Author: Will Hossack, The Univesrity of Edinburgh
"""

from vector import *
Gravity = Vector3d(0.0,0.0,-9.81)
#
#                 Classes to support particle simulation (mainly for
#                 planets and charged particle in fields)

#
#
class Time:
    """
    Class to hold time in an updatable way
    """
    # 
    def __init__(self,t=0.0):
       """
       Set time with initial optional value
       param t float initial time (defaults = 0.0)
       """
       self.t = t

    #
    #
    def __iadd__(self,dt):
        """
        Implment the += operator to add time or float(time)
        """
        if isinstance(dt,Time):
            self.t += dt.t
        else:
            self.t += dt
        return self
    
    #
    def __add__(self,dt):
        """
        Implement the newtime = time + dt opeator
        """ 
        if isinstance(dt,Time):
            nt = self.t + dt.t
        else:
            nt = self.t + dt
        return Time(nt)
    #
    #
    def __radd__(self,dt):
        """
        Implement the newtime = dt + time opeator
        """ 
        if isinstance(dt,Time):
            nt = self.t + dt.t
        else:
            nt = self.t + dt
        return Time(nt)


#
#                 Class to contain a impulse, being a force + time
#                 mainly used in dynamic simulatioms
#
class Impulse(Vector3d):
    #
    #             Constructor that takes
    #             x_or_v, x component, Vector3d or Impulse
    #             y_or_t, y component of t
    #             z component
    #             t component
    #
    def __init__(self,x_or_v,y_or_t = None ,z = None ,t = None):
        if isinstance(x_or_v,Impulse):
            Vector3d.__init__(self,x_or_v)
            self.t = x_or_v.t
        elif isinstance(x_or_v,Vector3d):
            Vector3d.__init__(self,x_or_v)
            if y_or_t == None:
                self.t = 1.0
            else:
                self.t = y_or_t 
        else:
            Vector3d.__init__(self,x_or_v,y_or_t,z)
            if t == None:
                self.t = 1.0
            else:
                self.t = t

    #
    #          return test string
    #
    def __str__(self):
        return "Impulse [{0:8.4e} , {1:8.4e}, {2:8.4e}] t: {3:8.4e}".\
            format(self.x,self.y,self.z,self.t)

    #       Return copy of current Impulse
    #       returns a copy of the current Impulse
    def copy(self):
        return Impulse(self)

    #       Get the force vector 
    #
    def getForce(self):
        return Vector3d(self.x,self.y,self.z)

    #        Redefine * to for intenration
    #
    def __mul__(self,b):
        t = self.t*b
        return Impulse(self.x,self.y,self.z,t)

#
#                 
class Particle(object):
    """
    Basic Particle class to hold information  about a single particle.
    """
   
    def __init__(self,position = Vector3d(),velocity = Vector3d(),mass = 0.0,charge = 0.0,radius = 0.0,time = 0.0,title = "Particle"):
        """
        Constuctor with
        param p Vector3d position (defaults to 0,0,0)
        param v Vector3d velocity (defaults to 0,0,0)
        param m mass (defaults to 0.0)
        param q change (defaults to 0.0)
        param radius (defaults to 0.0)
        param t time (defaults to 0.0)
        param title string name of particle (optional, defaults to "Particle")
        """
        self.position = Vector3d(position)         # Position
        self.velocity = Vector3d(velocity)         # Velocity
        self.mass = float(mass)                # Mass
        self.charge = float(charge)              # Charge
        self.radius = float(radius)              # radius (if needed)
        self.time = Time(time)
        self.title = title

    #
    #
    def __str__(self):
        """
        Implment str() method
        """
        return "Particle: Title {0:s} Mass: {1:6.3e} Charge: {2:6.3e} Radius: {3:6.3e}\n".\
            format(self.title, self.mass,self.charge, self.radius) +\
            "Position: " + str(self.position) + "\n" + "Velocity: " + str(self.velocity)

    #
    #           
    def copy(self):
        """
         Method to return a deep copy of the current Particle
        """
        return Particle(self.position,self.velocity,self.mass,\
                        self.charge,self.radius,self.time.t,self.title)
    #
    
    #
    def getKineticEnergy(self):
        """
        Get Kinetic Energy
        return the kinetic energy of the particle as float.
        """
        return 0.5*self.mass*self.velocity.absSquare()
    #
    #
    def setKineticEnergy(self,val):
        """
        Set Kinetic Energy by scaling velocity (direction not changed)
        param val float value that Kinetic Energy is set to
        """
        v = math.sqrt(2.0*val/self.mass)        # Required scalar velocity
        self.velocity.setLength(v)

    # 
    #
    def getLinearMomentum(self):
        """
        Get Linear Momentun
        return the linear momentum as a Vector3d.
        """
        return self.velocity * self.mass
    #
    #
    def getAngularMomentum(self, origin = Vector3d()):
        """
        Get angular momentum about origin
        return the angular moment as a Vector3d
        """
        r = self.position - origin
        return self.mass*r.cross(self.velocity)

    #           
    #
    def getGravitationalPotential(self,p):
        """
        Get gravitational potential at a specified position due to this particle.
        param p Vector3d position
        return gravitational potential at Vector3d p as a float
        """ 
        return  -self.mass/self.position.distance(p)

    
    def getElectrostaticPotential(self,p):
        """
        Get electrostatic potential as specified position due to this particle.
        param p Vector3d position where potential is calculated
        returns electrostatic potential as a float 
        """
        return  self.charge/self.position.distance(p)

    #
    #          
    def getGravitationalField(self,p):
        """
         Get Gravitational vector field due to this particle at specified position
        param p Vector3d position
        return gravitational field as Vector3d
        """
        r = p - self.position
        s = -self.mass/r.absCube()
        r *= s
        return r
    #
    #
    def getGravity(self):
        """
        Method to get Linear Gravitational Force, where Gravity is defined in the
        -z diretcion.
        """
        return Gravity*self.mass
    #
    #           
    def getElectrostaticField(self,p):
        """
        Get Electostatic vector field due to this paticle as specified position
        param  p Vector3d position
        param return electrostatic feild as Vector3d
        """
        r = p - self.position
        s = self.charge/r.absCube()
        r *= s
        return r

    #          Get the Electrostatic Potential Energy current
    #          Particle and specified particle
    #
    def getElectrostaticPotentialEnergy(self,p):
        """
        Get the Electrostatic Potential Energy between current Particle and specified particle.
        param p Particle, the specified Particle
        """
        if p == self:
            return 0.0
        else:
            v = self.getElectrostaticPotential(p.position)
            return p.charge*v

    #
    #
    def applyConstantElectrostaticField(self,ef,t):
        """
         Method to apply constant electrostatic field for specified time.
        param ef Vector3d constant electrostatic feild
        param t  time fields applield for.
        """
        accel = ef * (self.charge/self.mass)
        self.velocity += accel * t
        self.position += accel * (0.5*t*t)
        self.time += t

    #
    #          Method to get force linear from a spring on the current
    #          particle
    #          p point where spring is attacked
    #          k spring constant
    #          lam damping constant (defaults to 0.0)
    #          return force of particle to sping.
    #
    def getLinearForce(self,p,k,lam = 0.0):
        f = p - self.position
        al,fn = f.absNormalised()
        kx = k*al        # Required stength of force
        if lam!= 0.0:    # Dampint term
            vl = lam*fn.dot(self.velocity)
            kx -= vl
        fn *= kx
        return fn
    #
    #        Method to get a*cos(omega t)
    #        dir direction attractive point for force as Vector3d
    #        a amplitude of driving force
    #        omega angular frequeny of driving force
    #        t time 
    #        phase shift of force (default = 0.0)
    #        return driving force as Vector3d
    #
    def getHarmonicDrivingForce(self,dir,a,omega,phase = 0.0):
        f = dir.copy()
        amp = a*math.cos(omega*self.time.t + phase)
        f.setLength(amp)
        return f
    #
    #        Method to get Rayleigh turbulent drag which opposed the
    #        velocity and is given 0.5*rpo*cd*A*v**2
    #        where A is cross-sectional area.
    #        rpo mass density of fluid 
    #        cd drag coefficeint (default 0.47 for sphere) 
    #
    def getRayleighDrag(self,rpo,cd = 0.47):
        d = 0.5*rpo*cd*math.pi*self.radius*self.radius*\
            self.velocity.absSquare()
        dv = self.velocity.copy()
        dv.setLength(-d)
        return dv
    
    #        Method to get linear drag which opposed the
    #        velocity and is given by 6pi x mu r |v|
    #
    def getLinearDrag(self,mu):
        d = 6.0*math.pi*mu*self.radius*self.velocity.abs()
        dv = self.velocity.copy()
        dv.setLength(-d)
        return dv
        

    #       Method to apply an impulse to the current Particle
    #
    def applyImpulse(self,imp):
        self.position += self.velocity*imp.t     # Update position
        dv = imp.getForce()*(imp.t/self.mass)
        self.velocity += dv
        self.time += imp.t


    #       Method to get force due to Magnetic field
    #
    def getMagneticForce(self,bfield):
        f = self.velocity.cross(bfield)
        f *= self.charge
        return f

    #       Define += to apply impulse (for integration)
    #
    def __iadd__(self,imp):
        self.applyImpulse(imp)
        return self

    #     + to apply impulse and return new updated particle
    #
    def __add__(self,imp):
        p = self.copy()                         # Make copy
        p.applyImpulse(imp)
        return p


    def inelasticCollision(self, b):
        """
        Perform an inelastic collision between self and partile b.
        Afetr collision both stick togeter so have same vectocity
        """
        pt = self.getLinearMomentum() + b.getLinearMomentum()    # Total linear momentum
        mt = self.mass + b.mass          # Total mass
        v = pt/mt                        # New common vecocity
        self.velocity = v.copy()
        b.velocity = v.copy()


    def elasticCollision(self , b):
        """
        Perform an elastical collision between salf and particle b.
        After collision both particles will have new velocities.
        """
        mt = self.mass + b.mass              # Total mass
        vs = (self.mass - b.mass)*self.velocity + 2.0*b.getLinearMomentum()
        vs /= mt
        vb = (b.mass - self.mass)*b.velocity + 2.0*self.getLinearMomentum()
        vb /= mt
        self.velocity = vs.copy()
        b.velocity = vb.copy()

        


#              Class to hold a list of particles
#
class ParticleSystem(list):
    #
    #           Constructor to take partiels and add to list
    #
    def __init__(self,title = "Particles" ,*args):
        self.title = title
        for p in args:
            self.append(p)
    #
    def readFile(self,file):
        """
        Method to read in a particle system from a text file
        """
        
        if isinstance(file,str):
            try:
                file = open(file,"r")
            except:
                raise IOError("ParticleSystem.readFile: failde to open : " + filename)

        #            Read through file line at a time
        for line in file.readlines():

            if not line.startswith("#") and len(line) > 0:
                tokens = line.split()
                if tokens[0].startswith("title"):
                    self.title = str(tokens[1])
                else:
                    pos = Vector3d()
                    vec = Vector3d()
                    mass = 0.0
                    charge = 0.0
                    radius = 0.0
                    for i in range(0,len(tokens),2):

                        if tokens[i].startswith("pos"):
                            pos = Vector3d(eval(tokens[i+1]))
                        elif tokens[i].startswith("vel"):
                            vec = Vector3d(eval(tokens[i+1]))
                        elif tokens[i].startswith("mass"):
                            mass = float(eval(tokens[i+1]))
                        elif tokens[i].startswith("charge"):
                            charge = float(eval(tokens[i+1]))
                        elif tokens[i].startswith("radius"):
                            radius = float(eval(tokens[i+1]))
                        else:
                            print("Unknown token " + tokens[i])
                    part = Particle(pos,vec,mass = mass, charge = charge, radius = radius)
                    self.append(part)

        return self

    #
    #            Method to return while list as a string
    #
    def __str__(self):
        s = "Particle System: {0:s} Number of particles: {1:d}".format(self.title,len(self))
        for pt in self:
            s += "\n"                    # Add new line
            s += str(pt) 

        return s

    #            Method to get total kinetic energy
    #
    def getKineticEnergy(self):
        k = 0.0
        for p in self:
            k += p.getKineticEnergy()
        return k

    #           Method to get total linear momentum
    #
    def getLinearMomentum(self):
        m = Vector3d()
        for p in self:
            m += p.getLinearMomentum()
        return m


    def getCentreOfMass(self):
        """
        Get the centre of mass of the system
        """
        cm = Vector3d()
        tm = 0.0
        for p in self:
            cm += p.position*p.mass
            tm += p.mass

        cm /= tm
        return cm

    #          Method to get total angular momentum of the system about a specivied origin
    #
    def getAngularMomentum(self, origin = Vector3d()):
        m = Vector3d()
        for p in self:
            m += p.getAngularMomentum(origin)
        return m

    #           Get gravitational potential at a specified position due to system of  particle.
    #           p Vector3d position
    #           return gravitational potential as a float 
    #
    def getGravitationalPotential(self,p):
        g = 0.0
        for pt in self:
            g += pt.gerGravitationalPotential(p)

        return g

    #           Get electrostatic potential as specified position due to this system of particles
    #           p Vector3d position
    #           returns electrostatic potential as a float 
    #
    def getElectrostaticPotential(self,p):
        v = 0.0
        for pt in self:
            v += pt.getElectrostaticPotential(p)

        return v

    #
    #           Get Gravitational field due to this particle system at specified position
    #           p Vector3d position
    #           return gravitational feilds at Vector3d
    #
    def getGravitationalField(self,p):
        g = Vector3d()
        for pt in self:
            g += pt.getGravitationalField(p)
 
        return g
    #
    #           Get Electostatic feild due to this paticle system as 
    #           specified position
    #
    def getElectrostaticField(self,p):
        e = Vector3d()
        for pt in self:
            e += pt.getElectrostaticField(p)

        return e
    
    #
    #         Get the total electrostatic potential energy in the
    #         configuration. 
    #         return total electrostatic potential energy as a float.
    #
    def getElectrostaticPotentialEnergy(self):
        u = 0.0
        for i in range(len(self)):
            for j in range(i + 1, len(self)):
                u += self[i].getElectrostaticPotentialEnergy(self[j])

        return u


