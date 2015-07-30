"""
Set to classes to implement a general integrator with Euclid, RK2, Modfied Euclid and RK4 integrators.
See integrator examples on how to use these for the particle.py classes.

Author: Will Hossack, The University of Edinburgh
"""
import math
import copy
#                         
#
class Derivative(list):
    """
    Class to hold a Derivative (need to be extended)
    """

    #
    #  
    def __init__(self,*args):
        """
        Constuctor with an arbirary set of arguments that are appended to the Derivative list
        """
        for dy in args:
            self.append(dy)
    
#            
#
class Point(list):
    """
    Class to hold the current state of the integrator. 
    """
    #
    #
    def __init__(self,*args):
        """
        Constructor with an arbitray set of arguments that are appended to the Point list
        """
        for y in args:
            self.append(y)
        
    #
    #
    def copy(self):
        """
        Method to make a deep copy of the current Point.
        """
        return copy.deepcopy(self)

    #         
    def eulerStep(self,delta,step):
        """
        Method to apply a EulerStep to the current point. 
        This assumes that the elements of Derivative and Point match, the Derivative obets the *float()
        operator and Point obey the += operator.
        param delta Derivatives to be applied
        param step the step size.
        """
        if len(self) != len(delta):
            raise IndexError("Point.eulerStep: Point and Derivative length missmatch")
        for i in range(len(self)):
            self[i] += delta[i]*step

    #        

    def weightedStep(self, deltas, w, step):
        """
        Method to apply a weightedStep to the current point.
        This assumes that the elements of Derivative and Point match, the Derivative obets the *float()
        operatote and Point obey the += operator.
        param delta Derivatives to be applied
        param w float[] being the Derivatives weightings
        param step the step size.
        """
        for i in range(len(deltas)):
            self.eulerStep(deltas[i],w[i]*step)

    #      
    #        
    def trialStep(self,delta,step):
        """
        Method to implement a trial step and form a new Point with the updated vaules.
        This assumes that Derivaties obey the *float() operaror and that Point obeys the a = self + b operator.
        param delta Derivatives to be applied
        param step float the step to be applied.
        """
        if len(self) != len(delta):
            raise IndexError("Point.trialStep: Point and Derivative length missmatch")
        p = Point()
        for i in range(len(self)):
            n = self[i] + delta[i]*step
            p.append(n)
        return p


#
#           Equations class used to define the equations. This must be supplied
#           my the user by extening this class.
            


class Equations:
    """
    Abstract class to define equations, this must be extended to implement the integration.
    """

    #
    #
    def __init__(self):
        """
        Constrcuctor, to be defined
        """
        self.point = Point()

    #      
    #
    def derivative(self,pt):
        """
        Method to calculate and return the derivatives at a specified Point.
        Needs to be defined.
        """
        return Derivative()
        
    #
    #
    def monitor(self):
        """
        Method that is called on each step of the integration to monitor progress and/or record
        output. Needs to be defined.
        """
        print("Monitor needs to be defined")

    #      
    #
    def terminate(self):
        """
        Method used to specify a termination of the integration, defaults to False.
        Optional, if not dfeined the inregrator will complete the number of specified steps.
        """
        return False

#
#       
#
class Integrator:
    """
    The underlying Integrator class to do the integration
    """

    
    def __init__(self,eqn,step,maxstep):
        """
        Constructor with same parameters as inheriting classes.
        param eqn the equation set to be inregrated
        param step the integration step (fixed)
        param maxstep the max number of integeations
        """
        self.eqn = eqn
        self.step = step
        self.maxstep = maxstep
        self.time = 0.0

    #
    #   
    #
    def forwardStep(self,pt,step):
        """
         The "forwards step" method, overwritten in the actual integrators.
        """
        print("Integrator.forwardStep needs to be defined")
        return Point()
    #
    #    
    #
    def run(self,step = None, maxstep = None):
        """
         Method to run (or re-run) the Integrator.
        param step step size (defaults to value set in constructor)
        param maxstep maximum number of steps (defaults to value set in constrcuor)
        return True if Equation.terminate() becomes True, False if maxstep exceeded.
        """

        #        Work out what has been sent
        if step != None:
            self.step = step
        if maxstep != None:
            self.maxstep = maxstep
        
        #        Number of time loop executed
        n = 0
        self.eqn.monitor()                  # Trigger initial monitor to record 
        while True:
            self.eqn.point = self.forwardStep(self.eqn.point,self.step)
            self.time += self.step
            self.eqn.monitor()              # Trigger monitor to record 
            n += 1                          
            #
            #        Check termination conditions
            if self.eqn.terminate():
                return True                 # Completed sucessfully
            if n > self.maxstep:
                return False                # Failed to terminate


#
#        
class Euler(Integrator):
    """
    Basic fixedstep Euler method, not really very useful,  but good for testing.
    """

    

    def __init__(self,eqn,step = 1.0 ,maxstep = float('Inf')):
        """
        Constructor 
        param eqn the equations to be solved
        param step the fixed step size (defaults to 1.0)
        param maxstep maximum number of steps, defaults in float('inf')
        """        
        Integrator.__init__(self,eqn,step,maxstep)

    #
    #
    def forwardStep(self,pt,step):
        """
        Internal forwardStep used to define the mothod.
        """
        der = self.eqn.derivative(pt)      # Get the derivative at point
        pt.eulerStep(der,self.step)        # Do a single Euler step
        return pt
        
        

class RungeKuttaTwo(Integrator):
    """
    The RK2 integrator with the derivatives calculated at the centre of the step interval. 
    This is symplectic so conserves evergty in dynamic problems and is ofter a good inital scheme. 
    """
    #
    #
    def __init__(self,eqn,step = 1.0 ,maxstep = float('Inf')):
        """
        Constructor 
        param eqn the equations to be solved
        param step the fixed step size (defaults to 1.0)
        param maxstep maximum number of steps, defaults in float('inf')
        """        
        Integrator.__init__(self,eqn,step,maxstep)
        


    def forwardStep(self,pt,step):
        """
        Internal forwardStep used to define the mothod.
        """
        der = self.eqn.derivative(pt)           # Get the derivative at point
        midPt = pt.trialStep(der,0.5*self.step) # trial half step
        midDer = self.eqn.derivative(midPt)     # Derivative at midpoint
        pt.eulerStep(midDer,self.step)          # Do Euler step
        return pt

class ImprovedEuler(Integrator):
    """
    Improved Euler with calcualtion of the derivative at the start and end of the step interval.
    Slighly longer clcualtion times that RK2, with typically similar results.
    """

    def __init__(self,eqn,step = 1.0 ,maxstep = float('Inf')):
        """
        Constructor 
        param eqn the equations to be solved
        param step the fixed step size (defaults to 1.0)
        param maxstep maximum number of steps, defaults in float('inf')
        """        
        Integrator.__init__(self,eqn,step,maxstep)
        self.weight = [0.5,0.5]                 # Point weights

    def forwardStep(self,pt,step):
        """
        Internal forwardStep used to define the mothod.
        """
        startDer = self.eqn.derivative(pt)       # Get the derivative at point
        endPt = pt.trialStep(startDer,self.step) # trial full step
        endDer = self.eqn.derivative(endPt)      # Derivative at end
        delta = [startDer,endDer]
        pt.weightedStep(delta,self.weight,self.step)
        return pt

class RungeKuttaFour(Integrator):
    """
    The full RK4 integartor that makes 4 estimates for the derivatives and uses them in a weighted
    sum for the final step. Much slower that the Euler / RK2 since it needs 3 trial points  at each
    step, but is much more accurate for long planet calcualtions or tricky SHM cases where the other
    diverge.
    """
    
    def __init__(self,eqn,step = 1.0 ,maxstep = float('Inf')):
        """
        Constructor 
        param eqn the equations to be solved
        param step the fixed step size (defaults to 1.0)
        param maxstep maximum number of steps, defaults in float('inf')
        """        
        Integrator.__init__(self,eqn,step,maxstep)
        sixth = 1.0/6.0
        third = 1.0/3.0
        self.weight = [sixth,third,third,sixth]

    def forwardStep(self,pt,step):
        """
        Internal forwardStep used to define the mothod.
        """ 
        k1 = self.eqn.derivative(pt)        # Derivative at start
        firstmidpt = pt.trialStep(k1,0.5*self.step)   
        k2 = self.eqn.derivative(firstmidpt) # First mid pt estimate
        secondmidpt = pt.trialStep(k2,0.5*self.step)
        k3 = self.eqn.derivative(secondmidpt) # second mid point estimate
        endpt = pt.trialStep(k3,self.step)
        k4 = self.eqn.derivative(endpt)

        delta = [k1,k2,k3,k4]

        pt.weightedStep(delta,self.weight,self.step)
        return pt

                           
