"""
Set of classes to implement various types lenses, being lists of surfaces, with additional
method to extract the geomertic parameters in a simple way.
"""
import optics.surface as sur
import optics.ray
import optics.matrix as matrix
from vector import Vector3d,Unit3d
import optics.wavelength as wl
import optics.analysis as ana
import tio
from matplotlib.pyplot import plot


#      Define a current lens as a Global that defaults to a singlet.


def setCurrentLens(lens):
    """
    Function to set the Current Lens, it is held in global CurrentLens

    :param lens: the lens or a filename, if it is a str, it will try and open the DataBaseLens of this type.
    :type lens: Lens or str
    """
    global CurrentLens
    if isinstance(lens,str):
        CurrentLens = DataBaseLens(lens)
    else:
        CurrentLens = lens

def getCurrentLens():
    """ 
    Function to get the current default lens, if this is not
    overriddent it is initally set to the default SimpleSinglet

    :return: Current default lens
    """
    return CurrentLens

#   Global Current Angle (mainly used by GUI)
CurrentAngle = Unit3d(0.0,0.0,1.0)

def getCurrentAngle():
    return CurrentAngle

def setCurrentAngle(u):
    global CurrentAngle
    CurrentAngle = Unit3d(u)


#
class OpticalGroup(list):
    """
    Class to hold a list of surfaces in order they will be encourtered by a ray.
    this is the class typically used to represent a lens or optical system for for full ray tracing.

    The Group has a group point that defines the start of the group in global coordinates
    and the surfaces are located relative to this point. Moving the group point will
    move the whole group on-mass.

    :param group_pt: The Group reference point in global coordinates, (Defaults = 0,0,0)
    :type grout_pt: Vector3d or float
    :param \*args: list of Surfaces to be added to the OpticalGroup in order
    

    """

    def __init__(self,group_pt = Vector3d(),*args):
        """
        Constrtructor to make an OpticalGroup and add any Surfaces supplied in the argument list.
       
        """
        list.__init__(self)
        
        self.title = "Optical Group"        # Identify title (for user identification)
        self.aperture = None                # Aperture when added
        self.iris = None                    # Variable iris aperture when added
        self.paraxial = None                # associated paraxial group (auto added)
        self.wavelength = wl.Design         # Default wavelength 
        self.group = None                   # Allow to be member of an other group
        self.setPoint(group_pt)             # Us method that does tidy up as well
        #
        #              Add the surfaces in order.
        for s in args:
            self.add(s)


    def __str__(self):
        """
        Returm basic information by str() typically overwritten in extending classes.
        """
        return "pg: {0:s} t: {1:s} sn : {2:s}".format(str(self.point),self.title,str(len(self))) 
            
    #
    def __repr__(self):
        """
        Implement repr()
        """
        return "{0:s}: ".format(self.__class__.__name__) + str(self)
    
        

    def getInfo(self):
        """
        Get a full in for of the group + all its component surfaces.
        
        :return: Full detauils as formatted string.

        """
        st  = repr(self)
        for s in self:
            st += "\n{0:s}".format(repr(s))

        return st
    
    
    def setPoint(self,pt):
        """
        Method to set/reset the group point either with float, int or ray.Position

        :param pt: Set / reset the group point
        :type pt: Vector3d or float

        """
        if isinstance(pt,float) or isinstance(pt,int):
            self.point = Vector3d(0.0,0.0,float(pt))
        else:
            self.point = Vector3d(pt)
        self.paraxial = None          # Remove paraxial matrix since geometery changed.

    def getPoint(self):
        """
        Get the local point
        """
        return self.point

    def movePoint(self,delta):
        """
        Move the reference point by about delta

        :param delta: distance moved
        :type delta: Vector3d or float
        """
        if isinstance(delta,float) or insinstance(delta,int):
            self.point += Vector3d(0.0,0.0,float(delta))
        else:
            self.point += delta
        self.paraxial = None

    #
    def add(self,surface):
        """
        Method to add a Surface to the end of the OpticalGroup, it also updates the 
        .group variable is the Surface.

        :param surface: Surface to be appended.
        :type surface: OpticalSurface or 'OpticalGroup'

        Use this method rather than the underlying append since this it updates internal variables.

        If an OpticalGroup in added the induvudual surafecs are added, see makeStandAlone() for deatils

        """

        if isinstance(surface,OpticalGroup):        # Unpack OpticalGroup
            sl = list(surface)                      # Make local list copy
            for s in sl:
                self.add(s.makeStandAlone())
            
        else:

            self.append(surface)
            self.paraxial = None                         # Clear any old paraxial matrix.
            surface.group = self                         # Make surface joint group
            if isinstance(surface,sur.CircularAperture): # Record location of aperture
                self.aperture = surface

            if isinstance(surface,sur.IrisAperture):     # Record location of iris
                self.iris = surface

        return self

    #
    def scale(self,a):
        """
        Scale the whole group, scales all surfaces but NOT the group point, normally
        used to set the focal length.

        :param a: the scaling factor.
        :type a: float

        """
        for s in self:
            s.scale(a)
 
        self.paraxial = None
        return self
        
    #      
    def planePair(self,mag, xsize = 100.0, ysize = None , wave = wl.Design):
        """
        Get the object / image paraxial plane pair for a specified imaging magnification and object plane size.
        The underlying calcualtion uses paraxial matrix formultion to locate the planes.

        :param mag: the magification (shoukld be -ve for imaging system)
        :type mag: float
        :param xsize: horizontal size of the object plane (Default = 100mm)
        :type xsize: float
        :param ysize: vertical size of the object plane (Default = xsize)
        :type ysize: float
        :param wave: float wavelength (default = wl.Design)
        :type wave: float
        :return: object and image ImagePlane as a two element list.

        """
        if ysize == None:
            ysize = xsize
        pt = self[0].getPoint()
        pm = self.paraxialGroup(wave)
        pobj,pima = pm.planePair(1.0,mag)
        obj = sur.ImagePlane([pt.x,pt.y,pobj.inputPlane()],xsize,ysize)
        ima = sur.ImagePlane([pt.x,pt.y,pima.inputPlane()],xsize*abs(mag),ysize*abs(mag))
        return obj,ima

    def imagePoint(self,op,wave = wl.Design):
        """
        Method to give three-dimensional image of a point in object space in global coordinates using the ideal paxial
        formulas.

        :param op: object point, can be a Vector3d or Angle; if Angle it will assume an onfinite object.
        :type op: Vector3d or Angle
        :param wave: Wavelnegth used (Default = wl.Default)
        :type wave: float
        :return: Vector3d the locattion of the image point.
        """
        pm = self.paraxialGroup(wave)
        p =  pm.imagePoint(op)
        return Vector3d(p.x + self.point.x, p.y + self.point.y, p.z)

    
    def entranceAperture(self):
        """
        Get the entrance aperture, begin circular aperture at the edge of the first element.

        :return: CircularAperture with reference point in global coordinates.

        """
        if isinstance(self[0],list):
            s = self[0][0]
        else:
            s = self[0]

        pt = s.getPoint()                     # Reference pt of first surface in global coordinates
        if isinstance(s,sur.QuadricSurface):  # If curved surface, find the edge plane and add/subsract it
            pt.z += s.edgePlane()

        return sur.CircularAperture(pt,s.maxRadius)

    def exitAperture(self):
        """
        Get the exit aperture, being a circular aperture at the edge of the last element.

        :return: CircularAperture with reference point in globalcoordinates.

        """
        if isinstance(self[-1],list):
            s = self[-1][-1]
        else:
            s = self[-1]

        pt = s.getPoint()                      # Reference point of last surafce in global coordinates
        if isinstance(s,sur.QuadricSurface):   # If curved surafce, find surface edge plane and add/subtract
            pt.z +=  s.edgePlane()

        return sur.CircularAperture(pt,s.maxRadius)



    def paraxialMatrix(self,wave = wl.Design, first = 0, last = -1):
        """
        Get a the paraxial matrix for a subset of surfaces.

        :param wave: the wavelngth (Default = optics.wavelength.Design)
        :type wave: float
        :param first: first surface (Default = 0)
        :type first: int
        :param last: last surface (Default = -1)
        :type last: int

        """
        mat = matrix.ParaxialMatrix()           # Start with unit matrix
        if first >= len(self):
            raise IndexError("lens.OpticalGroup.paraxialMatrix, index out of range")
        
        if first == 0 :
            z = 0
        else:
            z = self[first].point.z               # Start point
        nl = wl.AirIndex().getValue(wave) 
        if first > 0:
            s = self[first-1]
            if s.type == 1:
                nl = s.refractiveindex.getValue(wave)
        
        if last == -1:
            last = len(self)

        for s in self[first:last]:
            if isinstance(s,sur.QuadricSurface) and s.type != 0:  # Curved surface and not clear
                zp = s.point.z
                mat += (zp - z)                                # propagate to surface
                if s.type == 1:
                    nr = s.refractiveindex.getValue(wave)         # index on right
                else:
                    nr = -nl
                m = matrix.DielectricMatrix(nl,nr,s.curvature)
                mat *= m                                       # add surface
                z = zp
                nl = nr
        #     
        return mat
        
            
        

    def paraxialGroup(self,wave = wl.Design):
        """
        Get the paraxial group pf this group at spectified wavelength. This will be
        remade if either first call, surface added, wavelength change or scale change.

        :param wave: the wavelength, if None it will use the defaults wavelngth for the Group
        :type wave: float
        :return: the :class:`optics.matrix.ParaxialGroup`

        """
        #              Make paraxial matrix if needed.
        #
        if self.paraxial == None or wave != self.wavelength: # make a new matrix
            self.wavelength = wave
            mat = self.paraxialMatrix(wave)
            en = self.entranceAperture()               # Get input/output heights
            ex = self.exitAperture()
            #  Make Paraxial Group with matching ref point, blank matrix cortect input/output heights
            self.paraxial = matrix.ParaxialGroup(self.point.z,mat,en.maxRadius,ex.maxRadius)
            
        
        return self.paraxial
            

    def draw(self):
        """
        Method to draw the surfaces   (but NOT the paraxial planes, see Lens.draw below for more useful method)
        """
        for s in self:
            s.draw()


#
class Lens(OpticalGroup):
    """
    Class to expend OpticalGroup with extra methods that assume that the Group hold a compound lens.
    This class and the expending classes Singlet / SimpleSinglet are the two main user classes for ray tracing.

    This class holds a list of surfaces in order they will be encourtered by a ray.
    this is the class typically used to represent a lens for for full ray tracing.

    The Lens has a group point that defines the start of the group in global coordinates
    and the surfaces are located relative to this point. Moving the group point will
    move the whole group on-mass.


    :param group_pt: the group point, Defauts = (0,0,0)
    :type group_pt: Vector3d or float
    :param \*args: list of OpticalSurfaces
    
    """
    
    def __init__(self,group_pt,*args):
        """
        Constucructor to make a Lens and add any Surfaces supplied in the argument list.
        """
        OpticalGroup.__init__(self,group_pt,*args)
    

    def cardinalPoints(self,wave = wl.Design):
        """
        Method to get the six cardinal point of the lens system in global coordinates as a list of Vector3d. 
        The z componts come from ParaxialGroup.cardinalPoints() while the x/y componets are x/y components of the OpticalGroup.
        
        :param wave: the wavelength, (Default = wl.Design)
        :type wave: float
        :return: list of six Vector3d as specified below,.

        Order is

        * 0: Front Focal Point
        * 1: Back Focal Point
        * 2: Front principal plane
        * 3: Back principal plane
        * 4: Front nodal point (front principal plane for system in air)
        * 5: Back nodal point (back principal plane for system in air)

        """
        pg = self.paraxialGroup(wave)
        card = pg.cardinalPoints()        # Get candinal point as a float list from ParaxialGroup
        cardinal = []
        for c in card:                    # Return them into a list of Vector3d by adding x/y from group point
            p = Vector3d(self.point.x,self.point.y,c)
            cardinal.append(p)

        return cardinal                   # Return list of Positions


    def entrancePupil(self,wave = wl.Design):
        """
        Get the entrance pupil, being and image of the aperture
        """
        if self.aperture == None:
            return self.entranceAperture()

        ia = self.index(self.aperture)
        matrix = self.paraxialMatrix(wave,0,ia)
        matrix.inverse()
        p = self.aperture.getPoint()                # Position of aperture in global
        p.z += matrix.thickness - matrix.B/matrix.D # Location 
        mag = abs(matrix.A - matrix.C*matrix.B/matrix.D) # Mag

        return sur.CircularAperture(p,mag*self.aperture.maxRadius)

        

    def exitPupil(self,wave = wl.Design):
        """
        Get the exit pupil, being the image of the apreture
        """ 
        if self.aperture == None:
            return self.exitAperture()

        ia = self.index(self.aperture)
        matrix = self.paraxialMatrix(wave,ia,-1)
        p = self.aperture.getPoint()                # Position of aperture in global
        p.z += matrix.thickness - matrix.B/matrix.D # Location 
        mag = abs(matrix.A - matrix.C*matrix.B/matrix.D) # Mag

        return sur.CircularAperture(p,mag*self.aperture.getRadius())

        


    def setIris(self,ratio):
        """
        Set the iris ratio if there is an IrisApeture in the group; if there is no iris this is ignored without message.

        :param ratio: Iris ratio between 0 -> 1
        :type ratio: float.

        """
        if self.iris != None:
            self.iris.ratio = ratio
        return self


    def backFocalLength(self, wave = wl.Design):
        """
        Method to get the back focal length calculated by paraxial matrix methods, which for positive lens will be +ve.
        
        :param wave:  specified wavelength (default is wl.Design).
        :type wave: float
        :return: the focal length.

        """
        pm = self.paraxialGroup(wave)
        return pm.backFocalLength()

    def frontFocalLength(self, wave = wl.Design):
        """
        Method to get the front focal length calculated by paraxial matrix mecthods, which for positive lens will be +ve.
        
        :param wave:  specified wavelength (default is wl.Default).
        :type wave: float
        :return: the focal length.

        """
        pm = self.paraxialGroup(wave)
        return pm.backFocalLength()

    
        
    def setFocalLength(self,f,wave = wl.Design):
        """
        Method to set the geometric focal length of the lens by scaling. This assumes that the lens is
        in air and we are setting the 'back focal length'. If requested focal length is -ve then
        the surface curvatures will be reversed  but the positions will not be altered. This may give odd results
        for compound lenses, use with care.

        :param f: the target focal length
        :type f: float
        :param wave: the wavelength (Default = wl.Design)
        :type wave: float

        """
        fc = self.backFocalLength(wave)
        self.scale(f/fc)
        return self


    def frontFocalPlane(self, wave = wl.Design):
        """
        Get the Front Focal Plane as an surface.OpticalPlane is global coordinates.

        :param wave: the wavelength (Default = optics.wavelength.Default)
        :type wave: float
        :return: surface.OpticalPlane in global coordinates.
        """
        pm = self.paraxialGroup(wave)
        zplane = pm.frontFocalPlane()
        return sur.OpticalPlane(Vector3d(self.point.x,self.point.y,zplane))

    def backFocalPlane(self, wave = wl.Design):
        """
        Get the back Focal Plane as an surface.OpticalPlane is global coordinates.

        :param wave: the wavelength (Default = wl.Design)
        :type wave: float

        return surface.OpticalPlane in global coordinates.
        """
        pm = self.paraxialGroup(wave)
        zplane = pm.backFocalPlane()
        return sur.OpticalPlane(Vector3d(self.point.x,self.point.y,zplane))


    def frontPrincipalPlane(self,wave = wl.Design):
        """
        Get the front principal place as an surface.OpticalPlane in gobal cordilanes

        :param wave: the wavelength (Default = wl.Design)
        :type wave: float
        :return: :py:class:`optics.surface.OpticalPlane` in Global coordinates

        """
        pm = self.paraxialGroup(wave)
        zplane = pm.frontPrincipalPlane()
        return sur.OpticalPlane(Vector3d(self.point.x,self.point.y,zplane))


    def backPrincipalPlane(self,wave = wl.Design):
        """
        Get the back principal place in gobal cordilanes as a surface.OpticalPlane in global coordinates

        :param wave: the wavelength (Default = wl.Design)
        :type wave: float
        :return: :py:class:`optics.surface.OpticalPlane` in global coordinates

        """
        pm = self.paraxialGroup(wave)
        zplane = pm.backPrincipalPlane()
        return sur.OpticalPlane(Vector3d(self.point.x,self.point.y,zplane))

    
    def frontNodalPoint(self,wave = wl.Design):
        """
        Get the front nodal point as ray.Position in gobal cordilanes

        :param wave: the wavelength (Default = wl.Design)
        :type wave: float
        :return: :py:class:`Vector3d` in Global coordinates

        """
        pm = self.paraxialGroup(wave)
        zplane = pm.frontNodalPoint()
        return Vector3d(self.point.x,self.point.y,zplane)

    def backNodalPoint(self,wave = wl.Design):
        """
        Get the front nodal point as ray.Position in gobal cordilanes

        :param wave: the wavelength (Default = optics.wavelength.Default)
        :type wave: float
        :return: Vector3d  in Global coordinates

        """
        pm = self.paraxialGroup(wave)
        zplane = pm.backNodalPoint()
        return Vector3d(self.point.x,self.point.y,zplane)


    def petzvalSum(self,wave = wl.Design):
        """
        Calcualte the Perzval field curvature sum assuming air on the left side.

        :return: Perzal sum as a float.
        """
        leftindex = wl.AirIndex()
        ps = 0.0
        for s in self:
            if s.type == 1:           # Refrating
                rightindex = s.refractiveindex
                nl = leftindex.getValue(wave)
                nr = rightindex.getValue(wave)
                ps += s.curvature*(nr - nl)/(nr*nl)
                leftindex = rightindex
        return ps

    

    
    def draw(self,planes = True, legend = False):
        """
        Method to draw the surfaces and add the paraxial planes is requested.

        :param planes: draw the paraxial planes (Default = True)
        :type planes: bool
        :param legend: draw the legend box (Default = False)
        """
        for s in self:
            s.draw()
        if planes:
            if self.paraxial == None:
                pg = self.paraxialGroup()
            self.paraxial.draw(legend)


class Singlet(Lens):
    """
    Class to implement a singlet lens with simpler interface than OpticalGroup or lens and additional methods to alter the 
    lens surfaces. The default is a 10mm bicovex lens of glass BK7 with focal length of 97mm.

    :param pt_or_z: the group point (Default =( 0.0,0.0,0.0))
    :type pt_or_z: Vector3d or float
    :param c1: curvature of front surface (default = 0.01)
    :type c1: float
    :param t: thickness at centre (Default = 5.0)
    :type t: float
    :param cr: curvature of back surface (Default = -0.01)
    :type cr: float
    :param radius: max radius (Defaults = 10.0)
    :param index: refractive index, may be RefrativeIndex or material key, (default = "BK7")
    :type index: str or RefractiveIndex

    """
    def __init__(self, pt_or_z = 0.0, cl = 0.01 , t= 5.0 ,cr = -0.01 ,rad = 10.0,index= "BK7"):
        """
        Basic constructor for Singlet
        
        """
        Lens.__init__(self,pt_or_z)
        
        if isinstance(index,str):                        # Look index is given string key
            index = wl.MaterialIndex(index)

        sl = sur.SphericalSurface(0.0,cl,rad,index)        # Front surface at 0.0 
        self.add(sl)
        sr = sur.SphericalSurface(t,cr,rad,wl.AirIndex())  # back surface at t
        self.add(sr)

        self.minThickness = 2.0                            # Sanity thickness

        
    def __str__(self):
        """
        Overwrite of str for a singlet.
        """
        return "pt: {0:s} bfl: {1:6.4f} radius: {2:6.4f}, bend: {3:6.4f} thickness: {4:6.4f}".\
            format(str(self.point),self.backFocalLength(),self.getRadius(),self.getBend(),self.getThickness())

    def getBend(self):
        """
        Get the bend parameter of the lens given by the back and front curvatures.

        :return: the bend of the lens as a float.

        """
        cl = self[0].curvature
        cr = self[1].curvature
        return (cl + cr)/(cl - cr)

    def setCurvatures(self,front= None,back = None):
        """
        Re-set the front and back curvatues of the lens buy keeps other parameters the same.
        
        :param front: the front curvature
        :type front: float
        :param back: the back curvature
        :type back: float

        """
        if front != None:
            self[0].curvature = front
        if back != None:
            self[1].curvature = back
        self.paraxial = None
        return self

    def invert(self):
        """
        Method to invert the lens by swapping left and right curvatures. Radius and location is not changed.

        :return: self but inverted. 
        """
        self.setCurvatures(-self[1].curvature,-self[0].curvature)
        return self

    
    def setFocalLength(self,fl, fixedradius = False, wave = wl.Design):
        """
        Set the focal length by scaling with the option to retain the current radius. This overwrites the method in Lens.

        :param fl: the new focal length
        :type fl: float
        :param fixedradius: flag to fix the radius (Default = False)
        :type fixedradius: bool
        :param wave: the wavelength (Defaul = wl.Design)
        :type wave: float

        """
        r = self.getRadius()
        Lens.setFocalLength(self,fl,wave)
        if fixedradius:
            self.setRadius(r)
        return self

    def getRadius(self):
        """
        Get the radius of the lens, being the maxRadius of the first surface.
        """
        return self[0].maxRadius


    def getFNo(self,wave = wl.Design):
        """
        Get the FNo, so focallength / diameter
        """
        return self.focalLength(wave)/(2.0*self.getRadius())

    def setRadius(self,r):
        """
        Sets max radus of the two surfaces
        """
        self[0].maxRadius = abs(r)
        self[1].maxRadius = abs(r)
        self.paraxial = None
        return self
        

    def setBend(self,bend = 0.0, fixedfocal = False ):
        """
        Set the bend of the lens of the lens by varying the curvatures
        param bend the bend parameter, this can be numerical of sting of "biconvex", "planoconvex" or "convexplano"
        param focal, if True (default) then lens will be scaled to retain current focal length at default wavelnegth
        Standard values are for positive lens are:

        :param bend: the bend, either numerical of str (Default = 0.0)
        :type bend: float or str
        :param fixedfocal: if True will scale to retain the local length and radius (Default = False)
        :type fixedfocal: bool

        """

        if isinstance(bend,str):
            bend = bend.lower().strip()
            if bend.startswith("biconvex"):
                beta = 0.0
            elif bend.startswith("planoconvex"):
                beta = -1.0
            elif bend.startswith("convexplano"):
                beta = 1.0
            else:
                print("Simple single with unknown shape parameter, setting to biconvex")
                beta = 0.0
        else:
            beta = float(bend)
        
        f = self.backFocalLength()
        c = self[0].curvature - self[1].curvature
        self.setCurvatures(0.5*c*(beta + 1.0),0.5*c*(beta - 1.0))
        self.paraxial = None
        if fixedfocal:                    # Scale focal length if required.
            self.setFocalLength(f,True)
        return self

    def getThickness(self):
        """
        Get the thickness of the lens as the centre.

        :return: the centre thickness as a float

        """
        return self[1].point.z - self[0].point.z

    def setCentreThickness(self, t = 0.0):
        """
        Method to set the centre thickness by moving the back surface. This will typically change the focal length.
        It will also check that the centre or edge thickness is not below the citera set in self.minThickness = 2mm
        
        :param t: the centre thickness, if set to 0.0 will make lens as thin as possible so that edge or central.
        :type t: float

        """
        t = max(t,self.minThickness)
        et = self.getEdgeThickness()
        ct = self.getThickness()
        move = t - ct
        if et + move < self.minThickness:    # Too thin at edge
            move = self.minThickness - et

        self[1].point.z += move
        self.paraxial = None
        return self

    def getEdgeThickness(self):
        """
        Get the edge thickness of the lens
        
        :return: the edge thickness in mm

        """
        front = self[0].point.z + self[0].edgePlane()
        back = self[1].point.z + self[1].edgePlane()
        return back - front

    def setEdgeThickness(self,t = 0.0):
        """
        Method to set the edge thickness by moving the second surface, i

        :param t: target edge thickness, if 0.0 will be self.minThickness.

        """
        t = max(t,self.minThickness)
        et = self.getEdgeThickness()
        ct = self.getThickness()
        move = t - et           
        if ct + move < self.minThickness: # Too thin in centre
            move = self.minThickness - ct

        self[1].point.z += move
        self.paraxial = None
        return self

    def setThickness(self,t = 0.0):
        """
        Set lens so that edge / centre is set to sepcified value. Called with t = 0 (or default)
        will set the current lens to the thinnest possible.

        :param t: target thickness
        :type t: float in mm
        :return: self

        """
        self.setEdgeThickness(t)
        self.setCentreThickness(t)
        return self



    def setParameters(self,focal,radius,thick = 0.0, wave = wl.Design):
        """
        Method to set the normal optical paramteters of the current lens by scaling, the bend is retained.
        
        :param focal: The target focal length in mm
        :type focal: float
        :param radius: the radius
        :type radius: float
        :param thick: the thickness, either centre of edge depending of whick is less, default = 0.0 is thinnest possible.
        :type thick: float
        :param wave: the  wavelength (default = wl.Default)
        :type wave: float

        Note this is itterative since that all repend on each other in a non-lienar way !
        """
        while abs(self.backFocalLength(wave) - focal) > abs(focal*1.0e-6) : 
            self.setFocalLength(focal,wave)   # scale to get focal length right
            self.setRadius(radius)            # Set correct radius
            self.setThickness(thick)          # Optimise thickness
        return self

    def setFromString(self,string):
        """
        Method to set parameters of a lens from a string with keywords. Tokens are processed in order.
        """

        string = string.strip()
        token = string.split()
        ntokens = len(token)
        next = 0
    
        #      Process tokens on order

        while next < ntokens:
            if token[next].startswith("point"):            # Deal with Point
                next += 1
                v = eval(token[next])
                next +=1
                self.setPoint(v)
            elif token[next].startswith("focal"):
                next += 1
                f = float(token[next])
                next += 1
                self.setFocalLength(f)
            elif token[next].startswith("radius"):
                next += 1
                r = float(token[next])
                next += 1
                self.setRadius(r)
            elif token[next].startswith("fno"):
                next += 1
                fno = float(token[next])
                next += 1
                r = 0.5*abs(self.backFocalLength())/fno
                self.setRadius(r)
            elif token[next].startswith("bend"):
                next += 1
                self.setBend(token[next])
                next += 1
            elif token[next].startswith("thick"):
                next += 1
                t = float(token[next])
                next += 1
                self.setThickness(b)   
            elif token[next].startswith("index"):
                next += 1
                index = wl.MaterialIndex(token[next])
                next += 1
                self[0].refractiveindex = index

                
            else:
                print("Unknown token")

        return self


    def draw(self,planes = True, legend = False):
        """
       Methoid to draw a single whith or without paraxial planes.

        :param planes: draw the paraxial planes (Default = True)
        :type planes: bool
        :param legend: draw legend panel (Default = False)
        :type legend: bool
        """

        fp = self[0].getPoint()
        bp = self[1].getPoint()
        self[0].draw()            # Front surface
        self[1].draw()            # Back surface
        zf = fp.z + self[0].edgePlane()
        zb = bp.z + self[1].edgePlane()
        y = fp.y + self.getRadius()
        plot([zf,zb],[y,y],"k",lw=2.0)  # Horizontal lines at edge
        plot([zf,zb],[-y,-y],"k",lw=2.0)

        
        if planes:                     # Add the planes if wanted
            if self.paraxial == None:
                pg = self.paraxialGroup()
            self.paraxial.draw(legend)

        
        

class SimpleSinglet(Singlet):        
    """
    Class to make a simple singlet with specified focal length, radius, type / bend 

    :param pt_or_z: location of lens, either ray.Position or location on z-axis (default = 0.0)
    :type pt_or_z: Vector3d or float
    :param f: focal length (default = 100.0)
    :type f: float
    :param r: radius  (default = 10.0)
    :type r: float
    :param bend: This can also be specified as string, being "biconvex", "planoconvex", or "convexplano" of numerically.
    :type bend: str or float
    :param index: refarctiveindex, can be glass key, (default = "BK7")
    :type index: RefractiveIndex or str

    """
    def __init__(self,pt_or_z = 0.0 ,f = 100.0 ,r = 10.0 ,bend = 0.0, index = "BK7"):
        """
        Simple lens:
       
        """
        if isinstance(bend,str):
            bend = bend.lower().strip()
            if bend.startswith("biconvex"):
                bend = 0.0
            elif bend.startswith("planoconvex"):
                bend = -1.0
            elif bend.startswith("convexplano"):
                bend = 1.0
            else:
                print("Simple single with unknown shape parameter, setting to biconvex")
                bend = 0.0
            
        n = 1.51                  #   Inital guess at parameters  and index and c
        c = 1.0/(2.0*f*(n - 1))   
        t = 0.25*r
        Singlet.__init__(self,pt_or_z,c,t,-c,r,index) # Make a biconvex lens
        self.setBend(bend)                            # Bend to right shape 
        self.setParameters(f,r,0.0)                   # Set the actual parameters.
        self.title  = "Simple Singlet" 

CurrentLens = SimpleSinglet()

class Doublet(Lens):
    """
    Class to implement a achromatic double lens with simpler interface than OpticalGroup and additional 
    method to alter the lens. The default is BK7 / F4 Franhoffer doublet with biconvex crown lens, flat back surface, 10mm diameter and focal length of 106mm

    :param pt_or_z: the group point (default = 0.0)
    :type pt_or_z: Vector3d or float
    :param c1:  curvature of front (Default = 0.0225)
    :type c1: float
    :param tf: thickness at centre first lens (default = 5.0)
    :type tf: float
    :param cm: curvature of common surface between lenses (default = -0.0225)
    :type cm: float
    :param ts: thickess at centre of second lens (default = 2.0)
    :type ts: float
    :param cr: float curcature of back element (default = 0.0)
    :type cr: float
    :param radius: max radius (defaults = 10.0)
    :type radius: float
    :param crownindex:  RefrativeIndex or first element material key, (default = "BK7")
    :type crownindex: RefractiveIndex or str
    :param flintindex: RefartiveIndex of second material (default = "F4")
    :type flintindex: RefratciveIndex or str

    """
    def __init__(self, pt_or_z = 0.0, cl = 0.0225 , tf= 5.0 ,cm = -0.0225, ts = 2.0, cr = 0.0 ,rad = 10.0,crownindex="BK7", \
                 flintindex = "F4"):
        """
        Basic constructor for Doublet
        
        """
        Lens.__init__(self,pt_or_z)
        
        if isinstance(crownindex,str):                        # Lookup  crownindex if given string key
            crownindex = wl.MaterialIndex(crownindex)
        if isinstance(flintindex,str):                        # Lookup flintindex if given string key
            flintindex = wl.MaterialIndex(flintindex)

        #           Add the three surfaces
        sl = sur.SphericalSurface(0.0,cl,rad,crownindex)        # Front surface at 0.0 
        self.add(sl)
        sm = sur.SphericalSurface(tf,cm,rad,flintindex)         # Common middle surface
        self.add(sm)
        sr = sur.SphericalSurface(tf + ts,cr,rad,wl.AirIndex())  # Back surface
        self.add(sr)

        self.minThickness = 2.0                            # Sanity thickness

        #
    def invert(self):
        """
        Method to invert the Doublet 
        
        """
        self[0].curvature,self[2].curvature = -self[2].curvature,-self[0].curvature
        self[1].curvature = -self[1].curvature
        self[0].refractiveindex,self[1].refractiveindex = self[1].refractiveindex,self[0].refractiveindex
        t = self[2].point - self[1].point
        self[1].setPoint(self[0].point + t)
        self.paraxial = None
        return self


    def setFocalLength(self,f,fixedradius = False, wave = wl.Design):
        """
        Set the focal length by scaling with the option to retain the radius.

        :param f: target focal length
        :type f: float
        :param fixedradius: option to retain the radius (Default = False)
        :type fixedradius: bool
        :param wave: the wavelnngth, (Default = optics.wavelength.Default)
        :type wave: float

        """
        r = self.getRadius()
        Lens.setFocalLength(self,f,wave)
        if fixedradius:
            self.setRadius(r)
        return self

    def getRadius(self):
        """
        Get the radius of the lens.

        :return: float, the radius of the lens.

        """
        return self[0].maxRadius


    def getFNo(self,wave = wl.Design):
        """
        Get the FNo, so focallength / diameter

        :param wave: the wavelength (Default = optics.wavelength.Default)
        :type wave: float
        :return: The FNo

        """
        return self.backFocalLength(wave)/(2.0*self.getRadius())

    def setRadius(self,r):
        """
        Sets max radus of the three surfaces

        :param r: the radius
        :type r: float

        """
        self[0].maxRadius = abs(r)
        self[1].maxRadius = abs(r)
        self[2].maxRadius = abs(r)
        self.paraxial = None
        return self

    def setCurvatures(self,front=None,centre=None,back=None):
        """
        Set or re-set the front, centre,  and back curvatues of the lens.

        :param: front: the front curvature. (Default =  None) (not changed)
        :param centre: the centre curvature (Default = None)
        :param back: the back curvature (Default = None)
        
        """

        if front != None:
            self[0].curvature = front
        if centre != None:
            self[1].curvature = centre
        if back != None:
            self[2].curvature = back
        self.paraxial = None
        return self


    def getSingletPair(self):
        """
        Get the two parts of the doublet as list of two Singlets

        :return: list of two Singlets, being the two component lenses of the doublet.

        """
        ft = self[1].point.z - self[0].point.z          # First thickness
        bt = self[2].point.z - self[1].point.z          # back thickness
        front = Singlet(self[0].point,self[0].curvature,ft,self[1].curvature,self.getRadius(),self[0].refractiveindex)
        bp = self[0].point + Vector3d(0,0,ft)
        back = Singlet(bp,self[1].curvature,bt,self[2].curvature,self.getRadius(),self[1].refractiveindex)
        return front,back


    def draw(self,planes = True, legend = False):
        """
        Draw a double at two singlets

        :param planes: Boolean to add the paraxial planes (Default = True)
        :type planes: bool
        :param legend: Add plane legend to plot (Default = False)
        :type legend: bool

        """

        front,back = self.getSingletPair()
        front.draw(False)
        back.draw(False)
        if planes:                     # Add the planes if wanted
            if self.paraxial == None:
                pg = self.paraxialGroup()
            self.paraxial.draw(legend)
        
        
#
class Eye(Lens):
    """
    Class to model the eye using values from Hyperphysics and guesses at Abbe Numbers
    """
    def __init__(self,pt, pixels = 0):
        """
        Param pt, the location of the first surface
        """
        Lens.__init__(self,pt)

        self.lensfrontcurvature = 0.11534          # Paramaters of crystaline lens
        self.lensfrontposition = 3.24
        self.lensbackcurvature = -0.15798
        self.lensbackposition = 8.22
        #
        #                      Add Cornea
        self.add(sur.SphericalSurface(0.0, 0.13774, 4.0, wl.CauchyIndex(1.376,50.0)))
        self.add(sur.SphericalSurface(0.45, 0.17606, 4.0, wl.CauchyIndex(1.336,53.0)))
        #                      Add pupil
        self.add(sur.IrisAperture(3.0,2.5))
        #                      Add cytstaline lens with default paramters
        self.add(sur.SphericalSurface(self.lensfrontposition, self.lensfrontcurvature, 3.0, \
                                      wl.GradedIndex([0.0,0.0], wl.CauchyIndex(1.406,50.0),[1.0,-1.5805e-3])))
        self.add(sur.SphericalSurface(self.lensbackposition, self.lensbackcurvature, 3.0, \
                                      wl.CauchyIndex(1.337,53.0)))
        #                      Find back focal plane at peak sensitivity 
        bfp = self.paraxialGroup(wl.PhotopicPeak).backFocalPlane() - self.point.z
        #                      Add in curved image plane as retina (curve of 1/12 mm)
        print("Back focal plane at : " + str(bfp))
        if pixels == 0:
            self.retina = sur.SphericalImagePlane(bfp,-8.333e-2,4.0)
        else:
            self.retina = ana.CurvedOpticalImage(bfp,-8.333e-1, 4.0, 4.9, pixels, pixels)
        self.add(self.retina)


    def accommodation(self,a):
        """
        Method to simulate accomodation by altering the thickness and curvatutes of the crystaline lens
        """
        t = self.lensbackposition - self.lensfrontposition             # Thickness of lens
        cratio = self.lensfrontcurvature/abs(self.lensbackcurvature)   # Ratio of curvatutes
        ft = cratio*t/(1.0 + cratio)                                   # front/back thcikness in ratio
        bt = t/(cratio * (1.0 + 1.0/cratio))
        centre = self.lensfrontposition + ft                           # centre of lens
        #
        #               form the new front and back positions by scaling the front/back surfaces
        #               while keeping the centre fixed
        frontposition = centre - a*ft                              
        backposition = centre + a*bt
        #
        #               Increase the curvatures by a**3 (impirical)
        frontcurve = self.lensfrontcurvature*a**3
        backcurve = self.lensbackcurvature*a**3
        #               Set the values of surface 3 and 4 (front/back of crystaline lens)
        self[3].point.z = frontposition
        self[3].curvature = frontcurve
        self[4].point.z = backposition
        self[4].curvature = backcurve

        self.paraxial = None
        return self


    def setNearPoint(self,distance):
        """
        Method to use the accommoation to set the NearPoint, this needs be a simple itterative
        scheme.
        param distance, near point distance from front Nodal point.
        return float the accommodation parameter.
        """
        #                     Fine obect point wrt to front nodal
        opt = Vector3d(self.frontNodalPoint(wl.PhotopicPeak) - Vector3d(0,0,distance))
    
        da = 0.1
        a = 1.0           # Initial a value
        while True:

            self.accommodation(a)
            ipt = self.imagePoint(opt,wl.PhotopicPeak)
            delta = ipt.z - self.retina.point.z
            if abs(delta) < 1e-4:
                break
                
            if delta > 0 and da < 0:    # Reverse on overshoot and / step by 3
                da = -da/3.0
            if delta < 0 and da > 0:
                da = -da/3.0
            
            a += da                     # Update 

        return a

#
class DataBaseLens(Lens):
    """
    Class to read lens from input file 
    """

    def __init__(self,fn = None):
        """
        Read in the lens from specified file.
        param fn the name of the lens file
        If the filename does not end in lens then ".lens" is appended

        If this is no filename given then the user will me prompted via tio openFile.
        """
        
        Lens.__init__(self,0.0)         # Create a blank lens

        #
        #         Open file, if None then prompt via tio interface
        #
        if fn == None:
            lensfile = tio.openFile("Lens file","r","lens")
        else:
            fn = tio.getExpandedFilename(fn)   # Sort out logicals
            if not fn.endswith("lens"):        # Append ".lens" if not given
                fn += ".lens"
            lensfile= open(fn,"r")             # open file

        #          read file and process one line at a time
        #
        for line in lensfile.readlines():
            line = line.strip()
            if not line.startswith("#") and len(line) > 0:   # Kill comments and blanks
                token = line.split()
                
                if token[0].startswith("title"):              # Deal with title
                    self.title = str(token[1])

                elif token[0].startswith("point"):            # Deal with Po
                    v = eval(token[1])
                    self.setPoint(v)                          # Set point

                elif token[0].startswith("iris"):             # Iris aperture
                    p = float(token[1])                       # z-position
                    r = float(token[3])                       # radius
                    s = sur.IrisAperture(p,r,1.0)
                    self.add(s)
                    
                elif token[0].startswith("circular"):         # Circular aperture
                    p = float(token[1])                       # z-poistion
                    r = float(token[3])                       # radius
                    s = sur.CircularAperture(p,r)
                    self.add(s)

                elif token[0].startswith("annular"):          # Annular aperture
                    p = float(token[1])                       # z-poistion
                    inner = float(token[3])                   # inner radius
                    outer = float(token[4])                   # outer radius
                    s = sur.AnnularAperture(p,inner,outer) 
                    self.add(s)

                elif token[0].startswith("spherical"):        # Spherical surface
                    p = float(token[1])                       # z-position
                    c = float(token[3])                       # curvature
                    r = float(token[5])                       # max radius
                    if token[6].startswith("index"):          # refrective index
                        index = wl.MaterialIndex(token[7])
                        s = sur.SphericalSurface(p,c,r,index) 
                    elif token[6].startswith("mirror"):       # else its a mirror
                        s = sur.SphericalSurface(p,c,r)
                    else:
                        raise IOError("unknow type")
                    self.add(s)

                elif token[0].startswith("parabolic"): # Parabolic
                    p = float(token[1])
                    c = float(token[3])
                    r = float(token[5])
                    if token[6].startswith("index"):
                        index = wl.MaterialIndex(token[7])
                        s = sur.ParabolicSurface(p,c,r,index)
                    elif token[6].startswith("mirror"):
                        s = sur.ParabolicSurface(p,c,r)
                    else:
                        raise IOError("unknow type")
                    self.add(s)

                elif token[0].startswith("quadric"): # Quadric
                    p = float(token[1])
                    c = float(token[3])
                    e = float(token[5])
                    r = float(token[7])
                    if token[8].startswith("index"):
                        index = mat.MaterialData().getIndex(token[9])
                        s = sur.QuadricSurface(p,c,r,e,index)
                    elif token[8].startswith("mirror"):
                        s = sur.QuadricSurface(p,c,e,r)
                    else:
                        raise IOError("inknow  index ")
                    self.add(s)

                else:
                    print("Unknown token : " + str(token[0]))
                    
            lensfile.close()             # close file


class OpticalSystem(Lens):
    """
    Class to hold an optical system consisting of a list of OpticalGroups.
    """

    def __init__(self,title = "System", *args):
        """
        Bacic constructor
        """
        list.__init__(self)
        self.title = title
        self.iris = None              # System Iris
        for g in args:
            self.add(g)

        

            
    def add(self,g):
        """
        Method to add a grpup, simple at the moment.

        
        """
        if hasattr(g,"iris"):
            if g.iris != None:
                self.iris = g.iris
        self.append(g)


    def draw(self,planes = True):
        """
        Method to draw the groups in order the plt.plot

        """
        for g in self:
            g.draw(planes)


    def entranceAperture(self):
        """
        Get the entrance aperture of the system which is the  entrance aperture to the first component.

        :return: a CircularAperture being the entrance aperture

        """
        return self[0].entranceAperture()

    def exitApeture(self):
        """
        Get the exit aperture of the system, being the exit apertrure of the last component.

        :return: a CircularAperture being the exit aperture

        """
        return self[-1].exitAperture()
    

    def setIris(self,ratio,group = None):
        """
        Method to set the system iris of the system, or the iris is a particular group

        :param ratio: the ratio between 0 and 1 
        :type ratio: float
        :param group: group with the iris (Default = None)
        :type group: int or None
        
        """
        if isinstance(group,int):       # Specific group
            self[group].setIris(ratio)
        elif self.iris != None:         # Sytsem iris set
            self.iris.setRatio(ratio)
        return self

    def setPoint(self,pt,group):
        """
        Method to set the reference point of a specified group
        """
        self[group].setPoint(pt)
        return self

    def movePoint(self,delta,group):
        """
        Move the reference point by about delta

        :param delta: distance moved
        :type delta: Vector3d or float
        """
        self[group].movePoint(delta)
        return self


    def paraxialGroup(self,wave = wl.Design):
        """
        Get the Paraxialmatrix for the whole system

        :param wave: the design wavelength (Default = wl.Default)
        :type wave: float
        :return: the ParaxialGroup for the whole system.
        """
        m = self[0].paraxialGroup(wave).copy()
        for g in self[1:]:
            m *= g.paraxialGroup(wave)

        return m
