"""
Set of classes to implement various types lenses, being lists of surfaces, with additional
method to extract the geomertic parameters in a simple way.
"""
import optics.surface as sur
from optics.matrix import ParaxialPlane,ParaxialMatrix,DielectricMatrix,ParaxialGroup
from vector import Vector3d,Unit3d,Angle
from optics.ray import IntensityRay
from optics.wavelength import getDesignWavelength,AirIndex,MaterialIndex,CauchyIndex,PhotopicPeak
import optics.analysis as ana
import tio
from matplotlib.pyplot import plot
import math



def setCurrentLens(lens):
    """
    Function to set the Current Lens or OpticalSystem, it is held in global CurrentLens
    This is mainly used for GUI interface. 
    
    Currentlens defaults to a default SimpleSinglet() 

    :param lens: the lens or a filename, if it is a str, it will try and open the DataBaseLens of this type.
    :type lens: OpticalGroup or extenting class.
    """
    global CurrentLens
    if isinstance(lens,str):
        CurrentLens = DataBaseLens(lens)
    else:
        CurrentLens = lens

def getCurrentLens():
    """ 
    Function to get the current default lens, if this is not
    overriddent it is initally set to the default SimpleSinglet()
    This is typically used in the GUI interface

    :return: Current default lens
    """
    return CurrentLens




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
        self.wavelength = getDesignWavelength()         # Default wavelength 
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
            self.point = Vector3d(pt) # Make local copy of vector
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
        if isinstance(delta,float) or isinstance(delta,int):
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
    def planePair(self,mag, xsize = 100.0, ysize = None , wave = None):
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


    def setWithPlanes(self,obj,ima,wave = None):
        """
        Set the position of the OpticalGroup  to match the supplied object and image
        planes. Also if  object plane has a sopecified height the image plane height
        will also be set to match the magmification.

        :param obj: the object plane
        :type obj: ImagePlane
        :param ima: the image plane
        :type ima: ImagePlane
        :param wave: the deign wavelengh, (Default = wl.Design)
        :type wave: float
        :return: magificantion as a float
        """
        xsize,ysize = obj.getSize()
        pobj = ParaxialPlane(obj.getPoint().z,ysize)      # Make paraxial planes
        pima = ParaxialPlane(ima.getPoint().z)
        pm = self.paraxialGroup(wave)                            # Paraxial matrix
        mag = pm.setWithPlanes(pobj,pima)                        # Set the planes
        pt = self.getPoint()                                     # get and update point to  move lens
        pt.z = pm.inputPlane()
        self.setPoint(pt)     
        xsize *= abs(mag)
        ysize *= abs(mag)
        ima.setSize(xsize,ysize)
        return mag

        
    def imagePoint(self,op,wave = None):
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



    def paraxialMatrix(self,wave = None , first = 0, last = -1):
        """
        Get a the paraxial matrix for a subset of surfaces.

        :param wave: the wavelngth (Default = optics.wavelength.Design)
        :type wave: float
        :param first: first surface (Default = 0)
        :type first: int
        :param last: last surface (Default = -1)
        :type last: int

        """
        wave = getDesignWavelength(wave)
        mat = ParaxialMatrix()           # Start with unit matrix
        if first >= len(self):
            raise IndexError("lens.OpticalGroup.paraxialMatrix, index out of range")
        
        if first == 0 :
            z = 0
        else:
            z = self[first].point.z               # Start point
        nl = AirIndex().getValue(wave) 
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
                m = DielectricMatrix(nl,nr,s.curvature)
                mat *= m                                       # add surface
                z = zp
                nl = nr
        #     
        return mat
        
            
        

    def paraxialGroup(self,wave = None):
        """
        Get the paraxial group pf this group at spectified wavelength. This will be
        remade if either first call, surface added, wavelength change or scale change.

        :param wave: the wavelength, if None it will use the defaults wavelngth for the Group
        :type wave: float
        :return: the :class:`optics.matrix.ParaxialGroup`

        """
        #              Make paraxial matrix if needed.
        #
        wave = getDesignWavelength(wave)
        if self.paraxial == None or wave != self.wavelength: # make a new matrix
            self.wavelength = wave
            mat = self.paraxialMatrix(wave)
            en = self.entranceAperture()               # Get input/output heights
            ex = self.exitAperture()
            #  Make Paraxial Group with matching ref point, blank matrix cortect input/output heights
            self.paraxial = ParaxialGroup(self.point.z,mat,en.maxRadius,ex.maxRadius)
            
        
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
    

    def cardinalPoints(self,wave = None):
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


    def entrancePupil(self,wave = None):
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

        

    def exitPupil(self,wave = None):
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


    def backFocalLength(self, wave = None):
        """
        Method to get the back focal length calculated by paraxial matrix methods, which for positive lens will be +ve.
        
        :param wave:  specified wavelength (default is wl.Design).
        :type wave: float
        :return: the focal length.

        """
        pm = self.paraxialGroup(wave)
        return pm.backFocalLength()

    def frontFocalLength(self, wave = None):
        """
        Method to get the front focal length calculated by paraxial matrix mecthods, which for positive lens will be +ve.
        
        :param wave:  specified wavelength (default is wl.Default).
        :type wave: float
        :return: the focal length.

        """
        pm = self.paraxialGroup(wave)
        return pm.backFocalLength()

    
        
    def setFocalLength(self,f,wave = None):
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


    def frontFocalPlane(self, wave = None):
        """
        Get the Front Focal Plane as an surface.OpticalPlane is global coordinates.

        :param wave: the wavelength (Default = optics.wavelength.Default)
        :type wave: float
        :return: surface.OpticalPlane in global coordinates.
        """
        pm = self.paraxialGroup(wave)
        zplane = pm.frontFocalPlane()
        return sur.OpticalPlane(Vector3d(self.point.x,self.point.y,zplane))

    def backFocalPlane(self, wave = None):
        """
        Get the back Focal Plane as an surface.OpticalPlane is global coordinates.

        :param wave: the wavelength (Default = wl.Design)
        :type wave: float

        return surface.OpticalPlane in global coordinates.
        """
        pm = self.paraxialGroup(wave)
        zplane = pm.backFocalPlane()
        return sur.OpticalPlane(Vector3d(self.point.x,self.point.y,zplane))


    def frontPrincipalPlane(self,wave = None):
        """
        Get the front principal place as an surface.OpticalPlane in gobal cordilanes

        :param wave: the wavelength (Default = wl.Design)
        :type wave: float
        :return: :py:class:`optics.surface.OpticalPlane` in Global coordinates

        """
        pm = self.paraxialGroup(wave)
        zplane = pm.frontPrincipalPlane()
        return sur.OpticalPlane(Vector3d(self.point.x,self.point.y,zplane))


    def backPrincipalPlane(self,wave = None):
        """
        Get the back principal place in gobal cordilanes as a surface.OpticalPlane in global coordinates

        :param wave: the wavelength (Default = wl.Design)
        :type wave: float
        :return: :py:class:`optics.surface.OpticalPlane` in global coordinates

        """
        pm = self.paraxialGroup(wave)
        zplane = pm.backPrincipalPlane()
        return sur.OpticalPlane(Vector3d(self.point.x,self.point.y,zplane))

    
    def frontNodalPoint(self,wave = None):
        """
        Get the front nodal point as ray.Position in gobal cordilanes

        :param wave: the wavelength (Default = wl.Design)
        :type wave: float
        :return: :py:class:`Vector3d` in Global coordinates

        """
        pm = self.paraxialGroup(wave)
        zplane = pm.frontNodalPoint()
        return Vector3d(self.point.x,self.point.y,zplane)

    def backNodalPoint(self,wave = None):
        """
        Get the front nodal point as ray.Position in gobal cordilanes

        :param wave: the wavelength (Default = optics.wavelength.Default)
        :type wave: float
        :return: Vector3d  in Global coordinates

        """
        pm = self.paraxialGroup(wave)
        zplane = pm.backNodalPoint()
        return Vector3d(self.point.x,self.point.y,zplane)


    def petzvalSum(self,wave = None):
        """
        Calcualte the Perzval field curvature sum assuming air on the left side.

        :return: Perzal sum as a float.
        """
        leftindex = AirIndex()
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
            self.paraxialGroup()    # Force it to be formed
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
            index = MaterialIndex(index)

        sl = sur.SphericalSurface(0.0,cl,rad,index)        # Front surface at 0.0 
        self.add(sl)
        sr = sur.SphericalSurface(t,cr,rad,AirIndex())  # back surface at t
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

    
    def setFocalLength(self,fl, fixedradius = False, wave = None):
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


    def getFNo(self,wave = None):
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



    def setParameters(self,focal,radius,thick = 0.0, wave = None):
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
                self.setThickness(t)   
            elif token[next].startswith("index"):
                next += 1
                index = MaterialIndex(token[next])
                next += 1
                self[0].refractiveindex = index

                
            else:
                print("Unknown token")

        return self


    def draw(self,planes = True, legend = False):
        """
       Method to draw a single whith or without paraxial planes.

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
            self.paraxialGroup()       # Force group to be formed
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

CurrentLens = SimpleSinglet()        # Set package default to simple signlet

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
            crownindex = MaterialIndex(crownindex)
        if isinstance(flintindex,str):                        # Lookup flintindex if given string key
            flintindex = MaterialIndex(flintindex)

        #           Add the three surfaces
        sl = sur.SphericalSurface(0.0,cl,rad,crownindex)        # Front surface at 0.0 
        self.add(sl)
        sm = sur.SphericalSurface(tf,cm,rad,flintindex)         # Common middle surface
        self.add(sm)
        sr = sur.SphericalSurface(tf + ts,cr,rad,AirIndex())  # Back surface
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


    def setFocalLength(self,f,fixedradius = False, wave = None):
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


    def getFNo(self,wave = None):
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
            self.paraxialGroup()
            self.paraxial.draw(legend)
        
        
#
class Eye(Lens):
    """
    Class to model the eye using values from Hyperphysics and guesses at Abbe Numbers
    for the eye lenses
    
    :param pt: group point of eye, front surface of cornea. (Defauylt = 0.0)
    :type pt: Vector3d or float.
    
    :
    """
    def __init__(self,pt = 0.0, pixels = 0):
        """
        Param pt, the location of the first surface
        """
        Lens.__init__(self,pt)

        self.pupilSize = 2.0                    # Pupil radius
        self.lensfrontcurvature = 0.11534          # Paramaters of crystaline lens at distance
        self.lensfrontposition = 3.24
        self.lensbackcurvature = -0.15798
        self.lensbackposition = 8.22
        #
        #                      Add Cornea with surface curvatures and index from Hyperphysics
        #
        self.add(sur.SphericalSurface(0.0, 0.13774, 4.0, CauchyIndex(1.376,50.0)))
        self.add(sur.SphericalSurface(0.45, 0.17606, 4.0, CauchyIndex(1.336,53.0)))
        
        
        #                      Add pupil as a iris aperture of 2.5 mm radius
        #
        self.add(sur.IrisAperture(3.0,self.pupilSize))
        
        #                      Add cytstaline lens with default paramters
        #
        self.add(sur.SphericalSurface(self.lensfrontposition, self.lensfrontcurvature, 3.0,\
                                      CauchyIndex(1.406,50.0)))
        self.add(sur.SphericalSurface(self.lensbackposition, self.lensbackcurvature, 3.0, \
                                      CauchyIndex(1.337,53.0)))
            
        # Record the maximum back focal length
        self.maxFocalLength = self.backFocalLength(PhotopicPeak)
                    
        #                      Find back focal plane at peak sensitivity 
        
    
        bfp = self.paraxialGroup(PhotopicPeak).backFocalPlane() - self.point.z
        
        #                      Add in curved image plane as retina (curve of 1/12 mm)
        if pixels == 0:
            self.retina = sur.SphericalImagePlane(bfp,-8.333e-2,4.0)
        else:
            self.retina = ana.CurvedOpticalImage(bfp,-8.333e-1, 4.0, 4.9, pixels, pixels)
        self.add(self.retina)


    def getRetina(self):
        """
        Get the eye retina
        
        """
        return self.retina


    def accommodation(self,a):
        """
        Method to simulate accomodation by altering the thickness and curvatutes of 
        the crystaline lens
        
        :param a: the accomodation parameter
        :type a: float
        
        
        """
        t = self.lensbackposition - self.lensfrontposition             # Thickness of lens
        cratio = self.lensfrontcurvature/abs(self.lensbackcurvature)   # Ratio of curvatutes
        ft = cratio*t/(1.0 + cratio)                                   # front/back thikkness in ratio
        bt = t/(cratio * (1.0 + 1.0/cratio))
        centre = self.lensfrontposition + ft                           # centre of lens
        #
        #               form the new front and back positions by scaling the front/back surfaces
        #               and moving the back surfacebackward
        frontposition = centre - ft                              
        backposition = centre + a*bt
        #
        #               Increase the curvatures by a**3 (impirical) give almost linear 
        #               realtion between a and focal length
        frontcurve = self.lensfrontcurvature*a**3
        backcurve = self.lensbackcurvature*a**3
        #               Set the values of surface 3 and 4 (front/back of crystaline lens)
        self[3].point.z = frontposition
        self[3].curvature = frontcurve
        self[4].point.z = backposition
        self[4].curvature = backcurve

        self.paraxial = None
        return self
    
    def setFocalLength(self,f):
        """
        Method to set the local length by chaning the accomodation using an impirical fit
        :param f: the focal length 
        :type f: float
        
        """
        if f > self.maxFocalLength:
            print("Cannon exceed maximum focal length")
            return self
        
        #    Use fitted polynomial to concered focal lenth to accomodatio.
        acc = lambda f :  -1.51431847e-04*f**3 + 9.43155993e-03*f**2 - 2.57653636e-01*f + 3.73705241e+00
        
        a = acc(f)
        self.accommodation(a)
        return self
        


    def setNearPoint(self,distance):
        """
        Method to use the accommoation to set the NearPoint, this needs be a simple itterative
        scheme since floal length an principal planes move.
        
        :param distance: near point distance from front Nodal point in mm.
        :type distance: float
        :return: float the accommodation parameter.
        
        """
        #                     Fine obect point wrt to front nodal
        opt = Vector3d(self.frontNodalPoint(PhotopicPeak) - Vector3d(0,0,distance))
    
        da = 0.1
        a = 1.1           # Initial guess
        while True:
            print("a valie is : " + str(a))
            self.accommodation(a)
            ipt = self.imagePoint(opt,PhotopicPeak)
            delta = ipt.z - self.retina.point.z
            if abs(delta) < 1e-4:
                break
                
            if delta > 0 and da < 0:    # Reverse on overshoot and / step by 3
                da = -da/3.0
            if delta < 0 and da > 0:
                da = -da/3.0
            
            a += da                     # Update 

        return a


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
                        index = MaterialIndex(token[7])
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
                        index = MaterialIndex(token[7])
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
                        index = MaterialIndex().getIndex(token[9])
                        s = sur.QuadricSurface(p,c,r,e,index)
                    elif token[8].startswith("mirror"):
                        s = sur.QuadricSurface(p,c,e,r)
                    else:
                        raise IOError("inknow  index ")
                    self.add(s)

                else:
                    print("Unknown token : " + str(token[0]))
                    
            lensfile.close()             # close file



class Prism(OpticalGroup):
    """
    Class to make a simple prism being an equilateral triangle with the 
    top of the prism being vertical and the base being parallel to the optical
    axis. 
    
    :param group_pt: the centre of the prism (Default = (0,0,0))
    :type group_pt: Vector3d or float
    :param angle: prism angle in degrees (Default = 60)
    :type angle: float
    :param height: height of prism from base to peak in mm (Default = 40)
    :type height: float
    :param index: the refactive index (Default = "BK7")
    :type index: RefractiveIndex or str
    
    """
    def __init__(self,group_pt = 0.0, angle = 60.0, height = 40.0, index = "BK7"):
        OpticalGroup.__init__(self,group_pt)

        # Make the prism from two FlatSurfaces at suitable positions and
        # surface normals

        if isinstance(index,str):        # Sort out index
            self.n = MaterialIndex(index)
        else:
            self.n = index
        
        #       Variable needed below
        self.height = float(height)
        self.angle = math.radians(angle)    # Hold angle in radians
        
        #        Make surface normals to front and back faces
        fn = Unit3d(Angle(2*math.pi - self.angle/2))
        bn = Unit3d(Angle(self.angle/2))
        
        #         Poistion of front and back surfaces
        p = height*math.tan(self.angle/2)/2
        
        #         Make the two surfaces and add then to the self (an OpticalGroup)
        #         Not need to specify type are refracting.
        front = sur.FlatSurface(-p,fn,type = 1,index = self.n)
        self.add(front)
        back = sur.FlatSurface(p,bn,type=1,index = AirIndex())
        self.add(back)
        
    def __str__(self):
        """
        The str function

        :return: inrmation on the prism
        """
        return "{0:s} a: {1:5.3f} h: {2:5.3f} n: {3:s}".format(str(self.getPoint()),\
                math.degrees(self.angle),self.height,self[0].refractiveindex.title)
        
    def minDeviation(self,wave = None):
        """
        Get the angle of minium deviation at specified wavelength.
        
        Note this needs to take into account the AirIndex used in the package.

        :param wave: the wavelength (Default = wl.Default)
        :type wave: float
        :return: the minimum deviatation as a float in radians.

        """
        a = AirIndex().getValue(wave)
        nval = self.n.getValue(wave)/a    # Correct for air index
        sa = nval*math.sin(self.angle/2)
        s = math.asin(sa)
        return 2*s - self.angle           # Return radians.
    
    
    def maxResolution(self,wave = None):
        """
        Get the maximum resolution at specifed wavelength at angle of 
        minimum deviation. This is given by d dn/d lambda where d is the maximium 
        pathlength difference given by the size of the prism.
        
        :param wave: the wavelength
        :return: maximum resolution lambda / d lambda as float
        """

        d = 2000.0*self.height*math.tan(self.angle/2) # Max pathlength in microns.
        dn = self.n.getDerivative(wave)          # dn/dy of materail
        return d*dn    # 
    
    
    def resolution(self, radius, wave = None):
        """
        Calcualte the resolution for a specifed input beam radius at 
        angle on minimum deviation. It assumes that the resolution is limited
        by the beam radius and not the size of the prism.
        
        :param radius: radius of input beam
        :type radius: float
        :param wave: wavelength
        :type wave: float
        :return: resolution lamdba / d lambda as a float
        
        """
        dev = self.minDeviation(wave) 
        alpha = dev/2 + self.angle/2 
        
        #      Form path difference between top and bottom of the beam
        d = 4*radius*math.sin(self.angle/2)/math.cos(alpha)
        dmax = 2.0*self.height*math.tan(self.angle/2) # Length of bottom of prism
        if d > dmax:
            d = dmax
            print("Resolution limited by size of prism")
        
        
        dn = self.n.getDerivative(wave)    # dn/d lambda
        return 1000*d*dn                   # scale to microms
        
    
        
    def getInputPoint(self):
        """
        Get the input point in the centre of the first face

        :return: Point on front face in global coordinates
        """
        return self[0].getPoint()
    
    def getOutputPoint(self):
        """
        Get tyhe ooutput point in the centre of the second surface
        
        :return: Point on the back face in global coordinates.
        """
        return self[1].getPoint()
    
    def getOutputAngle(self,inAngle,wavelength = None):
        """
        High level method to get the output angle of a beam given
        input angle and wavelength.
        
        Ths methods traces a ray through the system but hide the mechanics and
        complications of doing this.

        :param inAngle: input angle in radians (from the optical axis)
        :type inAngle: float
        :param wavelength: the wavelength (Default = 0.55)
        :type wavelength: float
        :return: output angle in radians (typically -ve) 
        """
        #        Get prism point and angle of input at Unit3d
        #
        pt = self.getInputPoint()
        u = Unit3d(Angle(inAngle))
        #         Make a ray and trace it
        ray = IntensityRay(pt,u,wavelength)
        ray *= self
        a = ray.getAngle()                 # The the ray Angle
        return a.theta*math.cos(a.psi)     # get in radians allowing for -ve
            
    
    def getWavelength(self, inAngle, outAngle, wavelengths = [0.25,1.0]):
        """
        Itterative methiod to get the wavelength given input and output angles.
        
        Note if impossible angles given or wavelength goes out of range, error message
        and NaN returned.
        
        :param inAngle: input angle in radians
        :type inAngle: float
        :param outAngle: output angle in radians
        :type outAngle: float
        :param wavelengths: range of allowable wavelengths [min,max]
        :type: wavelengths: list of length 2
        :return: the wavelength as a float, "NaN" if it fails to converge.

        """
    
        #        Get prism point and angle of input at Unit3d
        #
        pt = self.getInputPoint()
        u = Unit3d(Angle(inAngle))
    
        #         Guess at initial wavelngth
        wave = (wavelengths[1] - wavelengths[0])/2
        #         Make input ray at guess wavelength
        ray = IntensityRay(pt,u,wave)
    
        #       Parameters for seaerch
        delta = 0.1
        forward = True
        na = float("inf")   # New angle
    
        while abs(na - outAngle) > 1.0e-9/abs(outAngle) :
            nray = ray*self       #      New Ray through prism
            na = nray.getAngle()
            na = na.theta*math.cos(na.psi)    # In radians
            if na < outAngle:                       # Less that target
                wave += delta
                forward = True
            else:   
                if forward:                   # Half step
                    delta *= 0.5
                forward = False
                wave -= delta
            if wave < wavelengths[0] or wave > wavelengths[1]:
                print("Out of wavelength range :")
                return float("nan")
        
            ray.wavelength = wave             # Update the wavelength of ray
        
        return ray.getWavelength()             # End of loop, so success, return value
        
    def draw(self):
        """
        Method to draw the prism

        :return: None

        """
        pt = self.getPoint()     # Centre of prism
        
        #         Form top,left,right corners
        top = [pt.z, pt.y + self.height/2]
        d = self.height*math.tan(self.angle/2)
        left = [pt.z - d , pt.y - self.height/2]
        right = [pt.z + d, pt.y - self.height/2]
        
        #      Plot them out with plt.plot
        plot([top[0],left[0],right[0],top[0]],[top[1],left[1],right[1],top[1]],"k",lw=2.0)




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


    def paraxialGroup(self,wave = None):
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
