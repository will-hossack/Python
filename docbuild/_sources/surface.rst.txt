===============
Surface Classes
===============

Set of classes to respesent various types of optical planes and surfaces.
These class types are the main obects used in tracing of rays.


Surface Constants
=================

The following constants define surface types to

* Clear surfaces that do not change rays paths, these include image planes,
  apertures, input / output plane etc

* Refrating surfaces, flat or curved glass surfaces that refract rays.

* Reflecting surfaces, flat or curved mirrored surfaces.

* Constant Blocked to represent a blocked rays.

.. automodule:: optics.surface
   :members: Clear,Refracting,Reflecting,Blocked




Surface Class
=============

Abstract base class to give represent a general surface. This class defines

* Surface reference point being a Vector3d that defines the location of the surface.

* Optical Group that the surface belongs to, may be None

* Type of surface.

* Refractive index on the image side, which may be None.

.. autoclass:: optics.surface.Surface
   :members:


FlatSurface Class
=================

Class to represent a general Flat optical surface with specified surface
normal.

.. autoclass:: optics.surface.FlatSurface
   :members:


OpticalPlane Class
==================

Class to represent a flat optical plane pertendicular to the optical axis
so fixed surface normal.

.. autoclass:: optics.surface.OpticalPlane
   :members:

CircularAperture Class
======================

Class to represent a circular aperture of specifed fixed radius.

.. autoclass:: optics.surface.CircularAperture
   :members:


AnnularAperture Class
=====================

Class to represent a circular annular aperture of specifed fixed inner and
 outer radius.

.. autoclass:: optics.surface.AnnularAperture
   :members:

  
IrisAperture Class
======================

Class to represent a circular iris aperture of variable radius.

.. autoclass:: optics.surface.IrisAperture
   :members:


KnifeAperture Class
=======================

Class to implement a knife edge aperture for knife-edge test.

.. autoclass:: optics.surface.KnifeAperture
   :members:

   
ImagePlane Class
================

Class to represent an image plane, this is the same as OpticalPlane
but wuth additial size information. This affets the .draw() methodd.

.. autoclass:: optics.surface.ImagePlane
   :members:

QuadricSurface Class
====================

Class to represent a quadratic surface.

.. autoclass:: optics.surface.QuadricSurface
   :members:

SphericalSurface Class
======================

Class to represent a Spherical Surface

.. autoclass:: optics.surface.SphericalSurface
   :members:


ParabolicSurface Class
======================

Class to represent a Spherical Surface

.. autoclass:: optics.surface.ParabolicSurface
   :members:

SurfaceInteraction Class
========================

Class to return the interaction of a ray with a surface. Each surface an this object
via the method .getSurfaceInteraction(ray). It is then used to act on the ray to update
it.

%.. autoclass:: optics.surface.SurfaceInteraction
%   :members:

 There are no methods associated with this class, it is used to transfer information to
 update the ray.
