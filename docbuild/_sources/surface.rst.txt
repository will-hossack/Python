===============
Surface Classes
===============


Surface Constants
=================

The following constants define surafce types

.. automodule:: optics.surface
   :members: Clear,Refracting,Reflecting


Classes to represent optical surfaces

Surface Class
=============

Abstract base class to give represent a general surface. 

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

