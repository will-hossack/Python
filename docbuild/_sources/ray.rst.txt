===========
Ray Tracing
===========

Ray Class
=========

The base class in the Ray, which is abstract with the useful extensing
classes being ParaxialRay and Intensity ray.

.. autoclass:: optics.ray.Ray
   :members:

RayMonitor Class
================

A RayMonitor can be attached to each Ray that is automatically updated every time the ray propated. The abstarct RayMonitor
class is specifed as:

.. autoclass:: optics.ray.RayMonitor
   :members:

but is actually used via its extending clases.

.. autoclass:: optics.ray.PrintPath
   :members:

which prints the values of the ray as it is updated, and

.. autoclass:: optics.ray.RayPath
   :members:

which records the ray path in internal lists that can then be plotted via the .draw() method.
      

ParaxialRay Class
=================

This class definbes Paraxial Rays which work with the optics.matrix classes.

.. autoclass:: optics.ray.ParaxialRay
   :members:

Paraxial Ray Tracing
====================

.. code-block:: python

   import optics.matrix as m
   import optics.ray as r



Intensity Ray
=============

The IntensityRay is the main ray raytype for full ray traceing

.. autoclass:: optics.ray.IntensityRay
   :members:

SourcePoint
===========

Class to implement a source of ray, being a vector position with assoiated intensity or spectrum

.. autoclass:: optics.ray.SourcePoint
   :members:

RayPencil
=========

The most useful and powerful class for handling rays is the RayPencil, beinng a list or rays that can be propagated through surfaces with the
"\*=" operator. There are also a set of powerful method that allow the creation of various types or RayPencils.

.. autoclass:: optics.ray.RayPencil
   :members:


  
