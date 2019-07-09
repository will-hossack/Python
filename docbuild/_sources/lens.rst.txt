============
Lens Classes
============

Set of claases to handel full ray tracing of lenses.


OpticalGroup Class
==================

Basic class for lenes to hold list of OpticalSurfaces.

.. autoclass:: optics.lens.OpticalGroup
   :members:

Lens Class
==========

Class extending for lens with extra methods to make it s easier to call.

.. autoclass:: optics.lens.Lens
   :members:

Singlet Class
=============

Class for a general singlet lens with extra methods to set the parameters.

.. autoclass:: optics.lens.Singlet
   :members:

SimpleSinglet Class
===================

Simpler interface to Singlet to form a thin singlet lens specified by focal length, bend and radius.

.. autoclass:: optics.lens.SimpleSinglet
   :members:


Doublet Class
=============

Class to handle a doublet

.. autoclass:: optics.lens.Doublet
   :members:
      
DataBaseLens Class
==================

Class to read a lens from an input file.

.. autoclass:: optics.lens.DataBaseLens
   :members:
 
