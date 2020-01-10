============
Lens Classes
============

Set of claases to handel full ray tracing of lenses.


CurrentLens functions
=====================

There are two functions that control the default lens used in the package.
These are mainly to simply the GUI interface.

.. autofunction:: optics.lens.setCurrentLens

and the getter function

.. autofunction:: optics.lens.getCurrentLens

If there is no CurrectLens set, it defaults the default SimpleSinglet.

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
 
