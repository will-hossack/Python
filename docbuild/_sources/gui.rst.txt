===========
GUI Classes
===========

There is a GUI based on PyQt5 widgets, most of the graphics use Matplotlib
plotting in a Qt5 window.

PltMainWindow Class
===================

Class to form the main window with basic menu to read in lenses from database,
show lens, set wavelength, angle etc. This is the main class called by the more usable extending classes.

.. autoclass:: optics.gui.PltMainWindow
   :members:

LensViewer Class
================

 Class to display a lens with a collimated beam being imaged in the back focal
 plane. Menus allow change of lens, wavelength, angle and iris ratio.

 .. autoclass:: optics.gui.LensViewer
    :members:

AberrationViewer Class
======================

Class to display the aberration plot of a lens with collimated illumination
with menues to change lens, wavelngth, angle and isis ratio.

.. autoclass:: optics.gui.AbberationViewer
   :members:



SpotViewer Class
================

Class to implement an interactive spot diagram with geometric psf for a lens
with a collimated input beam. There are menus to change the lens, wavelength,
iris ratio and interactive movement of the spot plane.

.. autoclass:: optics.gui.SpotViewer
   :members:



 
