=================
Wavelength Module
=================

This module deal with the wavelnegth dependand features in particular refractive and
intensity spectrums.

Wavelength Value
=================

There are a set of globals defined covering the standard wavelength. All are in microns.

.. automodule:: optics.wavelength
   :members: Default,Red,Green,Blue,BlueLimit,RedLimit,BlueColourMatch,GreenColourMatch,RedColourMatch,Mercury_i,Mercury_h,Cadmium_F,Hydrogen_F


Default and Design Wavelength
=============================

The package had two default wavelength, this is the default wave length used in simulations and
the design wavelength which is the default for paraxial calcualtions. These are controlled
by six functions

.. autofunction:: optics.wavelength.getInitialDefaultWavelength

This function is automatically called on startup with the values being held in the global float,

- **optics.wavelength.Default**

This can also be accessed via the fuction

.. autofunction:: optics.wavelength.getDefaultWavelength

This is the advised route to access the default wavelnegth

This default can be changed with a call to

.. autofunction:: optics.wavelength.setDefaultWavelength

which updates the global variable.

There is an exaclty equivalent set of methods for the design wavelength, these being


.. autofunction:: optics.wavelength.getInitialDesignWavelength

which is automatically called on start up with the result being places in global variable.

- **optics.wavelength.Design**

This can also be accessed via the fuction

.. autofunction:: optics.wavelength.getDesignWavelength

which is the advised route to access the design  wavelength

This default can be changed with a call to

.. autofunction:: optics.wavelength.setDesignWavelength

which updates the global variable.

There is an aditional pair of functions mainly used in the GUI interface to handle the current wavelength, these being

.. autofunction:: optics.wavelength.getCurrentWavelength

which return the Current Wavelength, initially set the the Default Wavelength, and the corresponding fundtion to set the current wavelength.

.. autofunction:: optics.wavelength.setCurrentWavelength




WaveLength Class
================

This an abstarct class that defines the main method for for all the wavelength classes. 


.. autoclass:: optics.wavelength.WaveLength
   :members:

Refractive Index Class
======================

Refatcive Index is a abstarct class that extents Wavelength and adds methods to calcuate the refarctive
index at the standard wavelnegths and the Abbe numners.

.. autoclass:: optics.wavelength.RefractiveIndex
   :members:

FixedIndex Class
================

Class to implement a basic fixed index which is independent of wavelength.

.. autoclass:: optics.wavelength.FixedIndex
   :members:

AirIndex Class
==============

Class to implement the refative index of air, this can either be fixed or calcualte by InfoIndex with Cauchy
coefficience. 

.. autoclass:: optics.wavelength.AirIndex
   :members:

The switch beween fixed and variable for the whole package. Note for large scale calcualtions, for example
simulation of an image it is computationally useful to switch to fixed index.
      
.. autofunction:: optics.wavelength.setFixedAirIndex

CauchyIndex Class
=================

Class to implement the simple Cauchy Index formula with either 

- the three Cauchy paramters, (a,b,c)
- the Nd and Vd values
- the 6 digit type integer XXXYYY where Nd = 1.XXX and Vd = Y.YY

.. autoclass:: optics.wavelength.CauchyIndex
   :members:

Sellmeier Class
===============

 Class to implement a simple two parameter Sellmeier index with only alpha and lambda_0 terms, it is
 inplemented as a simplied case of InfoIindex

 .. autoclass:: optics.wavelength.Sellmeier
    :members:

MaterialIndex Class
===================

Class to implement a materail refractive index where the paramteers are looked up in an internal database in the same
format used in RefrativeIndex.info. All the standard class type are includes, see MaterialData for details.

.. autoclass:: optics.wavelength.MaterialIndex
   :members:
		  
InfoIndex class
===============

Class to implement the refractive index calculations in the format used by RefrativeIndex.info site, this is mainly interval
class used by the MaterialIndex, AirIndex and CauchyIndex to do the actual calculations.

.. autoclass:: optics.wavelength.InfoIndex
   :members:

Spectrum Class
==============

Class to implment a constant spectrum independant of wavelength.

.. autoclass:: optics.wavelength.Spectrum
   :members:



GuassianSpectrum Class
======================

Class to implement a single peak Gaussian spectum

.. autoclass:: optics.wavelength.GaussianSpectrum
   :members:

PhotopicSpectrum Class
======================

Class to implement the Photopic (high brightness) spectral response of the eye. Modelled as a Guassian.

.. autoclass:: optics.wavelength.PhotopicSpectrum
   :members:


ScotopicSpectrum Class
======================

Class to implement the Scotopic (dark adapted) spectral response of the eye. Modelled as a Guassian.

.. autoclass:: optics.wavelength.ScotopicSpectrum
   :members:


TriColourSpectrum Class
=======================

Class to implement a three colour spectrum with peaks as Red, Green and Blue. All peaks are the same width.

.. autoclass:: optics.wavelength.TriColourSpectrum
	       :members:



PlanckSpectrum Class
=====================

Class to implemnt the Plank temperture dependand specturm.

.. autoclass:: optics.wavelength.PlanckSpectrum
   :members:


Colour Support Functions
========================

There are two support functions used in graphm these being

.. autofunction:: optics.wavelength.WavelengthColour



.. autofunction:: optics.wavelength.RefractiveIndexColour



