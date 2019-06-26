=================
Wavelength Module
=================

This module deal with the wavelnegth dependand features in particular refractive and
intensity specstrums.

Wavelength Value
=================

There are a set of globals defined covering the standard wavelength. All are in microns.

.. automodule:: optics.wavelength
   :members: Red,Green,Blue,BlueLimit,RedLimit,BlueColourMatch,GreenColourMatch,RedColourMatch,Mercury_i,Mercury_h,Cadmium_F,Hydrogen_F


Default Wavelength
==================

The package default wavelength is controlled by two

.. autofunction:: optics.wavelength.getDefaultWavelength

which is automatically called in stareted up with the result being places in globar variable.

- optics.wavelength.Default

The package default can be changed with a call to

.. autofunction:: optics.wavelength.setDefaultWavelength

which updates the global variable.

Note that any objects created priod to this call will NOT be updated.


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
		  
