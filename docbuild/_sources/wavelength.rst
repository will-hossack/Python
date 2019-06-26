=================
Wavelength Module
=================

This module deal with the wavelnegth dependand features in particular refractive and
intensity specstrums.

Wavelength Value
=================

.. automodule:: optics.wavelength
   :members: Red,Green,Blue,BlueLimit,RedLimit,BlueColourMatch,GreenColourMatch,RedColourMatch
   :members: ScotopicPeak,ScotopicWidth,PhotopicPeak,PhotopicWidth


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
		  
