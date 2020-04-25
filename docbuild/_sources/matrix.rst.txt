==============
Matrix Methods
==============

The set of classes in optics.matrix implement maxtrix analysis of optical
systems.

ParaxialMatrix Class
====================

This class implements the basic paraxial matrix manipulations with the a
set of supporting methods to extract various optical parameters.


.. autoclass:: optics.matrix.ParaxialMatrix
   :members:

Extending Classes
=================

There are a set of extending classes that implement the various specific
ParaxialMatrices, these being

.. autoclass:: optics.matrix.PropagationMatrix
   :members:

.. autoclass:: optics.matrix.DielectricMatrix
   :members:

.. autoclass:: optics.matrix.ThinLensMatrix
   :members:

.. autoclass:: optics.matrix.ThickLensMatrix
   :members:

.. autoclass:: optics.matrix.DoubletMatrix
   :members:

.. autoclass:: optics.matrix.MirrorMatrix
   :members:


.. autoclass:: optics.matrix.CavityMatrix
   :members:


Matrix Algebra
==============

Matrix aglebra is implemented and the * operator that multiplies two matrices, so for example

   .. code-block:: python

      import optics.matrix as m
      a = m.ThinLensMatrix(100.0)
      b = m.ThinLensMatrix(30.0)
      c = a*b
      print("Combined focal length is : " + str(c.backFocalLength()))


This will create two lens one (a) of 100mm focal length, one (b) of
30mm and then form the paraxial matrix for (a) followed by (b). The properties of that
matrix can them be found with the method, for example focal length.

This can be extended for any number of matrices, for example
if we have the same two lenses as above, but separated by 20mm, then this can be be representated by:

     .. code-block:: python

      import optics.matrix as m
      a = m.ThinLensMatrix(100.0)
      d = m.PropagationMatrix(20.0)
      b = m.ThinLensMatrix(30.0)
      c = a*d*b
      print("Combined focal length is : " + str(c.backFocalLength()))

Note that matrix multiplication for NOT commute, so the order of the operants is significant.


The \*= operator implements the muiltiplication by, (in situe) so for example,

    .. code-block:: python

       import optics.matrix as m
       a = m.ThinLensMatrix(100.0)
       d = m.PropagationMatrix(20.0)
       b = m.ThinLensMatrix(30.0)
       s = m.ParaxialMatrix()
       s *= a                   # Through first lens
       s *= d                   # propagate distance 20mm
       s *= b                   # Through second lens
       print("Combined focal length is : " + str(s.backFocalLength()))


this will create a unit matrix, the multiply by lens (a),
distance of 20mm, then lens (b), so allowing the build up of complex optical systems.

Propagation a specified distance is such a commom operation that it is also implementred by the + and \+= operators
so that

    .. code-block:: python

       import optics.matrix as m
       a = m.ThinLensMatrix(100.0)
       b = m.ThinLensMatrix(30.0)
       s = m.ParaxialMatrix()
       s *= a                  # Though first lens
       s += 20.0               # Propagate 20 mm
       s *= b                  # Though second lens
       print("Combined focal length is : " + str(s.backFocalLength()))
      
repeats the above calculation with the \+= implementating a propagation of 20 mm.

Note that with the basic ParaxialMatrix classes focal and principle plane location are measure from the input and output
planes. 

ParaxialGroup Class
===================

The ParaxialGroup class in an etension of ParaxialMatrix that add input plane, output plane and maxium height
at input and output. The input and output planes are defined in global coordinates, so in mm along the optical axis.
The methods to obtain the imaging planes are defedined so that they also return their position in global coordinates.

.. autoclass:: optics.matrix.ParaxialGroup
   :members:


Extending Classes
=================

There are a set of extending classes to implement specific ParaxialGroup making these easier to call.

.. autoclass:: optics.matrix.ParaxialThinLens
   :members:


.. autoclass:: optics.matrix.ParaxialThickLens
   :members:


.. autoclass:: optics.matrix.ParaxialDoublet
   :members:


.. autoclass:: optics.matrix.ParaxialMirror
   :members:

.. autoclass:: optics.matrix.DataBaseMatrix
   :members:

.. autoclass:: optics.matrix.ParaxialPlane
   :members:

 
Example Code
============

The first example read the ParaxialGroup imformation from a lens from the lens data base, scale it focal length
to 40 mm, then calculates the object and images planes to give a magnification of -0.1.

.. code-block:: python

   import optics.matrix as m
   lens = m.DataBaseMatrix("$LENS/Tessar-100")    # Inport lens
   lens.setFocalLength(40)                        # Set focal length
   obj,ima = lens.planePair(300,-0.1)             # Form pair of planes
   print("Object plane at : " + str(obj.inputPlane()))
   print("Image plane at : " + str(ima.inputPlane())))


ParaxialGroup Algebra
=====================
 
The matrix algebra with the ParaxialGroups is a little more complex than above, there are two cases, the first is where
the ParaxialGroup is multiples by a ParaxialMatrix, here the ParaxialMatrix operates on the underlying matrix in the
ParaxalGroup, so for example

.. code-block:: python

   import optics.matrix as m
   lens = m.ParxialThinLens(40.0,100.0,radius=10.0)   # Thins lens
   lens += 10.0                         # Add a 10 mm propagation
   lens *= m.ThinLensMatrix(80.0)       # Add a 80 mm thin lens

This will create ParaxialGroup for two thin lenses serarated by 10mm. Note the second aregument cannot be a ParaxialGroup.

The second higher level is multiplying ParaxialGroups together, so for example the following piece of code.

.. code-block:: python
		
   import optics.matrix as m
   first = m.ParxialThinLens(40.0,100.0,radius=10.0)   # Thins lens at 40mm focal length 100 mm
   second = m.PaxaxialThinLens(60,-80.0,radius=10.0)   # Multiply by a lines at 60 mm focal lenght -80 mm
   system = first * second
   print("Focal length of system is : " + str(system.backFocalLength()))
 
 
 The First line creates a thin lens of focal length 100mm, radius 10 located at 40mm along the optical axis. The second
 a thin lens of focal length -80mm, radius 10, located at 60mm along the optical axis. These are then combined
 to form a single ParaxialGroup being the first lens followed by the second, so it will auotmatically put in the
 parpagation between the lenses. This allows the paraxial analysis of compound systems with multiple lenses.

 Note: the ParaxialGroups are assumed to be located along the optical axis so they do not overlap with the distance between
 then being positive. If this is not true then the calculations all work, but it will be the analysis of a physically impossible
 optical system.

Graphics
========

A garphical plot of the planes can be implemented in MatPlotLib with the .draw() method as show below.

.. code-block:: python

   import optics.matrix as m
   import matplotlib.pyplot as plt
   lens = m.ParxialThinLens(40.0,100.0,radius=10.0)   # Thins lens
   lens.draw(True)          # Draw plans with a ledgend
   plt.show()               # Display the plot

 This will create the PaxaialGroup of a thin lens of focal length 100mm located at 40mm along the optical axis, and then
 plot the the location of the input, output, focal and principle planes with an optional ledgend in the lower left of the
 plot.
 
