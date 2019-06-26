==============
Matrix Methods
==============

The set of classes in optics.matrix implement maxtrix analysis of optical
systems.

ParaxialMatrix class
====================

This class implements the basic paraxial matrix manipulations with the a
set of supporting methods


.. autoclass:: optics.matrix.ParaxialMatrix
   :members:

Extending Classes
-----------------

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
--------------

Matrix aglebra is implemented and the * operator that multiplies two matrices, so for example

   .. code-block:: python

      import optics.matrix as m
      a = m.ThinLensMatrix(100.0)
      b = m.ThinLensMatrix(30.0)
      c = a*b
      print("Combined focal length is : " + str(c.backFocalLength()))


This will create two lens one (a) of 100mm focal length, one (b) of 30mm and then form the paraxial matrix for
(a) followed by (b). The properties of that matrix can them be found with the method, for example focal length.

This can be extended for any number of matrices, for example if we have the same two lenses as above, but
separated by 20mm, then this can be be representated by:

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


this will create a unit matrix, the multiply by lens (a), distance of 20mm, then lens (b), so allowing the build up
of complex optical systems.

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

ParaxialGroup Class
===================

The ParaxialGroup class in an etension of ParaxialMatrix that add input plane, output plane and maxium height
at input and output. The input and output planes are defined in glabal coordinates, so in mm along the optical axis.
The methods to obtain the imaging planes are defedined so that they also return their position in global coordinates.

.. autoclass:: optics.matrix.ParaxialGroup
   :members:


Extending Classes
-----------------

There are a set of extending classes to implement specific ParaxialGroup making this erasiet to call.

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





