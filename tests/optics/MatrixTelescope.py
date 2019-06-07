"""
   Form a paraxial telescope with DataBase lenses

"""

import optics.matrix as m
import matplotlib.pyplot as plt
import tio as t


def main() :

   objective = m.DataBaseMatrix("$LENS/Doublet-100")
   objective.setFocalLength(250.0)
   t.tprint(objective.getInfo())

   eyepiece = m.DataBaseMatrix("$LENS/HygenEyepiece-100")
   eyepiece.setFocalLength(20.0)
   t.tprint(eyepiece.getInfo())

   bfp = objective.backFocalPlane()
   ffp = eyepiece.frontFocalPlane()

   eyepiece.input_plane += (bfp-ffp)

   telescope = m.ParaxialSystem(objective,eyepiece)

   
   fig,ax = plt.subplots()        # Get the axis
   telescope.draw()
   ax.set_aspect(1.0)             # set unit aspect
   plt.show()
    




main()
