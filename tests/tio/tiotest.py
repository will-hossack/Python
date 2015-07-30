import tio
import math


def main():

   f = tio.openFile("Give file","r","lens")

   for l in f.readlines():
      tio.tprint(l)

   
   z = tio.getComplex("Give complex",complex(1,1))
   tio.tprint("Complex is " + repr(z))

   v = tio.getVector3d("Vector")
   tio.tprint("vector is : " + repr(v))
   
   options = ["start","close","quit","restart"]
   n,nopt = tio.getOption("Option",options)

   tio.tprint("Option {0:d} chosen, name is {1:s}".format(n,nopt))

   x = tio.getFloat("float",3.0,0.0,5.0)
   tio.tprint(x)

   st = tio.getString("and a string")
   tio.tprint(st)
   

main()
