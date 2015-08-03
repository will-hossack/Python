import optics.surface as sur
import optics.ray as ray
import tio
import math

def main():

    s = sur.OpticalPlane(4)
    tio.tprint(repr(s))

    r = ray.IntensityRay([3,4,0],math.radians(10))
    tio.tprint(repr(r))

    inter = s.getSurfaceInteraction(r)

    tio.tprint(repr(inter))

main()
