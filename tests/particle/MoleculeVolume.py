import tio as t
import vector as c
import particle as p
import math

def main():

    r = 3.5
    part = p.Particle(c.Vector3d(1,2,3),radius = r)
    t.tprint("Initial Partile :", part)

    bb = p.BoundingBox(part)
    t.tprint("Box : ",bb)

    vol = bb.volume()
    t.tprint("Volume is : ",vol)

    maxPoint = 2000000
    i = 0
    inside = 0
    while i < maxPoint:
        pt = bb.randomPoint()
        if part.isInside(pt):
            inside += 1
        i += 1

    t.tprint("Number inside is ", inside)
    pv = vol*inside/maxPoint
    tv = 4.0*math.pi*r**3/3.0
    t.tprint("True vol is : ",tv, " Calcualted is : ",pv)


main()
