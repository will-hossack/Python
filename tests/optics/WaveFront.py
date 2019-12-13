import optics.wavefront as w
import optics.zernike as z
import vector as v

def main():

    pt = v.Vector2d(0.5,-0.6)
    wa = w.WavePoint(pt)
    print(repr(wa))
    ze = z.ZernikeExpansion(2.0,0.55,1.0,3.0,-4.0)
    print(repr(ze))

    wa.setWithZernike(ze)

    print(repr(wa))

main()
