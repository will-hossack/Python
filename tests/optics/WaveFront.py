import optics.analysis as a
import optics.zernike as z
import vector as v

def main():

    pt = v.Vector2d(0.5,-0.6)
    w = a.WavePoint(pt)
    print(repr(w))
    ze = z.ZernikeExpansion(2.0,0.55,1.0,3.0,-4.0)
    print(repr(ze))

    w.setWithZernike(ze)

    print(repr(w))

main()
