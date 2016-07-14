import tio as t
import vector as v

def main():
    
    u = v.Unit3d().random()

    vec = v.Vector3d().random(10.0)


    t.tprint("Random vector is : " + str(vec))

main()
