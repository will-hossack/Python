import tio as t
import cmath

def main():


    s = t.getString("Give a string") 
    t.tprint("Given string was : ",s)

    y = t.getFloat("Give a float",min=None,max=100)
    t.tprint("Given value was : ",y)

    i = t.getInt("Give and int",default = 30,max = None)
    t.tprint("int value is ",i)

    z = t.getComplex("Give complex",cmath.sqrt(-6))
    t.tprint("Complex value is :",z)


    filename = t.getFilename("File Name","pgm","image")
    t.tprint(filename)

    file = t.openFile("File","w","data","output")
    file.write("Hello\n")
    file.close()
main()
