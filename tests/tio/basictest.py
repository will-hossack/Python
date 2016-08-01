import tio as t
import cmath

def main():


    s = t.getString("Give a string") 
    t.tprint("Given string was : ",s)

    y = t.getFloat("Give a float",max=100)
    t.tprint("Given value was : ",y)

    i = t.getInt("Give and int",default = 30)
    t.tprint("int value is ",i)

    z = t.getComplex("Give complex",cmath.sqrt(-6))
    t.tprint("Complex value is :",z)


    file = t.openFile("File","w","data","output")
    file.write("Hello\n")
    file.close()
main()
