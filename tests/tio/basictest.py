import tio

def main():

    y = tio.getFloat("Give a float")
    tio.tprint("Given value was : " + str(y))

    i = tio.getInt("Give and int")
    tio.tprint("int value is ",i)

    z = tio.getComplex("Give complex")
    tio.tprint("Complex value is :",z)

main()
