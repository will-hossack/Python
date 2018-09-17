import tio as t

def main():

    while t.getBool("Continue",True):
        z = t.getComplex("Give a complex")
        t.tprint("Complex is : ",z)

main()
