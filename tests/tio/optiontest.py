import tio

def main():
    options = ["stop","go","reset","reformat"]
    i,opt = tio.getOption("Which option",options,1)
    tio.tprint("Option chosen was : ",i," being ",opt)

main()
