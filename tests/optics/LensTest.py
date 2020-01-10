from optics.lens import DataBaseLens,setCurrentLens,getCurrentLens

def main():

    newlens = DataBaseLens()

    print("Current is :" + str(getCurrentLens().getInfo()))
    setCurrentLens(newlens)
    print("New Lens is :" + str(getCurrentLens().getInfo()))

main()


