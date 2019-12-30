import optics.wavelength as w
import optics.gui as gui
from PyQt5.QtWidgets import *


def main():
    app = QApplication([])
    zw = gui.ZernikeOrderSetter()
    zw.show()
    
    app.exec()
    print("Order : " + str(gui.ZernikeOrder))
    print("Iris : " + str(gui.IrisRatio))

main()
