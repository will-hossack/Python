import optics.wavelength as w
from optics.gui import *
from PyQt5.QtWidgets import *


def main():
    print("Default : " + str(w.Default) + " Deign : " + str(w.Design))
    app = QApplication([])
    wave = WaveLengthSetter()
    wave.show()
    print("Default : " + str(w.Default) + " Deign : " + str(w.Design))
    app.exec()
    print("Default : " + str(w.Default) + " Deign : " + str(w.Design))

main()
