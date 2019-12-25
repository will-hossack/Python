import optics.wavelength as w
from optics.gui import *
from PyQt5.QtWidgets import *


def main():
    app = QApplication([])
    dire = DirectionSetter()
    dire.show()
    app.exec()
    print(str(CurrentAngle))

main()
