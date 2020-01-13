import optics.wavelength as w
from optics.gui import *
from PyQt5.QtWidgets import *


def main():
    app = QApplication([])
    p = PlaneSetter()
    p.show()
    app.exec()


main()
