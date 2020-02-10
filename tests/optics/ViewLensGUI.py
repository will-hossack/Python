from PyQt5.QtWidgets import *
from PyQt5.QtWidgets import QApplication

import optics.lens as l
from optics.gui import *


def main():
    app = QApplication([])
    lens = l.DataBaseLens("$LENS/Tessar-F2.8")
    ex = LensViewer()
    ex.show()
    app.exec_()

main()
