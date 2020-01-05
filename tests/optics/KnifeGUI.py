from PyQt5.QtWidgets import *
from PyQt5.QtWidgets import QApplication

import optics.lens as l
from optics.gui import *


def main():
    app = QApplication([])
    lens = l.SimpleSinglet()
    ex = KnifeViewer(lens) #LensViewer(lens)
    ex.show()
    app.exec_()

main()
