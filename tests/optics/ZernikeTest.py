"""      Test for Seidel wavefront
"""
from PyQt5.QtWidgets import QApplication
import optics.gui as gui
from optics.wavefront import ZernikeWaveFront

import tio as t
import matplotlib.pyplot as plt



def main():


    app  = QApplication([])
    w = gui.WaveFrontViewer()
    w.show()
    app.exec_()

main()
