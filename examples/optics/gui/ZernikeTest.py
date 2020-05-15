"""      Test for fitting wavefronts
"""
from PyQt5.QtWidgets import QApplication
import optics.gui as gui



def main():


    app  = QApplication([])
    w = gui.WaveFrontViewer()
    w.show()
    app.exec_()

main()
