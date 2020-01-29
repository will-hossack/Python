"""      Test for Seidel wavefront
"""
from PyQt5.QtWidgets import QApplication
import optics.gui as gui




def main():


    app  = QApplication([])
    w = gui.OTFViewer()
    w.show()
    app.exec_()

main()
