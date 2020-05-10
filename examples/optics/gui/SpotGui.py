"""
Example program to read in a lens from DataBase and display an interactiive stop diagram for
a collimated beam
"""
from PyQt5.QtWidgets import QApplication
from optics.gui import SpotViewer
from optics.lens import DataBaseLens



def main():
    app = QApplication([])
    lens = DataBaseLens("$LENS/Tessar-F6.3")
    ex = SpotViewer(lens) 
    ex.show()
    app.exec_()

main()
