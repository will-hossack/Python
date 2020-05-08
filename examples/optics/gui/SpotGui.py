
from PyQt5.QtWidgets import QApplication
from optics.gui import SpotViewer
from optics.lens import DataBaseLens



def main():
    app = QApplication([])
    lens = DataBaseLens("$LENS/Tessar-F6.3")
    ex = SpotViewer(lens) #LensViewer(lens)
    ex.show()
    app.exec_()

main()
