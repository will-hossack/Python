"""
Example program to use gui with PyQt5 to display a lens with collimated beam
"""

from PyQt5.QtWidgets import QApplication
from optics.gui import LensViewer
from optics.lens import DataBaseLens


def main():
    
    app = QApplication([])                    # Initialse PyQt5 
    lens = DataBaseLens("$LENS/Tessar-F2.8")  # Get a test lens
    ex = LensViewer(lens)                     # Make viewer
    ex.show()                               # Make viewer active
    app.exec()                              # run PyQt app

main()
