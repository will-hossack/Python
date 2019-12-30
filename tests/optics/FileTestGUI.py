import optics.wavelength as w
from optics.gui import *
from PyQt5.QtWidgets import *


def main():

        
    app = QApplication([])
    ls = LensSetter()
    ls.show()
    ls.close()
    
    app.exec()
    

main()
