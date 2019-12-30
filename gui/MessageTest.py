from PyQt5.QtWidgets import *
import optics.wavefront as w
import optics.lens as l
import optics.gui as gui

class MessageBox(QMessageBox):
    """
    Class to set the Default and Design wavelengths with spinners
    """
    def __init__(self, text = "None",dtext = "None", parent = None):
        super(MessageBox,self).__init__(parent)

        self.setIcon(QMessageBox.Critical)
        self.setText(text)
        self.setWindowTitle("Message")
        self.setDetailedText(dtext)
        self.setStandardButtons(QMessageBox.Ok )

def main():

    app = QApplication([])
    lens = l.DataBaseLens("$LENS/Tessar-F4.5")
    gui.CurrentLens = lens
    z = w.ZernikeWaveFront(1.0,0.65,1.0,2.0,-3.0,-4.5,6.45)
    s = MessageBox("Lens Components",lens.getInfo())
    s.resize(500,300)
    s.show()
    app.exec()

    

main()
