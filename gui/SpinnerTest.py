from PyQt5.QtWidgets import *

class spindemo(QWidget):
    def __init__(self, parent = None):
      super(spindemo, self).__init__(parent)

      layout = QVBoxLayout()
      self.lab = QLabel("Current Value : ")
      layout.addWidget(self.lab)
      self.spin = QDoubleSpinBox()
      self.spin.setValue(0.55)
      self.spin.setSingleStep(0.01)
      self.spin.setRange(0.3,0.8)
      layout.addWidget(self.spin)
      self.spin.valueChanged.connect(self.valuechange)
      exitbutton = QPushButton("Exit")
      exitbutton.clicked.connect(self.exit_button_clicked)
      layout.addWidget(exitbutton)
      self.setLayout(layout)
      self.setWindowTitle("SpinBox demo")
      self.show()

    def valuechange(self):
        self.lab.setText("Current Value : " + "{0:4.2f}".format(self.spin.value()))

    def exit_button_clicked(self):
        self.close()


def main():
    app = QApplication([])
    ex = spindemo()

    app.exec()

main()

