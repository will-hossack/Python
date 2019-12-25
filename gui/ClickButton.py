from PyQt5.QtWidgets import *

def cancel_button_clicked():
    alert = QMessageBox()
    print("cancel activated")
    alert.setText("Cancel Button pressed.")
    alert.exec()

def exit_button_clicked():
    alert = QMessageBox()
    alert.setText("Exit Button pressed.")
    alert.exec()

def main():

    app = QApplication([])
    window = QWidget()
    layout = QHBoxLayout()
    cancelbutton = QPushButton("Cancel")
    cancelbutton.clicked.connect(cancel_button_clicked)
    layout.addWidget(cancelbutton)
    exitbutton = QPushButton("Exit")
    exitbutton.clicked.connect(exit_button_clicked)
    layout.addWidget(cancelbutton)
    layout.addWidget(exitbutton)
    window.setLayout(layout)
    window.show()
    app.exec()

    

main()
