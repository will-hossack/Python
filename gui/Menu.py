import sys
from PyQt5.QtWidgets import *

class Example(QMainWindow):
    
    def __init__(self):
        super().__init__()
        
        self.initUI()
        
        
    def initUI(self):         
        
        self.statusbar = self.statusBar()
        statusMessage = QLabel("Hello")
        self.statusbar.addWidget(statusMessage)
        #self.statusbar.setNativeMenuBar(False)
        
        menubar = self.menuBar()
        menubar.setNativeMenuBar(False)
        viewMenu = menubar.addMenu('View')
        #menubar.setStatusTip("Menu Hover")
        
        viewStatAct = QAction('View statusbar', self, checkable=True)
        viewStatAct.setStatusTip('View statusbar')
        viewStatAct.setChecked(True)
        viewStatAct.triggered.connect(self.toggleMenu)
        
        viewMenu.addAction(viewStatAct)
        
        self.setGeometry(300, 300, 300, 200)
        self.setWindowTitle('Check menu')    
        self.show()
        
    def toggleMenu(self, state):
        
        if state:
            self.statusbar.show()
        else:
            self.statusbar.hide()
       
        
if __name__ == '__main__':
    
    app = QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())
    
