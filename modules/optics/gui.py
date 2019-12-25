"""
     Set of gui classes and methods
"""

import optics.wavelength as w
import optics.ray as ray
import vector as v
from PyQt5.QtWidgets import *

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt


CurrentLens = None
CurrentAngle = v.Unit3d(0.0,0.0,1.0)

class WaveLengthSetter(QWidget):
    """
    Class to set the Default and Design wavelengths with spinners
    """
    def __init__(self, parent = None):
        super(WaveLengthSetter,self).__init__(parent)

        layout = QGridLayout()     # Vertical box
        defaultLabel = QLabel("Default : ")
        designLabel = QLabel("Design : ")
        
        self.defaultSpin = QDoubleSpinBox()      # Default wave spinner
        self.defaultSpin.setValue(w.Default)
        self.defaultSpin.setSingleStep(0.01)
        self.defaultSpin.setRange(w.BlueLimit,w.RedLimit)
        self.defaultSpin.valueChanged.connect(self.defaultValueChange)

        self.designSpin = QDoubleSpinBox()       # Design wave spinner
        self.designSpin.setValue(w.Design)
        self.designSpin.setSingleStep(0.01)
        self.designSpin.setRange(w.BlueLimit,w.RedLimit)
        self.designSpin.valueChanged.connect(self.designValueChange)
        
        closeButton = QPushButton("Close")      # The close button
        closeButton.clicked.connect(self.closeButtonClicked)
        
        layout.addWidget(defaultLabel,0,0)      # Add the 5 item in Grid
        layout.addWidget(self.defaultSpin,0,1)
        layout.addWidget(designLabel,1,0)
        layout.addWidget(self.designSpin,1,1)
        layout.addWidget(closeButton,2,1)

        self.setLayout(layout)
        self.setWindowTitle("WaveLengthSetter")

    """
    Method to update the values from the spinners and close the widow
    """
        
    def defaultValueChange(self):
        v = self.defaultSpin.value()
        w.setDefaultWavelength(v)

    def designValueChange(self):
        v = self.designSpin.value()
        w.setDesignWavelength(v)
        
    def closeButtonClicked(self):      # Close the frame
        self.close()
        



class DirectionSetter(QWidget):
    """
    Class to set default direction in degress with spinners
    """
    def __init__(self, parent = None):
        super(DirectionSetter,self).__init__(parent)
        self.theta = 0.0
        self.psi = 0.0
        CurrentAngle.setPolarDegrees(self.theta,self.psi)
        layout = QGridLayout()     # Vertical box
        thetaLabel = QLabel("Theta in degress :")
        psiLabel = QLabel("Psi in degrees : ")
        self.unitLabel = QLabel(str(CurrentAngle))
        
        self.thetaSpin = QDoubleSpinBox()       # Theta wave spinner
        self.thetaSpin.setValue(self.theta)
        self.thetaSpin.setSingleStep(1.0)
        self.thetaSpin.setRange(0.0,90.0)
        self.thetaSpin.valueChanged.connect(self.valueChange)

        self.psiSpin = QDoubleSpinBox()       # Theta wave spinner
        self.psiSpin.setValue(self.theta)
        self.psiSpin.setSingleStep(1.0)
        self.psiSpin.setRange(-90.0,90.0)
        self.psiSpin.valueChanged.connect(self.valueChange)


        resetButton = QPushButton("Reset")
        resetButton.clicked.connect(self.resetButtonClicked)
        closeButton = QPushButton("Close")      # The close button
        closeButton.clicked.connect(self.closeButtonClicked)
        
        layout.addWidget(thetaLabel,0,0)
        layout.addWidget(self.thetaSpin,0,1)
        layout.addWidget(psiLabel,1,0)
        layout.addWidget(self.psiSpin,1,1)
        layout.addWidget(self.unitLabel,2,0)
        layout.addWidget(resetButton,3,0)
        layout.addWidget(closeButton,3,1)
        self.setLayout(layout)
        self.setWindowTitle("DirectionSetter")


    """
    Method to update the values from the spinners and close the widow
    """
        
    def valueChange(self):
        self.theta = self.thetaSpin.value()
        self.psi = self.psiSpin.value()
        CurrentAngle.setPolarDegrees(self.theta,self.psi)
        self.unitLabel.setText(str(CurrentAngle))        

    def resetButtonClicked(self):
        self.thetaSpin.setValue(0.0)
        self.psiSpin.setValue(0.0)
        self.valueChange()
        
    def closeButtonClicked(self):      # Close the frame
        self.close()




class LensViewer(QWidget):
    def __init__(self, lens = None , parent=None):
        super(LensViewer, self).__init__(parent)

        if lens == None:
            self.lens = CurrentLens
        else:
            self.lens = lens

        #     Setup components
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        toolbar = NavigationToolbar(self.canvas, self)

        plotButton = QPushButton('Plot')     # Plot button
        plotButton.clicked.connect(self.plot)
        waveButton = QPushButton("Change Wavelength")
        waveButton.clicked.connect(self.waveButtonClicked)
        angleButton = QPushButton("Set Angle")
        angleButton.clicked.connect(self.angleButtonClicked)
        


        #                  Set the layout
        layout = QVBoxLayout()
        layout.addWidget(toolbar)
        layout.addWidget(self.canvas)
        layout.addWidget(plotButton)
        layout.addWidget(waveButton)
        layout.addWidget(angleButton)
        self.setLayout(layout)


    def waveButtonClicked(self):
        """
        Wavelength setter
        """
        ws = WaveLengthSetter(parent=self)
        ws.move(50,50)
        ws.resize(200,100)
        ws.show()

    def angleButtonClicked(self):
        aset = DirectionSetter(parent=self)
        aset.move(50,100)
        aset.resize(200,200)
        aset.show()
        
    def plot(self):
        #                Do the plot
        self.figure.clear()
        panel = self.figure.add_subplot(111)
        panel.axis('equal')

        u = CurrentAngle
        pencil = ray.RayPencil().addCollimatedBeam(self.lens,u,"vl",\
                                                   wave=w.Default).addMonitor(ray.RayPath())

        #
        #        Set the output plane (being the back focal plane)
        op = self.lens.backFocalPlane()

        #         Propagate pencil through lens and one to back plane
        pencil *= self.lens       # Through lens
        pencil *= op         # To plane

        
        # plot data
        self.lens.draw()
        op.draw()
        pencil.draw()
        plt.grid()
        plt.xlabel("Optical Axis")
        plt.ylabel("Height")
        plt.title("Diagram of lens " + self.lens.title)
        # refresh canvas
        self.canvas.draw()
