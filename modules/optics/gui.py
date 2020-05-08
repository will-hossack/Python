"""
     Set of gui classes and methods
"""
from os import getenv
import sys
import math
from optics.wavelength import getCurrentWavelength,setCurrentWavelength,\
    getDesignWavelength,setDesignWavelength,getDefaultWavelength,BlueLimit,RedLimit
from optics.lens import setCurrentLens,getCurrentLens
from optics.wavefront import WaveFrontAnalysis,Interferometer
from optics.analysis import KnifeTest
from optics.ray import RayPencil,RayPath,getCurrentAngle,setCurrentAngle
from optics.psf import Psf,SpotDiagram,getReferencePointOption,setReferencePointOption,\
        getPlaneShift,setPlaneShift,incrementPlaneShift
from optics.surface import OpticalPlane
from vector import Unit3d
from PyQt5.QtWidgets import QWidget,QLabel,QDoubleSpinBox,QPushButton,QGridLayout,\
    QDial,QMessageBox,QMainWindow,QAction,QVBoxLayout,QFileDialog,QRadioButton
from PyQt5.QtCore import Qt

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt



"""
Globals to control the various plots, these are set by the various buttons
"""
PanelSize = (800,600)
IrisRatio = 1.0
ZernikeOrder = 4
CurrentWaveFront = None
Xtilt = 3.0
Ytilt = 0.0
CurrentKnife = 0.0
CurrentKnifeAngle = 0.0
CurrentKnifeShift = 0.0
CurrentWire = False


def getGlobals():
    return ZernikeOrder

class WaveLengthSetter(QWidget):
    """
    Class to set the Current and Design wavelengths with spinners. 

    :param parent: the calling frame of widget
    :param closeAction: function to be executed when panel closed.
    """
    def __init__(self, parent = None,closeAction = None):
        super(WaveLengthSetter,self).__init__(parent)

        self.closeAction = closeAction
        self.setAutoFillBackground(True)
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.white)
        self.setPalette(p)

        #      The default wavelength spinner
        currentLabel = QLabel("Current : ")
        self.currentSpin = QDoubleSpinBox()      # Default wave spinner
        self.currentSpin.setValue(getCurrentWavelength())
        self.currentSpin.setSingleStep(0.01)
        self.currentSpin.setRange(BlueLimit,RedLimit)
        self.currentSpin.valueChanged.connect(self.currentValueChange)

        #       The design wavelength spinner
        designLabel = QLabel("Design : ")
        self.designSpin = QDoubleSpinBox()       # Design wave spinner
        self.designSpin.setValue(getDesignWavelength())
        self.designSpin.setSingleStep(0.01)
        self.designSpin.setRange(BlueLimit,RedLimit)
        self.designSpin.valueChanged.connect(self.designValueChange)

        #       The close and rest buttons
        closeButton = QPushButton("Close")      # The close button
        closeButton.clicked.connect(self.closeButtonClicked)
        resetButton = QPushButton("Reset")
        resetButton.clicked.connect(self.resetButtonClicked)

        #      Use Grid layout 
        layout = QGridLayout()
        layout.addWidget(currentLabel,0,0)      # Add the 6 item in Grid
        layout.addWidget(self.currentSpin,0,1)
        layout.addWidget(designLabel,1,0)
        layout.addWidget(self.designSpin,1,1)
        layout.addWidget(resetButton,2,0)
        layout.addWidget(closeButton,2,1)
        self.setLayout(layout)
        self.setWindowTitle("Wavelength Setter")


        # Method to update the values from the spinners, close and rest buttons
        
    def currentValueChange(self):
        v = self.currentSpin.value()
        setCurrentWavelength(v)

    def designValueChange(self):
        v = self.designSpin.value()
        setDesignWavelength(v)

    def resetButtonClicked(self):     # Reset both wavelengths to Green
        c = getDefaultWavelength()
        self.currentSpin.setValue(c)
        self.designSpin.setValue(c)
        setCurrentWavelength(c)
        setDesignWavelength(c)
        
    def closeButtonClicked(self):     # Close and execute close action if given
        self.close()
        if self.closeAction != None:
            self.closeAction()
        

class DirectionSetter(QWidget):
    """
    Class to set default direction in degress with spinners, note the actual direction is help as Unit3d.

    :param parent: the calling frame
    :param closeAction: function to execute on clsoing the window
    """
    def __init__(self, parent = None,closeAction = None):
        super(DirectionSetter,self).__init__(parent)

        self.closeAction = closeAction
        
        self.setAutoFillBackground(True)
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.white)
        self.setPalette(p)

        #      Set up the lables.
        thetaLabel = QLabel("Theta in degress :")     # The labels
        psiLabel = QLabel("Psi in degrees : ")
        self.unitLabel = QLabel(str(getCurrentAngle()))    # Label in Unit3d format

        self.theta,self.psi = getCurrentAngle().getAngle().getDegrees()  # Current Theta/Psi
        self.thetaSpin = QDoubleSpinBox()       # Theta spinner
        self.thetaSpin.setRange(-90.0,90.0)
        self.thetaSpin.setValue(self.theta)
        self.thetaSpin.setSingleStep(1.0)
        self.thetaSpin.valueChanged.connect(self.valueChange)

        self.psiSpin = QDoubleSpinBox()       # Psi spinner
        self.psiSpin.setRange(-180.0,180.0)
        self.psiSpin.setValue(self.psi)
        self.psiSpin.setSingleStep(1.0)
        self.psiSpin.valueChanged.connect(self.valueChange)


        resetButton = QPushButton("Reset")      # Reset and close buttons
        resetButton.clicked.connect(self.resetButtonClicked)
        closeButton = QPushButton("Close")
        closeButton.clicked.connect(self.closeButtonClicked)

        layout = QGridLayout()                   # Set up layout as grid
        layout.addWidget(thetaLabel,0,0)
        layout.addWidget(self.thetaSpin,0,1)
        layout.addWidget(psiLabel,1,0)
        layout.addWidget(self.psiSpin,1,1)
        layout.addWidget(self.unitLabel,2,0,2,2)
        layout.addWidget(resetButton,3,0)
        layout.addWidget(closeButton,3,1)
        self.setLayout(layout)
        self.setWindowTitle("DirectionSetter")
        

    """
    Method to update the values from the spinners and close the widow
    """
        
    def valueChange(self):           # For either spinned
        self.theta = self.thetaSpin.value()
        self.psi = self.psiSpin.value()
        u = Unit3d().setPolarDegrees(self.theta,self.psi)
        setCurrentAngle(u)
        self.unitLabel.setText(str(getCurrentAngle()))        

    def resetButtonClicked(self):
        self.thetaSpin.setValue(0.0)
        self.psiSpin.setValue(0.0)
        self.valueChange()
        
    def closeButtonClicked(self):      # Close the frame
        self.close()
        if self.closeAction != None:
            self.closeAction()


class IrisSetter(QWidget):
    """
    Class to set default direction in degress with spinners
    """
    def __init__(self, parent = None,closeAction = None):
        super(IrisSetter,self).__init__(parent)

        global IrisRatio
        
        self.closeAction = closeAction
        self.setAutoFillBackground(True)
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.white)
        self.setPalette(p)

        
        self.irisRatio = IrisRatio
        self.irisLabel = QLabel("Iris : " + str(self.irisRatio))
        self.irisDial = QDial()
        self.irisDial.setRange(0,100)
        self.irisDial.setValue(int(100*self.irisRatio))
        self.irisDial.valueChanged.connect(self.valueChange)
        closeButton = QPushButton("Close")
        closeButton.clicked.connect(self.closeButtonClicked)
        resetButton = QPushButton("Reset")
        resetButton.clicked.connect(self.resetButtonClicked)
        

        layout = QGridLayout()
        layout.addWidget(self.irisLabel,0,1)
        layout.addWidget(self.irisDial,1,0,1,2)
        layout.addWidget(resetButton,2,0)
        layout.addWidget(closeButton,2,1)
        self.setLayout(layout)
    
    def valueChange(self):
        global IrisRatio
        self.irisRatio = self.irisDial.value()/100.0
        self.irisLabel.setText("Iris : " + str(self.irisRatio))
        IrisRatio = self.irisRatio
        getCurrentLens().setIris(self.irisRatio)


    def resetButtonClicked(self):
        self.iris = 1.0
        self.irisDial.setValue(100)
        self.valueChange()
        
    def closeButtonClicked(self):      # Close the frame
        self.close()
        if self.closeAction != None:
            self.closeAction()



class PlaneSetter(QWidget):
    """
    Class for plane setter with dial, this update using setPlaneShift and
    incrementPlaneShift
    
    """
    def __init__(self, scale = 0.01 , parent = None,closeAction = None,changedAction = None):
        super(PlaneSetter,self).__init__(parent)
        
        self.closeAction = closeAction
        self.changedAction = changedAction
        self.setAutoFillBackground(True)
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.white)
        self.setPalette(p)

        self.shift = getPlaneShift()
        self.scale = scale
        
        self.planeLabel = QLabel("Plane : " + "{0:5.3f}".format(self.shift))
        self.planeDial = QDial()
        self.planeDial.setRange(-100,100)
        self.planeDial.setValue(int(self.shift/self.scale))
        self.planeDial.valueChanged.connect(self.valueChange)
        closeButton = QPushButton("Close")
        closeButton.clicked.connect(self.closeButtonClicked)
        resetButton = QPushButton("Reset")
        resetButton.clicked.connect(self.resetButtonClicked)
        

        layout = QGridLayout()
        layout.addWidget(self.planeLabel,0,1)
        layout.addWidget(self.planeDial,1,0,1,2)
        layout.addWidget(resetButton,2,0)
        layout.addWidget(closeButton,2,1)
        self.setLayout(layout)
    
    def valueChange(self):
        global PlaneShift
        self.shift = self.planeDial.value()*self.scale
        self.planeLabel.setText("Plane : " + "{0:5.3f}".format(self.shift))
        setPlaneShift(self.shift)
        if self.changedAction != None:
            self.changedAction()


    def resetButtonClicked(self):
        self.shift = 0.0
        self.planeDial.setValue(0)
        self.valueChange()
        
    def closeButtonClicked(self):      # Close the frame
        self.close()
        if self.closeAction != None:
            self.closeAction()


class LensSetter(QWidget):
    """
    Class to set the current lens with a dialogue box.
    """
    def __init__(self, parent = None,closeAction = None):
        super(LensSetter,self).__init__(parent)

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        dir = getenv("LENS")
        fileName, _ = QFileDialog.getOpenFileName(self,"Lens Files",dir,\
                                                  "Lens Files (*.lens)", options=options)
        if fileName:
            setCurrentLens(fileName)

        if closeAction != None:      # Close the frame
            closeAction()
            
        self.close()
        

class ZernikeOrderSetter(QWidget):
    """
    Widget to set the Zernike Expansion order.
    """
    def __init__(self, parent = None,closeAction = None):
        super(ZernikeOrderSetter,self).__init__(parent)

        global ZernikeOrder

        self.closeAction = closeAction
        self.setAutoFillBackground(True)
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.white)
        self.setPalette(p)

        self.zernikeFourButton = QRadioButton("Fourth")
        if ZernikeOrder == 4:
            self.zernikeFourButton.setChecked(True)
        self.zernikeFourButton.order = 4
        self.zernikeFourButton.clicked.connect(lambda: self.buttonClicked(self.zernikeFourButton))
        self.zernikeSixButton = QRadioButton("Sixth")
        if ZernikeOrder == 6:
            self.zernikeSixButton.setChecked(True)
        self.zernikeSixButton.order = 6
        self.zernikeSixButton.clicked.connect(lambda: self.buttonClicked(self.zernikeSixButton))
        self.zernikeEightButton = QRadioButton("Eighth")
        if ZernikeOrder == 8:
            self.zernikeEightButton.setChecked(True)
        self.zernikeEightButton.order = 8
        self.zernikeEightButton.clicked.connect(lambda: self.buttonClicked(self.zernikeEightButton))

        
        closeButton = QPushButton("Close")
        closeButton.clicked.connect(self.closeButtonClicked)
        resetButton = QPushButton("Reset")
        resetButton.clicked.connect(self.resetButtonClicked)

        layout = QGridLayout()
        layout.addWidget(self.zernikeFourButton,0,0)
        layout.addWidget(self.zernikeSixButton,1,0)
        layout.addWidget(self.zernikeEightButton,2,0)
        layout.addWidget(resetButton,4,0)
        layout.addWidget(closeButton,4,1)
        self.setLayout(layout)
        self.setWindowTitle("Zernike Order")

        #self.setRadioButtons()

    #     The action buttons


    def buttonClicked(self,button):
        global ZernikeOrder
        ZernikeOrder = button.order
        print("Order " + str(ZernikeOrder))
        
    def resetButtonClicked(self):
        global ZernikeOrder
        ZernikeOrder = 4
        self.zernikeFourButton.setChecked(True)

        
    def closeButtonClicked(self):      # Close the frame
        self.close()
        if self.closeAction != None:
            self.closeAction()


class ReferenceOptionSetter(QWidget):
    """
    Widget to set the reference option with Radio bottons.
    """
    def __init__(self, parent = None,closeAction = None):
        super(ReferenceOptionSetter,self).__init__(parent)

        refopt = getReferencePointOption()
        self.closeAction = closeAction
        self.setAutoFillBackground(True)
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.white)
        self.setPalette(p)

        #     Setup the the three option buttons
        self.paraxialButton = QRadioButton("Paraxial")
        if refopt == 0:
            self.paraxialButton.setChecked(True)
        self.paraxialButton.option = 0
        self.paraxialButton.clicked.connect(lambda: self.buttonClicked(self.paraxialButton))
        self.inplaneButton = QRadioButton("In plane")
        if refopt == 1:
            self.inplaneButton.setChecked(True)
        self.inplaneButton.option = 1
        self.inplaneButton.clicked.connect(lambda: self.buttonClicked(self.inplaneButton))
        self.optimalButton = QRadioButton("Optimal")
        if refopt == 2:
            self.optimalButton.setChecked(True)
        self.optimalButton.option = 2
        self.optimalButton.clicked.connect(lambda: self.buttonClicked(self.optimalButton))

        # Setup the stardard reset and close buttons.
        closeButton = QPushButton("Close")
        closeButton.clicked.connect(self.closeButtonClicked)
        resetButton = QPushButton("Reset")
        resetButton.clicked.connect(self.resetButtonClicked)

        #   Set grid layout
        layout = QGridLayout()
        layout.addWidget(self.paraxialButton,0,0)
        layout.addWidget(self.inplaneButton,1,0)
        layout.addWidget(self.optimalButton,2,0)
        layout.addWidget(resetButton,4,0)
        layout.addWidget(closeButton,4,1)
        self.setLayout(layout)
        self.setWindowTitle("Reference Option")

    #     The action buttons


    def buttonClicked(self,button):
        setReferencePointOption(button.option)
    
        
    def resetButtonClicked(self):
        setReferencePointOption(1)
        self.inplaneButton.setChecked(True)

        
    def closeButtonClicked(self):      # Close the frame
        self.close()
        if self.closeAction != None:
            self.closeAction()

class TiltSetter(QWidget):
    """
    Set the interferometer tilts
    """
    def __init__(self, parent = None,closeAction = None):
        super(TiltSetter,self).__init__(parent)

        global Xtilt, Ytilt

        self.closeAction = closeAction
        self.setAutoFillBackground(True)
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.white)
        self.setPalette(p)

        xtiltLabel = QLabel("X tilt")
        self.xtiltSetter = QDoubleSpinBox()
        self.xtiltSetter.setValue(Xtilt)
        self.xtiltSetter.setSingleStep(0.1)
        self.xtiltSetter.valueChanged.connect(self.xtiltSetterClicked)

        ytiltLabel = QLabel("Y tilt")
        self.ytiltSetter = QDoubleSpinBox()
        self.ytiltSetter.setValue(Ytilt)
        self.ytiltSetter.setSingleStep(0.1)
        self.ytiltSetter.valueChanged.connect(self.ytiltSetterClicked)
        
        closeButton = QPushButton("Close")
        closeButton.clicked.connect(self.closeButtonClicked)
        resetButton = QPushButton("Reset")
        resetButton.clicked.connect(self.resetButtonClicked)
        

        layout = QGridLayout()     # Vertical box
        layout.addWidget(xtiltLabel,0,0)
        layout.addWidget(self.xtiltSetter,0,1)
        layout.addWidget(ytiltLabel,1,0)
        layout.addWidget(self.ytiltSetter,1,1)
        layout.addWidget(resetButton,2,0)
        layout.addWidget(closeButton,2,1)
        self.setLayout(layout)



    def xtiltSetterClicked(self):
        global Xtilt
        x = self.xtiltSetter.value()
        Xtilt  = float(x)

    def ytiltSetterClicked(self):
        global Ytilt
        y = self.ytiltSetter.value()
        Ytilt  = float(y)
        

    def resetButtonClicked(self):
        global Xtilt
        global Ytilt
        Xtilt = 3.0
        self.xtiltSetter.setValue(Xtilt)
        Ytilt = 0.0
        self.ytiltSetter.setValue(Ytilt)
        
    def closeButtonClicked(self):      # Close the frame
        self.close()
        if self.closeAction != None:
            self.closeAction()
        


class KnifeSetter(QWidget):
    """
    Widget to set the knife paramteers
    """
    def __init__(self, parent = None,closeAction = None):
        super(KnifeSetter,self).__init__(parent)

        # global CurrentKnife, CurrentKnifeAngle,CurrentKnifeShift

        self.closeAction = closeAction
        self.setAutoFillBackground(True)
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.white)
        self.setPalette(p)

        knifeLabel = QLabel("Knife Position")
        self.knifeSetter = QDoubleSpinBox()
        self.knifeSetter.setValue(CurrentKnife)
        self.knifeSetter.setSingleStep(0.01)
        self.knifeSetter.valueChanged.connect(self.knifeSetterClicked)
        angleLabel = QLabel("Knife Angle")
        self.knifeAngleSetter = QDoubleSpinBox()
        self.knifeAngleSetter.setValue(math.degrees(CurrentKnifeAngle))
        self.knifeAngleSetter.setSingleStep(1.0)
        self.knifeAngleSetter.setRange(-90.0,90.0)
        self.knifeAngleSetter.valueChanged.connect(self.knifeAngleSetterClicked)
        shiftLabel = QLabel("Axial Shift")
        self.shiftSetter = QDoubleSpinBox()
        self.shiftSetter.setSingleStep(0.1)
        self.shiftSetter.setRange(-20.0,20.0)
        self.shiftSetter.setValue(CurrentKnifeShift)
        self.shiftSetter.valueChanged.connect(self.shiftSetterClicked)
        self.wireSetter = QRadioButton("Wire")
        self.wireSetter.clicked.connect(self.wireSetterClicked)
        self.wireSetter.setChecked(CurrentWire)
                
        closeButton = QPushButton("Close")
        closeButton.clicked.connect(self.closeButtonClicked)
        resetButton = QPushButton("Reset")
        resetButton.clicked.connect(self.resetButtonClicked)

        layout = QGridLayout()     # Vertical box
        layout.addWidget(knifeLabel,0,0)
        layout.addWidget(self.knifeSetter,0,1)
        layout.addWidget(angleLabel,1,0)
        layout.addWidget(self.knifeAngleSetter,1,1)
        layout.addWidget(shiftLabel,2,0)
        layout.addWidget(self.shiftSetter,2,1)
        layout.addWidget(self.wireSetter,3,0)
        layout.addWidget(resetButton,4,0)
        layout.addWidget(closeButton,4,1)
        self.setLayout(layout)

    #         Fun to set the buttons
    #
    def knifeSetterClicked(self):
        global CurrentKnife
        v = self.knifeSetter.value()
        CurrentKnife = float(v)


    def knifeAngleSetterClicked(self):
        global CurrentKnifeAngle
        v = self.knifeAngleSetter.value()
        CurrentKnifeAngle = math.radians(v)

    def shiftSetterClicked(self):
        global CurrentKnifeShift
        v = self.shiftSetter.value()
        CurrentKnifeShift = float(v)

    def wireSetterClicked(self):
        global CurrentWire
        CurrentWire = not CurrentWire
        self.wireSetter.setChecked(CurrentWire)
        
        

    def resetButtonClicked(self):
        global CurrentKnife
        global CurrentKnifeAngle
        global CurrenntKnifeShift
        CurrentKnife = 0.0
        self.knifeSetter.setValue(CurrentKnife)
        CurrentKnifeAngle = 0.0
        self.knifeAngleSetter.setValue(CurrentKnifeAngle)
        CurrentKnifeShift = 0.0
        self.shiftSetter.setValue(CurrentKnifeShift)
        
        
    def closeButtonClicked(self):      # Close the frame
        self.close()
        if self.closeAction != None:
            self.closeAction()
        

class MessageBox(QMessageBox):
    """
    Class to display messages
    """
    def __init__(self, itext = "None",dtext = "None", parent = None):
        super(MessageBox,self).__init__(parent)

        self.setIcon(QMessageBox.Information)
        self.setText(getCurrentLens().title)
        self.setInformativeText(itext)
        self.setWindowTitle("Message")
        self.setDetailedText(dtext)
        self.setStandardButtons(QMessageBox.Ok )

class PltMainWindow(QMainWindow):
    """ Class to set up a main window with the central pane being a Matplotlib pane
    
    :param lens: The lens to be shown or manipulated. (Default = None)
    :type lens: OpticalGroiup or extending class
    :param parent: the parent window if there is ine (Default = None)
    """
    def __init__(self, lens = None, parent=None):
        super(PltMainWindow, self).__init__(parent)

        # Set size of panel
        self.resize(PanelSize[0],PanelSize[1])
        #  If lens given then set current lens, if not use default lens.
        if lens != None:
            setCurrentLens(lens) 

        #     Setup components for plt window
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        toolbar = NavigationToolbar(self.canvas, self)
        self.menubar = self.menuBar()
        self.menubar.setNativeMenuBar(False)
        
        #      Set up basic menu items 
        fileMenu = self.menubar.addMenu("File")
        fileAction = QAction("New Lens",self)
        fileAction.triggered.connect(self.fileButtonClicked)
        fileMenu.addAction(fileAction)
        plotAction = QAction("Replot",self)
        plotAction.triggered.connect(self.plot)
        fileMenu.addAction(plotAction)
        lensAction = QAction("Show Lens",self)
        lensAction.triggered.connect(self.lensPlot)
        fileMenu.addAction(lensAction)
        infoAction = QAction("Lens Info",self)
        infoAction.triggered.connect(self.infoButtonClicked)
        fileMenu.addAction(infoAction)
        exitAction = QAction("Exit",self)
        exitAction.triggered.connect(sys.exit)
        fileMenu.addAction(exitAction)

        #      Setup central panel
        
        layout = QVBoxLayout()
        layout.addWidget(toolbar)
        layout.addWidget(self.canvas)
        panel = QWidget()
        panel.setLayout(layout)
        self.setCentralWidget(panel)

        #      Setup basic options menu (common to all intrefaces)
        optionMenu = self.menubar.addMenu("Options")
        waveAction = QAction("Wavelength",self)
        waveAction.triggered.connect(self.waveButtonClicked)
        optionMenu.addAction(waveAction)
        angleAction = QAction("Angle",self)
        optionMenu.addAction(angleAction)
        angleAction.triggered.connect(self.angleButtonClicked)
        irisAction = QAction("Iris",self)
        optionMenu.addAction(irisAction)
        irisAction.triggered.connect(self.irisButtonClicked)
        referenceAction = QAction("Reference Point",self)
        referenceAction.triggered.connect(self.referenceButtonClicked)
        optionMenu.addAction(referenceAction)
        
        # Do a default plot of the lens if one exits
        if getCurrentLens() == None:
            self.fileButtonClicked()
        else:
            self.lensPlot()

        
    #        The basic buttons


    def fileButtonClicked(self):
        """
        New File button. 
        """
        fs = LensSetter(parent=self,closeAction=self.lensPlot)

    def infoButtonClicked(self):
        """
        Lens Information button.
        """
        m = MessageBox("Lens Information",getCurrentLens().getInfo(),parent = self)
        m.setWindowTitle("Information")
        m.show()
    
    def waveButtonClicked(self):
        """
        Wavelength setter clicked
        """
        ws = WaveLengthSetter(parent=self,closeAction=self.plot)
        ws.move(50,50)
        ws.resize(200,100)
        ws.show()

    def angleButtonClicked(self):
        """
        Ray angle button
        """
        aset = DirectionSetter(parent=self,closeAction=self.plot)
        aset.move(50,50)
        aset.resize(300,150)
        aset.show()

    def irisButtonClicked(self):
        """
        Iris setter button
        """

        ir = IrisSetter(parent=self,closeAction=self.plot)
        ir.move(50,50)
        ir.resize(200,200)
        ir.show()

    def referenceButtonClicked(self):
        """
        Reference setter
        """
        w = ReferenceOptionSetter(parent=self,closeAction=self.plot)
        w.move(50,50)
        w.resize(200,200)
        w.show()

    #     The plot method

    def plot(self):
        """    
        Method to plot the actual image
        """
        self.figure.clear()      # Clear current figure
        #       Use a subPlot() method to do the actual work 
        self.subPlot()
        #refresh the canvas
        self.canvas.draw()

   
    def lensPlot(self):
        """
        The default lens plot
        """
        self.figure.clear()      # Clear current figure
        panel = self.figure.add_subplot(111)
        panel.axis('equal')

        # plot data
        getCurrentLens().draw(planes = False)
        plt.grid()
        plt.xlabel("Optical Axis")
        plt.ylabel("Height")
        plt.title("Diagram of lens " + getCurrentLens().title)
        self.canvas.draw()
        
class LensViewer(PltMainWindow):
    """
    Class to plot a lens the panel with a collimated beam of rays.
    
    :param lens: the lens to plot (Default = None)
    :type lens: OpticalGroup or extending class.
    :param parent: Parent window if there is one (Default = None)
    
    """
    def __init__(self, lens = None , parent=None):
        super(LensViewer, self).__init__(lens,parent)

        if getCurrentLens() != None:
            self.plot() 
        
        
    def subPlot(self):
        """
        Method to do the actual plot, automatically called
        """
        panel = self.figure.add_subplot(111)
        panel.axis('equal')

        u = getCurrentAngle()
        pencil = RayPencil().addBeam(getCurrentLens(),u,"vl",\
                                               wave=getCurrentWavelength()).addMonitor(RayPath())

        #
        #        Set the output plane (being the back focal plane)
        op = getCurrentLens().backFocalPlane()

        #         Propagate pencil through lens and one to back plane
        pencil *= getCurrentLens()       # Through lens
        pencil *= op         # To plane

        
        # plot data
        getCurrentLens().draw(True,True)
        op.draw()
        pencil.draw()
        plt.grid()
        plt.xlabel("Optical Axis")
        plt.ylabel("Height")
        plt.title("Diagram of lens " + getCurrentLens().title)



class AbberationViewer(PltMainWindow):
    """
    Class to plot the three aberration plots for a lens at choice of angles
    and wavelength
    
    :param lens: the lens to plot (Default = None)
    :type lens: OpticalGroup or extending class.
    :param parent: Parent window if there is one (Default = None)    
    
    """
    def __init__(self, lens = None , parent=None):
        super(AbberationViewer, self).__init__(lens,parent)

        if getCurrentLens() != None:
            self.plot()

    def subPlot(self):
        """
        Method to do the actual plot (called automatically
        """
        wa = WaveFrontAnalysis(getCurrentLens(),getDesignWavelength())
        wa.drawAberrationPlot(getCurrentAngle(),getCurrentWavelength())

        
        
        
class WaveFrontViewer(PltMainWindow):
    """
    Class to plot a lens the panel with rays
    """
    def __init__(self, lens = None , parent=None):
        super(WaveFrontViewer, self).__init__(lens,parent)

        waveMenu = self.menubar.addMenu("Wave")
        orderAction = QAction("Zerkike Order",self)
        orderAction.triggered.connect(self.orderButtonClicked)
        waveMenu.addAction(orderAction)
        zernikeInfoAction = QAction("Zernike Details",self)
        zernikeInfoAction.triggered.connect(self.zernikeButtonClicked)
        waveMenu.addAction(zernikeInfoAction)
        tiltAction = QAction("Tilts",self)
        tiltAction.triggered.connect(self.tiltButtonClicked)
        waveMenu.addAction(tiltAction)
        print("Wave menu added")
        self.interferometer = Interferometer()


    def subPlot(self):
        global CurrentWaveFront
        wa = WaveFrontAnalysis(getCurrentLens(),getDesignWavelength())
        CurrentWaveFront = wa.fitZernike(getCurrentAngle(),getDefaultWavelength(),ZernikeOrder,ReferencePointOption)
        self.interferometer.setWaveFront(CurrentWaveFront)
        self.displayPlot()


    def displayPlot(self):
        self.interferometer.setTilt(Xtilt,Ytilt)
        self.interferometer.draw()
        #CurrentWaveFront.plotImage(xtilt=Xtilt,ytilt=Ytilt)
        print("Xtilt is " + str(Xtilt) + " and Ytilt is " + str(Ytilt))


    #     Wave button aclions

    def orderButtonClicked(self):
        """
        Wavelength setter
        """
        zs = ZernikeOrderSetter(parent=self,closeAction=self.plot)
        zs.move(50,50)
        zs.resize(200,100)
        zs.show()

    def tiltButtonClicked(self):
        """
        The tilt button
        """
        tb = TiltSetter(parent=self,closeAction=self.plot)
        tb.move(50,50)
        tb.resize(200,100)
        tb.show()
        
        
    def zernikeButtonClicked(self):
        m = MessageBox("Zernike Expansion w: {0:4.2f}".format(getDefaultWavelength()),repr(CurrentWaveFront),parent = self)
        m.setWindowTitle("Information")
        m.show()

class OTFViewer(WaveFrontViewer):
    """
    Class extending WaveFrontViewer to display OFT 
    """
    def __init__(self,lens = None, parent = None):
        super(OTFViewer,self).__init__(lens,parent)

    def displayPlot(self):
        CurrentWaveFront.plotOTF(128,"b")

        
class KnifeViewer(PltMainWindow):
    """
    Class to plot a lens the panel with rays
    """
    def __init__(self, lens = None , parent=None):
        super(KnifeViewer, self).__init__(lens,parent)

        knifeMenu = self.menubar.addMenu("Knife")
        knifeAction = QAction("Knife Parameters",self)
        knifeAction.triggered.connect(self.knifeButtonClicked)
        knifeMenu.addAction(knifeAction)


    def subPlot(self):
        kt = KnifeTest(getCurrentLens(),getCurrentAngle(),getCurrentWavelength(),getDesignWavelength())     # New knife test
        kt.setKnife(CurrentKnife,CurrentKnifeAngle,CurrentKnifeShift)   # Set knife
        kt.setWire(CurrentWire)
        kt.setReference(getReferencePointOption())
        kt.getImage().draw()                        # make and plot image
        plt.title(getCurrentLens().title) 
        
    #   Additional buttons in menue

    def knifeButtonClicked(self):
        w = KnifeSetter(parent = self, closeAction = self.plot)
        w.move(50,50)
        w.resize(200,200)
        w.show()


class SpotViewer(PltMainWindow):
    """
    Class to plot spot diagrams of a lens with a collimnated beam with options 
    to change angle, wavelength and plane of the spots.
    
    :param lens: the lens to plot (Default = None)
    :type lens: OpticalGroup or extending class.
    :param parent: Parent window if there is one (Default = None)
    """
    def __init__(self, lens = None , parent=None):
        super(SpotViewer, self).__init__(lens,parent)


        self.delta = 0.1

        #  Make the additioal menu to control the position of the stop plane
        spotMenu = self.menubar.addMenu("Plane")
        plusAction = QAction("Plus",self)
        minusAction = QAction("Minus",self)
        controlAction = QAction("Variable",self)
        plusAction.triggered.connect(self.plusClicked)
        minusAction.triggered.connect(self.minusClicked)
        controlAction.triggered.connect(self.variableClicked)

        spotMenu.addAction(plusAction)
        spotMenu.addAction(minusAction)
        spotMenu.addAction(controlAction)

    def plusClicked(self):
        """
        Move plane forward along optical axis
        """
        incrementPlaneShift(self.delta)
        self.updatePlane()

    def minusClicked(self):
        """
        Move the plane back along the optical axis
        """        
        incrementPlaneShift(-self.data)
        self.updatePlane()


    def variableClicked(self):
        """
        Update the plane osition with the dial
        """
        p = PlaneSetter(scale = 0.01, parent = self , closeAction = None , changedAction = self.updatePlane)
        p.move(50,50)
        p.resize(200,200)
        p.show()
        
    def subPlot(self):
        """   
        Method to set up the spot analysis, and trace the rays. This need to be called
        if any of the geometry of the sytsem changes.
        """
        pencil = RayPencil().addBeam(getCurrentLens(),getCurrentAngle(),"array",wave=getCurrentWavelength())
        bf = getCurrentLens().backFocalPlane(getDesignWavelength())

        pencil *= getCurrentLens()
        pencil *= bf
        
        refopt = getReferencePointOption()

        if refopt == 0:
            self.pt = getCurrentLens().imagePoint(getCurrentAngle(),getDesignWavelength())    # Paraxial point location
        elif refopt == 1:
            self.pt = Psf().setWithRays(pencil,bf)              # Centre of PSF in image plane
        elif refopt == 2:
            self.pt = Psf().optimalArea(pencil,bf)              # Optimal area PSF, not in image plane
        else:
            print("Illegal ref type")


        self.spot = SpotDiagram(pencil)
        self.updatePlane()

    def updatePlane(self):
        """
        Method to update the position of the spot plane only.

        """
        self.figure.clear()
        pos = self.pt.z + getPlaneShift()
        plane = OpticalPlane(pos)
        self.spot.draw(plane)
        plt.title(getCurrentLens().title)
        plt.xlabel("Shift : {0:5.4f} Plane at : {1:7.4f}".format(getPlaneShift(),pos))
        self.canvas.draw()



        

