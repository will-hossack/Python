"""
     Set of gui classes and methods
"""
from os import getenv
import sys
import math
from optics.wavelength import Green,getDefaultWavelength,setDefaultWavelength,getDesignWavelength,setDesignWavelength
from optics.lens import setCurrentLens,getCurrentLens
from optics.wavefront import WaveFrontAnalysis
from optics.analysis import KnifeTest
import optics.ray as ray
import optics.psf as psf
from optics.surface import OpticalPlane
from vector import Unit3d
from PyQt5.QtWidgets import *
from PyQt5.QtCore import Qt

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt


CurrentAngle = Unit3d(0.0,0.0,1.0)
IrisRatio = 1.0
ZernikeOrder = 4
ReferencePointOption = 1
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
    Class to set the Default and Design wavelengths with spinners
    """
    def __init__(self, parent = None,closeAction = None):
        super(WaveLengthSetter,self).__init__(parent)

        self.closeAction = closeAction
        self.setAutoFillBackground(True)
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.white)
        self.setPalette(p)

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
        resetButton = QPushButton("Reset")
        resetButton.clicked.connect(self.resetButtonClicked)


        
        layout.addWidget(defaultLabel,0,0)      # Add the 5 item in Grid
        layout.addWidget(self.defaultSpin,0,1)
        layout.addWidget(designLabel,1,0)
        layout.addWidget(self.designSpin,1,1)
        layout.addWidget(resetButton,2,0)
        layout.addWidget(closeButton,2,1)

        self.setLayout(layout)
        self.setWindowTitle("WaveLengthSetter")

    """
    Method to update the values from the spinners and close the widow
    """
        
    def defaultValueChange(self):
        v = self.defaultSpin.value()
        setDefaultWavelength(v)

    def designValueChange(self):
        v = self.designSpin.value()
        setDesignWavelength(v)

    def resetButtonClicked(self):
        self.defaultSpin.setValue(Green)
        self.designSpin.setValue(Green)
        setDefaultWavelength(Green)
        setDesignWavelength(Green)
        
    def closeButtonClicked(self):      # Close the frame
        self.close()
        if self.closeAction != None:
            self.closeAction()
        



class DirectionSetter(QWidget):
    """
    Class to set default direction in degress with spinners
    """
    def __init__(self, parent = None,closeAction = None):
        super(DirectionSetter,self).__init__(parent)

        self.closeAction = closeAction
        self.setAutoFillBackground(True)
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.white)
        self.setPalette(p)
        
        thetaLabel = QLabel("Theta in degress :")     # The labels
        psiLabel = QLabel("Psi in degrees : ")
        self.unitLabel = QLabel(str(CurrentAngle))

        self.theta,self.psi = CurrentAngle.getAngle().getDegrees()  # Current Theta/Psi
        self.thetaSpin = QDoubleSpinBox()       # Theta wave spinner
        self.thetaSpin.setRange(-90.0,90.0)
        self.thetaSpin.setValue(self.theta)
        self.thetaSpin.setSingleStep(1.0)
        self.thetaSpin.valueChanged.connect(self.valueChange)

        self.psiSpin = QDoubleSpinBox()       # Theta wave spinner
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


class LensSetter(QWidget):
    """
    Class to set default direction in degress with spinners
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
    Widget to set the reference option
    """
    def __init__(self, parent = None,closeAction = None):
        super(ReferenceOptionSetter,self).__init__(parent)

        global ReferencePointOption

        self.closeAction = closeAction
        self.setAutoFillBackground(True)
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.white)
        self.setPalette(p)

        self.paraxialButton = QRadioButton("Paraxial")
        if ReferencePointOption == 0:
            self.paraxialButton.setChecked(True)
        self.paraxialButton.option = 0
        self.paraxialButton.clicked.connect(lambda: self.buttonClicked(self.paraxialButton))
        self.inplaneButton = QRadioButton("In plane")
        if ReferencePointOption == 1:
            self.inplaneButton.setChecked(True)
        self.inplaneButton.option = 1
        self.inplaneButton.clicked.connect(lambda: self.buttonClicked(self.inplaneButton))
        self.optimalButton = QRadioButton("Optimal")
        if ReferencePointOption == 2:
            self.optimalButton.setChecked(True)
        self.optimalButton.option = 2
        self.optimalButton.clicked.connect(lambda: self.buttonClicked(self.optimalButton))

        
        closeButton = QPushButton("Close")
        closeButton.clicked.connect(self.closeButtonClicked)
        resetButton = QPushButton("Reset")
        resetButton.clicked.connect(self.resetButtonClicked)

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
        global ReferencePointOption
        ReferencePointOption = button.option
    
        
    def resetButtonClicked(self):
        global ReferncePointOption
        ZernikeOrder = 1
        self.inplaneButton.setChecked(True)

        
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
    """
    def __init__(self, lens = None, parent=None):
        super(PltMainWindow, self).__init__(parent)

        if lens != None:
            setCurrentLens(lens) 

        #     Setup components for plt window
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        toolbar = NavigationToolbar(self.canvas, self)
        self.menubar = self.menuBar()
        self.menubar.setNativeMenuBar(False)
        
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
        
    
        if getCurrentLens() == None:
            self.fileButtonClicked()
        else:
            self.lensPlot()

        
    #        The basic buttons


    def fileButtonClicked(self):
        fs = LensSetter(parent=self,closeAction=self.lensPlot)

    def infoButtonClicked(self):
        m = MessageBox("Lens Information",getCurrentLens().getInfo(),parent = self)
        m.setWindowTitle("Information")
        m.show()
    
    def waveButtonClicked(self):
        """
        Wavelength setter
        """
        ws = WaveLengthSetter(parent=self,closeAction=self.plot)
        ws.move(50,50)
        ws.resize(200,100)
        ws.show()

    def angleButtonClicked(self):
        aset = DirectionSetter(parent=self,closeAction=self.plot)
        aset.move(50,50)
        aset.resize(300,150)
        aset.show()

    def irisButtonClicked(self):
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
        """    Method to plot the actual image
        """
        self.figure.clear()      # Clear current figure
        #       Use a subPlot() method to do the actual work 
        self.subPlot()
        #refresh the canvas
        self.canvas.draw()

    """       Default lens plot always available
    """
    def lensPlot(self):
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
    Class to plot a lens the panel with rays
    """
    def __init__(self, lens = None , parent=None):
        super(LensViewer, self).__init__(lens,parent)

        if CurrentLens != None:
            self.plot() 
        
    def subPlot(self):
        panel = self.figure.add_subplot(111)
        panel.axis('equal')

        u = CurrentAngle
        pencil = ray.RayPencil().addCollimatedBeam(getCurrentLens(),u,"vl",\
                                                   wave=getDefaultWavelength()).addMonitor(ray.RayPath())

        #
        #        Set the output plane (being the back focal plane)
        op = getCurrentLens().backFocalPlane()

        #         Propagate pencil through lens and one to back plane
        pencil *= CurrentLens       # Through lens
        pencil *= op         # To plane

        
        # plot data
        getCurrentLens().draw()
        op.draw()
        pencil.draw()
        plt.grid()
        plt.xlabel("Optical Axis")
        plt.ylabel("Height")
        plt.title("Diagram of lens " + getCurrentLens().title)



class AbberationViewer(PltMainWindow):
    """
    Class to plot a lens the panel with rays
    """
    def __init__(self, lens = None , parent=None):
        super(AbberationViewer, self).__init__(lens,parent)

        if getCurrentLens() != None:
            self.plot()

    def subPlot(self):
        wa = WaveFrontAnalysis(getCurrentLens(),getDesignWavelength())
        wa.drawAberrationPlot(CurrentAngle,getDefaultWavelength())

        
        
        
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


    def subPlot(self):
        global CurrentWaveFront
        wa = WaveFrontAnalysis(getCurrentLens(),getDesignWavelength())
        CurrentWaveFront = wa.fitZernike(CurrentAngle,getDefaultWavelength(),ZernikeOrder,ReferencePointOption)
        self.fringePlot()


    def fringePlot(self):
        CurrentWaveFront.plotImage(xtilt=Xtilt,ytilt=Ytilt)


    #     Wave button aclions

    def orderButtonClicked(self):
        """
        Wavelength setter
        """
        zs = ZernikeOrderSetter(parent=self,closeAction=self.plot)
        zs.move(50,50)
        zs.resize(200,100)
        zs.show()
        
        
    def zernikeButtonClicked(self):
        m = MessageBox("Zernike Expansion w: {0:4.2f}".format(getDefaultWavelebgth()),repr(CurrentWaveFront),parent = self)
        m.setWindowTitle("Information")
        m.show()



        
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
        kt = KnifeTest(getCurrentLens(),CurrentAngle,getDefaultWavelength(),getDesignWavelength())     # New knife test
        kt.setKnife(CurrentKnife,CurrentKnifeAngle,CurrentKnifeShift)   # Set knife
        kt.setWire(CurrentWire)
        kt.getImage(ReferencePointOption).draw()                        # make and plot image
        plt.title(getCurrentLens().title) 
        
    #   Additional buttons in menue

    def knifeButtonClicked(self):
        w = KnifeSetter(parent = self, closeAction = self.plot)
        w.move(50,50)
        w.resize(200,200)
        w.show()


class SpotViewer(PltMainWindow):
    """
    Class to plot a lens the panel with rays
    """
    def __init__(self, lens = None , parent=None):
        super(SpotViewer, self).__init__(lens,parent)

        spotMenu = self.menubar.addMenu("Plane")

        plusAction = QAction("Plus",self)
        minusAction = QAction("Minus",self)
        plusAction.triggered.connect(self.plusClicked)
        minusAction.triggered.connect(self.minusClicked)

        self.delta = 0.0

        spotMenu.addAction(plusAction)
        spotMenu.addAction(minusAction)

    def plusClicked(self):
        print("Plus")
        self.delta += 0.2
        self.updatePlane()

    def minusClicked(self):
        self.delta -= 0.2
        self.updatePlane()

    def subPlot(self):
        """        Do a full plot 
        """
        pencil = ray.RayPencil().addCollimatedBeam(getCurrentLens(),CurrentAngle,"array",wave=getDefaultWavelength())
        bf = getCurrentLens().backFocalPlane(getDesignWavelength())

        pencil *= getCurrentLens()
        pencil *= bf

        if ReferencePointOption == 0:
            self.pt = getCurrentLens().imagePoint(CurrentAngle,getDesignWavelength())    # Paraxial point location
        elif ReferencePointOption == 1:
            self.pt = psf.Psf().setWithRays(pencil,bf)              # Centre of PSF in image plane
        elif ReferencePointOption == 2:
            self.pt = psf.Psf().optimalArea(pencil,bf)              # Optimal area PSF, not in image plane
        else:
            print("Illegal ref type")


        self.spot = psf.SpotDiagram(pencil)
        self.updatePlane()

    def updatePlane(self):
        self.figure.clear()
        plane = OpticalPlane(self.pt.z + self.delta)
        self.spot.draw(plane)
        self.canvas.draw()



        

