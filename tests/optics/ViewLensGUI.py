from PyQt5.QtWidgets import *
from PyQt5.QtGui import QIcon

import optics.lens as l
from optics.gui import *


from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt

#import optics.lens as len
import optics.ray as ray
import vector as v
import optics.wavelength as w 



class LensPlot(QWidget):
    def __init__(self, lens, parent=None):
        super(LensPlot, self).__init__(parent)

        self.lens = lens
        # a figure instance to plot on
        self.figure = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        # Just some button connected to `plot` method
        self.button = QPushButton('Plot')
        self.button.clicked.connect(self.plot)

        # set the layout
        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        layout.addWidget(self.button)
        self.setLayout(layout)
        #self.plot()

    def plot(self):

        # instead of ax.hold(False)
        self.figure.clear()

        # create an axis
        panel = self.figure.add_subplot(111)
        panel.axis('equal')
        # discards the old graph
        # ax.hold(False) # deprecated, see above


        angle = 2.0
        u = v.Unit3d(v.Angle().setDegrees(angle))       # Unit3d of ray pencil.
        print("wave is : " + str(w.Default))
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



    

def main():
    app = QApplication([])
    #al = DirectionSetter()
    #al.show()
    #wl = WaveLengthSetter()
    #wl.show()
    lens = l.DataBaseLens("$LENS/Tessar-F2.8")
    ex = LensViewer(lens)
    ex.show()
    app.exec_()

main()
