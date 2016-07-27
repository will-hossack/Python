#
#           Poitential contour example
#
import tio as t
import particle as p
import vector as v
import matplotlib.pyplot as plt
import numpy as np

def main():
    ps = p.ParticleSystem().readFile(t.openFile("File",defaulttype="part"))
    t.tprint(ps)
    
    delta = 0.025
    x = np.arange(-3.0, 3.0, delta)
    y = np.arange(-2.0, 2.0, delta)
    X, Y = np.meshgrid(x, y)

    
