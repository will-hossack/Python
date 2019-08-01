import tio as t
import optics.analysis as anal
import matplotlib.pyplot as plt
import vector as v


def main():

    target = anal.TargetPlane([20.0,30.0,100.0],100,100)

    print(repr(target))

    t = v.Vector2d()
    target.addGrid(11,5)

    fig = plt.figure()
    panel = fig.add_subplot(1,1,1)
    panel.axis('equal')

    target.draw()
    plt.show()

main()
