"""       A Vector3d class using np.array for internal storage
"""

import math
import numpy as np
import random

class Vector3d(object):
    """
    Class to impolement a Vector3d
    """
    def __init__(self,x = 0.0, y = 0.0, z = 0.0):
        self.data = np.array([x,y,z])

def main():
    a = Vector3d()
    print(a.data)

main()

