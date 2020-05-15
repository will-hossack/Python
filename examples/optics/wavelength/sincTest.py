#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 16:39:47 2020

@author: wjh
"""

import numpy as np
from scipy.special import sinc,j1
import matplotlib.pyplot as plt
import math

def jinc(x):
    if x == 0.0:
        return 1.0
    else:
        return 2*j1(x)/x
    
    
    
def main():
    
    x = np.linspace(-10,10,300)
    
    plt.plot(x,sinc(x/math.pi))
    plt.plot(x,(2*j1(x)/x))
    plt.show()
    
    
main()