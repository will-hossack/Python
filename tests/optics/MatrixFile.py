"""
   Test of Matrix File class
"""
import optics.matrix as m
import vector as v
import tio as t
import  matplotlib.pyplot as plt

def main():
    
    pm =  m.ParaxialGroup().readFile()
    t.tprint(pm) 
    t.tprint("Focal length is ",pm.backFocalLength())

main()
