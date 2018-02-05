"""
   Basic use of the file handing methods

Author: Will Hossack, The University of Edinburgh.
"""
import tio as t
import math

def main():

    n = t.getFilename("File Name","csv","$HOME/mydata")
    t.tprint("Name is : " ,n)

main()
