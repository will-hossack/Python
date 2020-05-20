
"""
    Test of the default /degign wavelength
"""

from optics.wavelength import Default,getDefaultWavelength,setDefaultWavelength,MaterialIndex
import tio as t




#from optics.ray import RayPencil,Disc
import os
from importlib.resources import open_text,path
from optics.lens import DataBaseLens



def getW(w = None):
    if w == None:
        w = getDefaultWavelength()
        
    return w

def main():
    
    """
    while True:
        try:
            fileName = t.getFilename("Lens","lens")
            print(fileName)
            file = open_text("optics.lenses",fileName)
            break
        except:
            t.tprint("No suc lens" + fileName)
            
        
    lines = file.readlines()
    file.close()
    for l in lines:
        print(l,end="")
    """
    
    with path("optics","lenses") as p:
        lenspath = p
        
    print(str(lenspath))
    lensName = t.getFilename("Lens","lens")
    print(lensName) 
    
    filename = os.path.join(lenspath,lensName)
    
    print(filename)
    
    lens = DataBaseLens(filename)
    print(lens.getInfo())
    
    
    
    #p = resources.path("optics.lenses","materials.data")
    #print(str(p.exists()))
    
    #cont = resources.contents("optics.lenses")
    #for s in cont:
    #    print(s)
    
    #       n = MaterialIndex("BK7")
    #    print(repr(n))
    #   n = MaterialIndex("F2")
    #   print(repr(n))
    
    #directory,filename = os.path.split(__file__)
    #print("Director is : " + directory)
    
    
    
    #   print("Ray is " + repr(pencil))
    
main()