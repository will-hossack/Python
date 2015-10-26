"""
Classes to handle the refractive index properties of materails, mainly glasses, but
also gases and liquids where are held in the asci databases using the same 
conventions as in RefractiveIndex.org.

Author:   Will Hossack, The University of Edinburgh
""" 
from tio import getExpandedFilename,getOption
import wavelength as wl

#                    
#

DataBaseFile = "$LENS/materials.data" 
DataBase = None

class MaterialData(object):
    """
    Class to hold the Material Data file
    """
    #
    #
    def __init__(self,filename = DataBaseFile):
        """
        Constructor to open the materials definition file and read in the 
        data to self.data as a list of lines. 
        The database is only loaded once.
        param filename the database (defaults to package default)
        """
        global DataBase
        if DataBase == None :
            try :
                filename = getExpandedFilename(filename)   # Sort out logicals
                filestream = open(filename,"r")
                DataBase = filestream.readlines()
                filestream.close()
            except :
                print("materail.MaterialData() unable to open data file {0:s}".format(filename))
                raise IOError("optics: Unable to run without materail data")
    #
    #      
    def getList(self):
        """
        Method to a list of the materials keys in the database as a list of strings.
        """
        key = []
        for line in DataBase:
            if not line.startswith("#") or len(line.strip()) == 0:
                token = line.split()
                key.append(token[0].strip())
        return key
        
    #    
    #
    def getMaterial(self,key = None):
        """
        Method to get a material from the loaded database.
        parm key name, (usually glass name) if None, then it will be prompted for 
        via tio.getOption()
        """
        if key == None:
            options = self.getList()
            i,key = getOption("Material",options)
        #
        #        Scan through the data looking for the key thsi will crash if the format of the database is wrong.
        try:
            for line in DataBase:
                if not line.startswith("#") or len(line.strip()) == 0:
                    token = line.split()
                    if token[0].strip() == key:                # Got whole key

                        if token[1].strip().lower() == "formula":
                            formula = int(token[2].strip())    # Formula type
                        else:
                            raise IOError("Material.getMaterial: failed to find formula i {0:s}".format(line))

                        if token[3].strip().lower() == "range": 
                            ragn = []
                            rl = float(token[4].strip())       # Lower range
                            ragn.append(rl)
                            rh = float(token[5].strip())       # Upper range
                            ragn.append(rh)
                        else:
                            raise IOError("Material.getMaterail: failed to finf range in {0:s}".format(line))

                        if token[6].strip().lower() == "coef": # Start of coef
                            coef = []
                            for c in token[7:]:                # May be any numner
                                cf = float(c.strip())
                                coef.append(cf)
                                                               # Success, return
                            return Material(key,formula,ragn,coef)
                        else:
                            raise IOError("Material.getMaterail: failed to finf range in {0:s}".format(line))
        except (IOError,ValueError):
            raise SyntaxError("MaterialData.getMaterial: syntax error on line [{0:s}]".format(line))
            
    
        raise ValueError("MaterialData.getMaterial: Material {0:s} not found.".format(key))


    #
    #
    def getIndex(self,key = None ):
        """
        Method to get the index of glass specified by name key name  
        glass in MaterialData, if None it will be prompted for via tio.getOption()
        param key the glass key name
        return RefrativeIndex, None of key not found.
        """
        
        #           Trap special case of "air"
        if key != None and key.lower().strip() == "air":
            return wl.AirIndex()
        else:
            material = self.getMaterial(key)
            return material.getIndex()

            

#           
class Material(object):
    """
    Class to hold a meging type being name, formula and coefficents in the form held in 
    RefractiveIndex.org
    """

    #       
    
    def __init__(self,name,formula,wrange,coef):
        """
        Constructor to for a material
        param name string name of material, usually glass key
        param formula int formula type, only 1 and 2 suppoted at the moment
        param coef list of coefficeint is same syntax ar RefrativeIndex.info
        """
        self.name = name
        self.formula = int(formula)
        self.wrange = list(wrange)
        self.coef = list(coef)
    #    
    #
    def __repr__(self):
        """
        Implement repr() to return string of informatiom
        """
        s = "Material: name : {0:s} formula : {1:d} range: {2:s} coef: {3:s}\n".\
            format(self.name,self.formula,str(self.range),str(self.coef))
        return s

    #
    #
    def getIndex(self):
        """
        Method to get the RefrativeIndex for the current Material.
        Suppoted:      Air as special case of formula 6
                       Sellmeir formula 1 and 2
                       Polynomial formula 3
                       Cauchy formula 5
                       GasIndex formula 6
    
        """
        if  self.name.lower().strip() == "air":     # Trap air as special case
            return wl.AirIndex()
        else:
            return wl.InfoIndex(self.formula,self.wrange,self.coef,self.name)



