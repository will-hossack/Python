"""
Classes to handle of materail database, mainly glasses
"""
from tio import getExpandedFilename,getOption
from importlib import resources

#         Open the "material.data" file and hold in global
try:
    file =resources.open_text("optics.lenses","materials.data")
    DataBase = file.readlines()
    file.close()
except:
    raise OSError("MaterialData() unable to open Material file -- PANIC STOP")
    
#                    
#
#import os
#directory,filename = os.path.split(__file__)
#print("Director is : " + directory)

#DataBaseFile = "$LENS/materials.data" 
#DataBase = None

class MaterialData(object):
    """
    Class to manage the Materail database. This is preloaded from within the package at
    startuo but subsequent calle to this class can be used to change this database.
    """
    #
    #
    def __init__(self,filename = None):
        """
        

        :param filename:  the new database (defaults to package default) 
        :type filename: str

        If the database is missing an OSError is raised
        """
        global DataBase
        if filename != None:            # Assume a new database
            try :
                filename = getExpandedFilename(filename)   # Sort out logicals
                filestream = open(filename,"r")
                DataBase = filestream.readlines()          # Read in the database
                filestream.close()
            except :
                raise OSError("MaterialData() unable to open data file --{0:s}-- PANIC STOP".format(filename))
    
          
    def getList(self):
        """
        Method to get a list of the materials keys in the database as a list of strings.

        :retrurn: a list of avaialable material (glass) keys
        """
        key = []
        for line in DataBase:
            if not line.startswith("#") or len(line.strip()) == 0:
                token = line.split()
                key.append(token[0].strip())     # Key is first token
        return key
        
    
    
    def getMaterial(self,key = None):
        """
        Method to get a material from the loaded database by key.

        :param key: name, (usually glass name) if None, then it will be prompted for via tio.getOption()
        :type key: str
        :return: The Material object

        """
        if key == None:
            options = self.getList()
            i,key = getOption("Material",options)
        #
        #        Scan through the data looking for the key; this  will crash if the format of the database is wrong.
        try:
            for line in DataBase:
                if not line.startswith("#") or len(line.strip()) == 0:
                    token = line.split()
                    if token[0].strip() == key:                # Got the key (not key MUST be exact)

                        if token[1].strip().lower() == "formula":
                            formula = int(token[2].strip())    # Formula type
                        else:
                            raise OSError("MaterialData.getMaterial: failed to find formula in {0:s}".format(line))

                        if token[3].strip().lower() == "range": 
                            ragn = []
                            rl = float(token[4].strip())       # Lower range
                            ragn.append(rl)
                            rh = float(token[5].strip())       # Upper range
                            ragn.append(rh)
                        else:
                            raise OSError("MaterialData.getMaterail: failed to find range in {0:s}".format(line))

                        if token[6].strip().lower() == "coef": # Start of coef
                            coef = []
                            for c in token[7:]:                # May be any number of coefficeints.
                                cf = float(c.strip())
                                coef.append(cf)
                                                               
                        else:
                            raise OSError("MaterialData.getMaterial: failed to find coef in {0:s}".format(line))

                        return Material(key,formula,ragn,coef)  #   Success found material
                    
            return Material("NotValid",0,[0,0],[0])                                         
                        
        except (OSError):
            raise SyntaxError("MaterialData.getMaterial: syntax error on line [{0:s}]".format(line))
            
    
                  

#           
class Material(object):
    """
    Class to hold a material type being name, formula and coefficents in the form held the same form as 
    https://refractiveindex.info/

    :param name: string name of material typically the key
    :type name: str
    :param formula: formula type, note that formula 0 = invalid
    :type formula: int
    :param wrange: range of validity of formula
    :type wrange: list[low,high]
    :param coef: list of float coefficeint is same syntax as RefrativeIndex.info
    :type coef: list[float]

    """
    
    def __init__(self,name,formula,wrange,coef):
        """
        Constructor to for a material
       
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


