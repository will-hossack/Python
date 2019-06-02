"""
Classes to handle of materail database, mainly glasses, but
also gases and liquids where are held in the asci databases using the same 
conventions as in RefractiveIndex.org. These class just handes the data bases of the materials, the refreatcive 
index is handeled by MaterialIndex() and InfoIndex() in the wavelangth module.

Author:   Will Hossack, The University of Edinburgh
"""
from sys import exit
from tio import getExpandedFilename,getOption
#                    
#

DataBaseFile = "$LENS/materials.data" 
DataBase = None

class MaterialData(object):
    """
    Class to readin the material index database and lookup materials.
    """
    #
    #
    def __init__(self,filename = None):
        """
        Constructor to open the materials definition file and read in the  data to self.data as a list of lines. 
        This is not normally called by the user but auotomatically in the background
        when the first referative index is looked up.
        The database is only loaded once for efficiency.
        param filename the database (defaults to package default) 
        If the database is missing an OSError is raised.
        """
        global DataBase
        if DataBase == None :                # Database not loaded, so this is first call
            if filename == None:             # Blank name, so use package default.
                filename = DataBaseFile
            try :
                filename = getExpandedFilename(filename)   # Sort out logicals
                filestream = open(filename,"r")
                DataBase = filestream.readlines()          # Read in the database
                filestream.close()
            except :
                raise OSError("MaterialData() unable to open data file --{0:s}-- PANIC STOP".format(filename))
    #
    #      
    def getList(self):
        """
        Method to get a ist of the materials keys in the database as a list of strings.
        """
        key = []
        for line in DataBase:
            if not line.startswith("#") or len(line.strip()) == 0:
                token = line.split()
                key.append(token[0].strip())     # Key is first token
        return key
        
    #    
    #
    def getMaterial(self,key = None):
        """
        Method to get a material from the loaded database by key.
        parm key name, (usually glass name) if None, then it will be prompted for 
        via tio.getOption()
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
                    
            #raise UserWarning("MaterialData.getMaterial: Material {0:s} not found.".format(key)) # Did not find the key  
            return Material("NotValid",0,[0,0],[0])                                         
                        
        except (OSError):
            raise SyntaxError("MaterialData.getMaterial: syntax error on line [{0:s}]".format(line))
            
    
                  

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


