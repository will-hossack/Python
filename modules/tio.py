""" 
The tio package is a set of terminal input methods to read the  built in Python types plus Vector2d and Vector3d classes.\
The packages formats prompts, give defaults, does range checking and general sanity checking. All will re-prompt on errors.

There is also a simple internal command handler and journal facility which is being developed.

Author:   Will Hossack, The University of Edinburgh 
"""

from os import getenv
from os.path import expanduser
from datetime import datetime
import math
import cmath
from vector import Vector2d,Vector3d
import sys

#               Internal globals 
__journalFile = None
__tioinput = sys.stdin
__tiooutput = sys.stdout
__tioerr = sys.stderr
__tioesc = "%"



def getString(prompt,default = None):
    """
    Read a string from the terminal with no processing or evaluation. 
    param prompt string the prompt to be displayed.
    param default the default string (defaults to None)
    returns string with no processing, but with leading and trailing white space removed.

    Note: will allow a zero length string to be returned.
    """
    while True:
        val = __getInput(prompt,default)
        try:
            s = str(val)                         # Force to be a string
            return s.strip()                     # Strip off leading and trailing white space
        except:                                  # Catch everything
            __tioerr.write("tio.getString.error: string conversion failed, try again\n")


def getFloat(prompt,default = None ,min = float("-Inf") ,max = float("Inf")):
    """
    Read a float from the terminal with optional default and range checking.
    param prompt the prompt string to be displayed
    param default float (defaults to None)
    param min float min value accepted (defaults to -Inf)
    param max float max value accepted (defaults to Inf)
    returns float in specified range.

    Note: if response is a string it will be evaluated.
    """
    if default != None:
        default = float(default)
    while True:  
        val = __getInput(prompt,default)
        try:                                             # Work out what happened
            if isinstance(val,str):                      # If str eval what been given
                val = eval(val)
            fval = float(val)                            # force to float
            if fval >= min and fval <= max:              # test range
                return fval
            else:
                __tioerr.write("tio.getFloat.error: {0:6.3e} outside range {1:6.3e} to {2:6.3e}\n".\
                format(fval,min,max))
        except:
           __tioerr.write("tio.getFloat.error: conversion of {0:s} failed\n".format(str(val)))


#        
def getInt(prompt,default = None ,min = -sys.maxint - 1 ,max = sys.maxint):
    """
    Read an int from the terminal with optional default and range checking. This will accept decimal,
    but also also binary (prefix 0b) , oct (prefix 0o) , hex prefix (0x)
    param prompt string the prompt string to be displayed.
    param default int (defaults no None)
    param min int min value accepted (defaults to -sys.maxint - 1) 
    paramfloat max value, (defaults to -sys.maxin
    return int in specified range.

    Note: if response is a string it will be evaluated.
    """
    while True:  
        val = __getInput(prompt,default)                 # Get input
        try:                                             # Work out what happened
            if isinstance(val,str):                      # If str eval what been given
                val = eval(val)
            if isinstance(val,float):
                if round(val) == val:                    # It is a int
                    val = int(round(val))
            if isinstance(val,int):                      # Check we have an int
                if val >= min and val <= max:              # test range
                    return val
                else:
                    __tioerr.write("tio.getInt.error: {0:d} outside range {1:d} to {2:d}\n".\
                                   format(ival,min,max))
            else:
                __tioerr.write("tio.getInt.error: conversion of {0:s} not an integer.\n".format(str(val)))
        except:
            __tioerr.write("getInt.error: conversion of {0:s} failed\n".format(str(val)))
#
#       
def getBool(prompt,default = None):
    """
    Read a logical boolen from the terminal with checking. 
    It will accept: yes / no , true / false in lower and upper case,  1 / 0 or any logical expression
    prompt string the prompt to be displayed
    default boolean logical (may be None)
    returns boolean True/False

    Note: if response is a string it will be evaluated.
    """
    while True:
        val = __getInput(prompt,default)
        try:
            if isinstance(val,bool):
                return val
            if isinstance(val,str):                      # if str 
                bval = val.lower().strip()               
                if bval.startswith("yes") or bval.startswith("true"):
                    return True
                if bval.startswith("no") or bval.startswith("false"):
                    return False

            val = eval(val)                               # do eval to for expression
            bval = bool(val)                              # Try and convert to bool
            return bval
        except:
            __tioerr.write("tio.getBool.error: conversion of {0:s} failed\n".format(str(val)))
#
#       
#
def getComplex(prompt, default = None, maxabs = float("Inf")):
    """
    Read a complex from the terminal with optional defaul and range checking.
    prompt string the prompt to be displayed.
    param default complex the default (may be None)
    param maxabs float maximum abs (defaults to float("Inf"))
    returns complex in specified range (if imaginary part is zero, it will still retun a complex)
    """
    while True:
        val = __getInput(prompt,default)                       # Get input
        try:                                                   # Work out what happened
            if isinstance(val,str):                            # If str eval what been given
                val = eval(val)
            if isinstance(val,list) or isinstance(val,tuple): # Convert from list or tuple if needed
                cval = complex(val[0],val[1])
            else:
                cval = complex(val)                            # Convert to complex
            if abs(cval) <= maxabs:              # test range
                return cval
            else:
                __tioerr.write("getComplex.error: abs {0:6.3e} outside range {1:6.3e}\n".\
                      format(abs(cval),maxabs))
        except:
            __tioerr.write("getComplex.error: conversion of {0:s} failed\n".format(str(val)))

#
# 
def getVector3d(prompt, default = None, maxabs = float("Inf")):
    """
    Read a Vector3d form the module vector from the terminal with checking.
    Format from terminal may be 'x,y,z'   OR   '[x,y,z]',  also each componet will be evaluated.
    prompt string the prompt to be displayed
    default Vector3d the default value (may be None)
    maxabs float maximum absolutle value of the Vector3d, defaults to float("Inf")
    returns a Vector3d
    """
    while True:
        val = __getInput(prompt,default)
        try:
            if isinstance(val,str):          # Its a string
                val = eval(val)              # Eval list
                vec = Vector3d(val)
            else:
                vec = Vector3d(val)          # default vector3d

            if abs(vec) <= maxabs:
                return vec                  #  Success
            else:
                __tioerr.write("tio.getVector3d.error: abs {0:6.3e} outside range {1:6.3e}\n".\
                               format(abs(vec),maxabs))
        except:
            __tioerr.write("tio.getVector3d.error: conversion of {0:s} failed\n".format(str(val)))

#
#
def getVector2d(prompt, default = None, maxabs = float("Inf")):
    """
    Read a Vector2d form the module vector from the terminal with checking.
    Format from terminal may be 'x,y'   OR   '[x,y]',  also each componet will be evaluated.
    prompt string the prompt to be displayed
    default Vector2d the default value (may be None)
    maxabs float maximum absolutle value of the Vector2d, defaults to float("Inf")
    returns a Vector2d
    """
    #
    while True:
        val = __getInput(prompt,default)
        try:
            if isinstance(val,str):          # Its a string
                val = eval(val)         # Eval list
                vec = Vector2d(val)
            else:
                vec = Vector2d(val)  # default vector3d

            if abs(vec) <= maxabs:
                return vec                  #  Success
            else:
                __tioerr.write("tio.getVector2d.error: abs {0:6.3e} outside range {1:6.3e}\n".\
                               format(abs(vec),maxabs))
        except:
            __tioerr.write("getVector2d.error: conversion of {0:s} failed\n".format(str(val)))

#
#       
def getExpandedFilename(name):
    """
    Method to expand a filename and process environmental variable in $env or ~username prefix to a filename. 
    param name string with original name, assumed to contains NO leading white spaces
    returns string the expanded filename.

    Typical input is $HOME/data.data or ~fred/data.dat, where $HOME is env name and fred is username.
    """
    if name.startswith("$") or name.startswith("~") :
        i = name.find("/")               # Name seperator
        if i < 0 :
            i = len(name)                # No seperator
        postname = name[i:]              # name after prefix, may be null
        if name.startswith("$"):         # Env name
            envname = name[1:i]
            prename = getenv(envname)
            if prename == None :        # No name, try and send back what we can
                return postname[1:]
        else:
            user = name[0:i]            # user supplied in ~username
            prename = expanduser(user)  # try and expand it
        return prename + postname       # Return processed name
    else:
        return name                     # No pre characters, just return name
#
#
def openFile(prompt,key = "r",defaulttype = None, defaultname = None):
    """
    Method to open a text file with sanity checking with re-prompt on failure.
    param prompt string the prompt to be displayed.
    param key string the mode as used by open, default is "r" (read)
    param defaultype deafult extension which will be added if not supplied.
    param  defaultname string the default filename (default to None)
    returns the opened filestream

    The file names is processded to expand environmental variable and user names\
    so for example $ENV/dir/file.data or ~user/dir/file.data are expanded
    """
    while True:
        val = getString(prompt,defaultname)
        filename = getExpandedFilename(val)                       # Expand to process env/user
        if defaulttype != None:
            if not filename.endswith(defaulttype):
                filename += "." + defaulttype
        try:
            filestream = open(filename,str(key))                  # try and open
            return filestream
        except IOError:                                           # Error
            __tioerr.write("openFile.error: failed to open {0:s} with key {1:s}\n".\
                           format(filename,str(key)))

def tprint(*args):
    """
    Simply alternative to print that will print to the sysout and also journal if there is a journal file open.
    Output to the journal file be prefixed with a comment character
    param, argument list, each will be conveterd to a str() and concatinated to a single string. 

    Also  "\n" will be appended if not present and the print buffer will be flushed.
    """

    #               Form output string by appending str() of each argument
    string = ""
    for a in args:
        string += str(a)

    if not string.endswith("\n"):       # Add automatic newline if needed
        string += "\n"

    __tiooutput.write(string)           # Write string and flush output buffer
    __tiooutput.flush()          
    
    if __journalFile != None:           # Journal file open, so write string, but with prefix of "# "
        __journalFile.write("# " + string)
    

#
#
def getOption(prompt,options,default = None):
    """
    Method to get a choice of options from a supplied list with error checking.
    para prompt string the prompt to be displayed.
    param options [] the options list assumed to be list of strings.
    default int the default option (int in range 0 -> < len(options)) (default to None)
    returns list [] of (opt,option[opt])
    
    Each option is tested for existance and uniquness. It will fail an re-prompt as required. There is also\\
    a simple internal 'help' option that pints the list of options.
    """
    if default == None:
        defaultOption = None
    else:
        if default >= 0 and default < len(options):
            defaultOption = options[default]
        else:
            defaultOption = None
    while True:
        val = getString(prompt,defaultOption)       # Get input
        #
        #                 Search for option (also check for uniquness)
        opt = -1
        i = 0
        for o in options:
            if o.startswith(val):                            # found option
                if opt == -1:                                # first find
                    opt = i
                else:
                    opt = -2                                 # Not unique
            i += 1

        if opt >= 0:                                         # success, one found return
            return opt,options[opt]

        #
        #                Deal with help.
        if val.startswith("help"):
            pl = ""
            for o in options:
                pl += " [{0:s}]".format(o)
            
            tprint("Options are : {0:s}".format(pl))
        elif opt == -1:                                     # Unknown option
            __tioerr.write("tio.getOption.error: Invalid option {0:s}, help for list of options.\n".format(val))
        elif opt == -2:                                     # Non-unique option
            __tioerr.write("tio.getOption.error: Non-unique option {0:s}, help for list of options.\n".format(val))
                  
#
#
def __formatPrompt(prompt,default = None):
    """
    Internal method to format the prompt and add a default if included.
    """
    if default != None:                     # Add a default if given
        if isinstance(default,float):       # For float
            prompt += " (default : {0:5.3g}) : ".format(default)
        elif isinstance(default,int):       # for int
            if isinstance(default,bool):    # for bool (subclass of int)
                prompt += " (default : {0} ) : ".format(default)
            else:
                prompt += " (default : {0:d}) : ".format(default)
        elif isinstance(default,complex):   # For complex
            prompt += " (default : {0.real:5.3e} {0.imag:+5.3e}j) : ".\
                      format(default)
        else:                               # assume is string or obeys str
            prompt += " (default : {0:s}) : ".format(str(default))
            
    #             Add a " : " at end of there is not at least one ":"
    #
    if prompt.find(":") < 0 :
        prompt += " : "

    return prompt

#
#
def __getInput(prompt,default):
    """
   Internal method to get the response from the terminal and apply default if given.\
    This method will strip comments denoted by #  but no other processing.

    This used a .readline() from a input steeam and not input() or raw_input() so will work with both P2 and P3
    """
    p = __formatPrompt(prompt,default)
    #
    while True:
        __tiooutput.write(p)
        val = __tioinput.readline()
        i = val.find("#")              # is there a comment
        if (i >= 0):                   # comment found
            val = val[0:i]             # Kill comment
        val = val.strip()              # kill white space
        if val.startswith(__tioesc):    # process tio commands
          __tiocommand(val)
        elif len(val) > 0:
            break                  # have something valid
        else:
            if default != None:
                val = default      # Take default
                break
    
    if __journalFile != None:      # Journal to stream open 
        __journalFile.write("{0:s}    # {1:s}\n".format(str(val),p))
    return val
        

#
#
def setJournal(filename = None):
    """
    Method of open / close a journal file that records prompts and commands typed to a text file.
    param filename string the name of the journal file, if None will close any current open journal file.
    """
    global __journalFile

    if filename == None:             # No file give.
        if __journalFile != None:    # Close Journal file if open
            __journalFile.write("#             closed at {0:s}\n".format(str(datetime.now()))) 
            __journalFile.close()
            __journalFile = None     # Null journal file
            tprint("tio.info: Journal off.")
        return                       # All finished

    #                                File given, so try and open it
    fn = getExpandedFilename(filename)
    if not fn.endswith("tio"):
        fn += ".tio"
    try:
        __journalFile = open(fn,"w")
        __journalFile.write("#             tio Journal file\n")
        __journalFile.write("#             opened at {0:s}\n".format(str(datetime.now()))) 
        tprint("tio.info: Journal on.")
    except IOError:
        __tioerr.write("setJournal.error: file open of {0:s} failed\n".format(fn))
        if getBool("Manually open journal file",False):
            __journalFile = openFile("Journal File","w","tio")
        else:
            __journalFile = None     

       
        
#   
#
def __tiocommand(cmd):
    """ Internal command handler, limited use at the moment, but will be expanded
    """
    cmd = cmd[1:].strip().lower()          # Remove % and clean up
    if cmd.startswith("beep"):
        __tiooutput.write("\a")
    elif cmd.startswith("exit"):
        sys.exit(0)
    elif cmd.startswith("journal"):
        tokens = cmd.split()
        filename = tokens[1].strip()
        setJournal(filename)
    elif cmd.startswith("nojournal"):
        setJournal()
    else: 
        __tioerr.write("toi.command error: unknown command {0:s}, ignored.\n".format(cmd))
 

