""" Functions to read and write CSV file to / from numpy arrays.
Author: Will Hossack, to support Junior Hons Lab.
"""

import csv
import numpy as np

comment = "#"         # Comment delimeter


def readCSV(file,cols = None,separ=",",headerskip = 0):
    """
    Read a comma delimited cvs files of floats with secified columns if supplied.
   
    :param file: the csv file to be read
    :type file: str or file
    :param cols: truple of locicals specifiying which columns to be read Default = None (all)
    :type cols: list[bool] or None
    :param separ: field seperator, (Default = comma)
    :param headerskip: number of lines in headed to skip (Default = 0)
    :return:   two dimensional np.array of values.

    """

    #     Open the file if a name was given as a str.
    if isinstance(file,str):
        file = open(file,"r",newline='')
        
    #      open the reader with specified field separator
    reader = csv.reader(file,delimiter=separ)
    
    #       Create list and read line one at time
    data = []
    i = 0   # Nubmber of lines read
    for line in reader:
        # Skip if in of lenth 0, in header or starts with comment
        if len(line) > 0 and i >= headerskip and not line[0].startswith(comment) :  
            vals = []
            if cols == None:           # If no cols given so read all
                for t in line:
                    vals.append(float(t))
            else:
                for t,c in zip(line,cols):
                    if c:
                        vals.append(float(t))
                    
            data.append(vals)    # Add read line to data
        i += 1
    file.close()                 # Close the file
    # Now transpose and convert to np arrays.
    return np.array([*zip(*data)])



def writeCSV(file,data,cols = None):
    """
    Write CSV file with data typically supplied as a list or np.arrays.

    :param file: file or name of file, 
    :type file: file or str
    :param data: list of np.arrays.
    :type data: list of np.ndarrays
    :param cols: boolean list specifying which colums of data are to be written, defaults to None, which means all cols.
    :type cole: truple of booleans or None
    """
    #     Open the file if a name was given as a str.
    if isinstance(file,str):
        file = open(file,"w",newline='')
        
    writer = csv.writer(file)           # Make writer
    #      Read through data writing out what is required.
    for j in range(0,data[0].size):     # Number of elements
        s = []
        for i in range(0,len(data)):    # Number of cole
            if cols == None or cols[i] :
                s.append("{0:12.5e}".format(data[i][j]))

        writer.writerow(s)

    file.close()
    return data[0].size           # Number writer
