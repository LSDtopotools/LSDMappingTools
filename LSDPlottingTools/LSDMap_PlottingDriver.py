# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 10:20:18 2017

@author: smudd
"""


import LSDMap_OSystemTools as LSDOst
import os
import sys


class LSDMap_PlottingDriver(object):
    
    # The constructor: it needs a filename to read    
    def __init__(self,FileName):
        """This is the LSDMap_PlottingDriver object. 
        
        It manages plotting functions via a parameter file. 
        The parameter file syntax is similar to LSDTopoTools driver file syntax       
        
        Args:
            Filename (str): The name of the parameter file (with path and extension).
            
            Parameter files are made up of keywords followed by semicolons. 
            The keyword file_prefix *MUST* be included.
            
        Author: SMM
        """
        
        print("The filename is: "+FileName)         
        
        # This gets the path of the parameter file
        file_path = LSDOst.GetPath(FileName)
        self.FilePath = file_path        
        print("The path is: "+file_path)        


        #See if the parameter file exists
        if os.access(FileName,os.F_OK):
            this_file = open(FileName, 'r')
            lines = this_file.readlines()

            param_dict = {}
            for line in lines:            
                this_line =LSDOst.RemoveEscapeCharacters(line)

                # get the parameters
                split_line = this_line.split(':')
            
                param_dict[split_line[0]]=split_line[1] 

            self.parameter_dict = param_dict
        else:
            print("Parameter file does not exit! I cannot proceed. Check your filenname.")
            sys.exit()

        # Now get out the file prefix
        self.FilePrefix = "None"
        if "file_prefix" in param_dict:
            this_prefix = param_dict["file_prefix"]          
            self.FilePrefix = LSDOst.RemoveWhitespace(this_prefix)
            
            print("The file prefix is: "+self.FilePrefix)
        else:
            print("Parameter file does does not contain the keyword 'file_prefix': Check your parameter file. I must exit now.")
            sys.exit()            
        
        self.hs_fname = file_path+os.sep+self.FilePrefix+"_hs.bil"
        self.chi_csv_fname = file_path+os.sep+self.FilePrefix+"_MChi_segmented.csv"
        self.basin_csv_fname = file_path+os.sep+self.FilePrefix+"_AllBasinsInfo.csv"        
