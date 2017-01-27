# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 10:17:29 2017

@author: smudd
"""

from __future__ import print_function
import sys
import os
from glob import glob
import LSDPlottingTools as LSDPT

def check_directory(directory):
    """This checks to see if a directory exists.
    
    Args: 
        directory (str): the path to the directory
        
    Return:
        bool: True if the directory exists
    """
    
    this_bool = False
    if not os.access(directory,os.F_OK):
        print("I can't find the directory! You wanted: "+directory)
    else:
        print("I found the directory: "+directory)
        this_bool = True
        
    return this_bool
    
def check_file_name(file_name):
    """This checks to see if a file exists.
    
    Args: 
        file_name (str): the name of the file with extension and path
        
    Return:
        bool: True if the directory exists
    """
    
    this_bool = False
    if not os.access(file_name,os.F_OK):
        print("I can't find the file! You wanted: "+file_name)
    else:
        print("I found the directory: "+file_name)
        this_bool = True
        
    return this_bool    

def read_parameter_file(filename):
    """This reads a parameter file
    
    Args:
        filename (str): The name of the parameter file (with path and extension)
    
    Author: SMM
    """
    
    file_there = check_file_name(filename)
    
    if file_there:
        
        this_file = open(filename, 'r')
        lines = this_file.readlines()

        param_dict = {}
        # loop through the file harvesting the thge keys to the dicts
        for line in lines:
            # get rid of the control characters
            this_line =LSDPT.RemoveEscapeCharacters(line)
        
            # get the parameters
            split_line = this_line.split(':')
            
            param_dict[split_line[0]]=split_line[1] 
            
        print(param_dict)
        
    else:
        print("I did not find your parameter file so I can't read it.")

    print         
    

#=============================================================================
# This is the main function that runs the whole thing 
#=============================================================================
def main(argv):
 
    # If there are no arguments, send to the welcome screen
    if not len(sys.argv) > 1:
        print_welcome()
        sys.exit()

    # Get the arguments
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-wd", "--working_directory",type=str, default=os.path.dirname(os.getcwd()), 
                        help="The working directory where you have your files.")
    parser.add_argument("-fp", "--file_prefix",type=str, default="Raster", 
                        help="The prefix of your files.") 
    parser.add_argument("-pf", "--parameter_file",type=str, default="Params.param", 
                        help="The name, with extension but not path, of your parameter file.")    
    args = parser.parse_args()

    working_dir = args.working_directory
    working_dir = LSDPT.AppendSepToDirectoryPath(working_dir)
    full_paramfile = working_dir+args.parameter_file
    
    print("The full parameter file is: "+full_paramfile)
    read_parameter_file(full_paramfile)
    
#=============================================================================
    
    
#=============================================================================    
# This is just a welcome screen that is displayed if no arguments are provided.
#=============================================================================
def print_welcome():

    print("\n\n=======================================================================")
    print("Hello there, I am the going to help you plot LSDTopoTools data!")
    print("You will need to tell me where the files are.")
    print("Use the -wd flag to define the working directory. If you dont do this I will assume the data is in the same directory as this script.")
    print("=======================================================================\n\n ")
#=============================================================================



if __name__ == "__main__":
    main(sys.argv[1:])