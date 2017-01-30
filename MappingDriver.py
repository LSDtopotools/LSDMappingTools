# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 10:17:29 2017

@author: smudd
"""

from __future__ import print_function
import sys
import os
import LSDPlottingTools as LSDPT


    

#=============================================================================
# This is the main function that runs the whole thing 
#=============================================================================
def main(argv):
 
    # If there are no arguments, send to the welcome screen
    if not len(sys.argv) > 1:
        full_paramfile = print_welcome()
        #sys.exit()
    else:

        # Get the arguments
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument("-pf", "--parameter_file",type=str, default="Params.param", 
                            help="The name, with extension and path, of your parameter file.")    
        args = parser.parse_args()

        full_paramfile = args.parameter_file
    
    print("The full parameter file is: "+full_paramfile)
    
    # Now make a plotting driver object
    PD = LSDPT.LSDMap_PlottingDriver(full_paramfile)
    
    print(PD.FilePrefix)
    print(PD.FilePath)
    
    print("The plotting switches are: ")
    print(PD.plotting_switches)
    
    # Now make some plots!!
    PD.plot_data()
    
#=============================================================================
    
    
#=============================================================================    
# This is just a welcome screen that is displayed if no arguments are provided.
#=============================================================================
def print_welcome():

    print("\n\n=======================================================================")
    print("Hello there, I am the going to help you plot LSDTopoTools data!")
    print("You will need to tell me where the parameter file is.")
    print("Use the -wd flag to define the working directory.")
    print("If you dont do this I will assume the data is in the same directory as this script.")
    print("=======================================================================\n\n ")

    from Tkinter import Tk
    from tkFileDialog import askopenfilename

    Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file
    
    return filename    
   
    
    
#=============================================================================


if __name__ == "__main__":
    main(sys.argv[1:])