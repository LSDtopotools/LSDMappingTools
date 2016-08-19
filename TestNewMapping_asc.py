#==============================================================================
# These are some scripts for testing the functionality of LSDMappingTools
#==============================================================================
# -*- coding: utf-8 -*-
"""
Created on August 18th 2016

@author: smudd
"""

import numpy as np
import LSDPlottingTools as LSDP

def TestNewMappingTools_asc():
    DataDirectory = "T:\\analysis_for_papers\\Beaches\\"
    Filename1 = "BedThickness_050.asc"
    Filename2 = "BedThickness_100.asc"

    ThisFile = DataDirectory+Filename1


    yo = LSDP.GetRasterExtent(ThisFile)
    
    print "raster extent is: " 
    print yo
    
    
    LSDP.BasicMassBalance(DataDirectory, Filename1, Filename2)
    


if __name__ == "__main__":
    #fit_weibull_from_file(sys.argv[1]) 
    #TestNewMappingTools2() 
    TestNewMappingTools_asc()
    #FloodThenHillshade()
    #FixStupidNoData()
    