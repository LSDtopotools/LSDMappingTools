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
    
    
        
    
    #MB = LSDP.BasicMassBalance(DataDirectory, Filename1, Filename2)
    #print "Mass balance between these two time steps is: " + str(MB) + " cubic metres"

    #Mean1 = LSDP.RasterMeanValue(DataDirectory, Filename1)
    #Mean2 = LSDP.RasterMeanValue(DataDirectory, Filename2)

    #print "The mean values of the two rasters are: " + str(Mean1) +" and "+ str(Mean2)
    
    # now try the swath plotting
    axis = 0
    LSDP.SwathPlot(DataDirectory, Filename1, axis)
    
    axis = 1
    LSDP.SwathPlot(DataDirectory, Filename1, axis)

if __name__ == "__main__":
    #fit_weibull_from_file(sys.argv[1]) 
    #TestNewMappingTools2() 
    TestNewMappingTools_asc()
    #FloodThenHillshade()
    #FixStupidNoData()
    