#==============================================================================
# These are some scripts for testing the functionality of LSDMappingTools
#==============================================================================
# -*- coding: utf-8 -*-
"""
Created on Tue May 05 14:08:16 2015

@author: smudd
"""

import numpy as np
import LSDPlottingTools as LSDP

def TestNewMappingTools2():
    DataDirectory = "T://analysis_for_papers//Manny_idaho//"
    Filename = "TestIdaho.bil"
    NewFilename = "TestIdaho_after2.bil"
    ThisFile = DataDirectory+Filename
    NewFilename = DataDirectory+NewFilename
    
    #FigFormat = 'svg'
    #FigFileN= 'Sorbas_chi.svg'
    #FigFileName= DataDirectory+FigFileN
    

    #Plot the basin over the elevation    
    #tcmapcolorbarlabel = "Elevation (m)"
    #tcmap = 'jet'
    #tcmapcolorbarlabel='Chi'
    #clim_val = (50,300)
    #LSDP.BasicDensityPlotGridPlot(ThisFile,tcmap,tcmapcolorbarlabel,clim_val)
    
    #get the nodata values
    NData = LSDP.getNoDataValue(ThisFile)
    
    print "NoData is: " + str(NData)
    
    # get the data as an array
    array = LSDP.ReadRasterArrayBlocks(ThisFile)
    
    driver_name = "ENVI"
    LSDP.array2raster(ThisFile,NewFilename,array,driver_name, -9999)
    
    #get the nodata values
    NData2 = LSDP.getNoDataValue(NewFilename)
    
    print "NoData 2 is: " + str(NData2)    



if __name__ == "__main__":
    #fit_weibull_from_file(sys.argv[1]) 
    TestNewMappingTools2() 