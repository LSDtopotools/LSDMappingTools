# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 17:32:40 2015

@author: smudd
"""

#==============================================================================
# These are some scripts for testing the functionality of LSDMappingTools
#==============================================================================
# -*- coding: utf-8 -*-
"""
Created on Tue May 05 14:08:16 2015

@author: smudd
"""

import numpy as np
import LSDMappingTools as LSDmt

def TestMappingTools():
    #DataDirectory = "C://basin_data//Chile//lat26p0//"
    #Filename = "Site_lat26p0_UTM19_DEM_FILL.bil"
    #DataDirectory = "C://basin_data//Model_results//June2015_Results//HighK//" 
    #Filename = "InitialForCRN.asc"
    #Filename = "CRNvariable_long_0_0_1var128.asc"
    DataDirectory = "C://basin_data//Babault2//"
    #Filename = "Julien_DEM_HS.bil"    
    Filename = "Julien_DEM_Q.bil"
    
    DrapeFileName = "Julien_DEM_Q.bil"
    ThisFile = DataDirectory+Filename
    DrapeFile =DataDirectory+DrapeFileName 
    
    #data = LSDmt.ReadRasterArrayBlocks(ThisFile)
    #print "Data is: "
    #print data
    
    #print data.shape
    
    NDV, xsize, ysize, GeoT, Projection, DataType = LSDmt.GetGeoInfo(ThisFile)
    print "NDV: " + str(NDV)
    print "xsize: " + str(xsize)
    print "ysize: " + str(ysize)
    print "GeoT: " + str(GeoT)
    print "Projection: " + str(Projection)
    print "DataType: " + str(DataType)
    
    CellSize,XMin,XMax,YMin,YMax = LSDmt.GetUTMMaxMin(ThisFile)
    print "CellSize: " + str(CellSize)
    print "XMin: " + str(XMin)
    print "XMax: " + str(XMax)
    print "YMin: " + str(YMin)
    print "YMax: " + str(YMax)    
    
    
    tcmap = 'autumn'
    tcmapcolorbarlabel='Discharge in $m^3$ per year'
    clim_val = (0,0)
    LSDmt.BasicDensityPlot(ThisFile,tcmap,tcmapcolorbarlabel,clim_val)
    LSDmt.LogStretchDensityPlot(ThisFile,tcmap,tcmapcolorbarlabel,clim_val)
    
    #LSDmt.DrapedPlot(ThisFile,DrapeFile)

if __name__ == "__main__":
    #fit_weibull_from_file(sys.argv[1]) 
    TestMappingTools() 