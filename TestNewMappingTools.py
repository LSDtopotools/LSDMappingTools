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

def TestNewMappingTools():
    #DataDirectory = "C://basin_data//Chile//lat26p0//"
    #Filename = "Site_lat26p0_UTM19_DEM_FILL.bil"
    #DataDirectory = "C://basin_data//Model_results//June2015_Results//HighK//" 
    #DataDirectory = "T://test_clone//topodata//"
    DataDirectory = "T://analysis_for_papers//Cosmo_paper//Tibet//for_plotting//"
    #Filename = "SanBern.bil"
    Filename = "SpawnedBasin_07C13-(Q8)_SH.bil"
    #Filename = "CRNvariable_long_0_0_1var128.asc"
    #DrapeFileName = "CRNvariable_long_0_0_1var128_erosion.asc"
    DrapeFileName = "SpawnedBasin_07C13-(Q8)_BASINS.bil"
    ThisFile = DataDirectory+Filename
    DrapeFile =DataDirectory+DrapeFileName 

    
    NDV, xsize, ysize, GeoT, Projection, DataType = LSDP.GetGeoInfo(ThisFile)
    print "NDV: " + str(NDV)
    print "xsize: " + str(xsize)
    print "ysize: " + str(ysize)
    print "GeoT: " + str(GeoT)
    print "Projection: " + str(Projection)
    print "DataType: " + str(DataType)
    
    CellSize,XMin,XMax,YMin,YMax = LSDP.GetUTMMaxMin(ThisFile)
    print "CellSize: " + str(CellSize)
    print "XMin: " + str(XMin)
    print "XMax: " + str(XMax)
    print "YMin: " + str(YMin)
    print "YMax: " + str(YMax)    
    
    
    tcmap = 'autumn'
    tcmapcolorbarlabel='Topographic shielding'
    clim_val = (0.72,1)
    #LSDP.BasicDensityPlot(ThisFile,tcmap,tcmapcolorbarlabel,clim_val)
    #LSDP.DrapedPlot(ThisFile,DrapeFile)
    
    LSDP.BasicDensityPlotGridPlot(ThisFile,tcmap,tcmapcolorbarlabel,clim_val)

if __name__ == "__main__":
    #fit_weibull_from_file(sys.argv[1]) 
    TestNewMappingTools() 