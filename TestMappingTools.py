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
    DataDirectory = "C://basin_data//Model_results//" 
    Filename = "CRNvariable_long_0_0_1var128.asc"
    DrapeFileName = "CRNvariable_long_0_0_1var128_erosion.asc"
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
    
    #LSDmt.BasicDensityPlot(DrapeFile)
    LSDmt.DrapedPlot(ThisFile,DrapeFile)

if __name__ == "__main__":
    #fit_weibull_from_file(sys.argv[1]) 
    TestMappingTools() 