#==============================================================================
# These are some scripts for testing the functionality of LSDMappingTools
#==============================================================================
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 14:08:16 2016

@author: smudd
"""

import numpy as np
import LSDPlottingTools as LSDP

def TestMappingToolsPoints():
    #DataDirectory = "T://analysis_for_papers//Test_map_chi_gradient//results//"
    DataDirectory = "T://students//Wainwright//"
    #Filename = "Mandakini_OutletList.csv"
    Filename = "Bids_DD.csv"
    
    #ExportName = "ShapeTest"
    
    fname = DataDirectory+Filename
    #Exp_fname = DataDirectory+ExportName

    thisPointData = LSDP.LSDMap_PointData(fname)   

    thisPointData.GetParameterNames(True)
    thisPointData.GetLongitude(True)
    
    print "Hey buddy, the province is: "
    thisPointData.QueryData("province",True)
    
    print "Hey buddy, the gma is: "
    thisPointData.QueryData("gma",True)    
    
    thisPointData.TranslateToReducedShapefile(fname)
    thisPointData.TranslateToReducedGeoJSON(fname)
    
def TestMappingToolsLassoCSV(): 
    
    #DataDirectory = "T://analysis_for_papers//Test_map_chi_gradient//results//"
    DataDirectory = "T://analysis_for_papers//Indus//NewChi//"
    LSDP.ConvertAllCSVToGeoJSON(DataDirectory)




if __name__ == "__main__":
    #TestMappingToolsPoints() 
    TestMappingToolsLassoCSV()
    