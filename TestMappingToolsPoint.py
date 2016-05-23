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
    DataDirectory = "T://analysis_for_papers//Test_map_chi_gradient//results//"
    Filename = "Mandakini_OutletList.csv"
    ExportName = "ShapeTest"
    
    fname = DataDirectory+Filename

    thisPointData = LSDP.LSDMap_PointData(fname)   

    thisPointData.GetParameterNames(True)
    thisPointData.GetLongitude(True)
    
    
    





if __name__ == "__main__":
    TestMappingToolsPoints() 
    