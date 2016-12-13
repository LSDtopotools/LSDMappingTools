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

def ChiMappingToolsTest():
    DataDirectory = "/home/smudd/SMMDataStore/analysis_for_papers/Meghalaya/chi_analysis/"
    #DataDirectory = "T:\\analysis_for_papers\\Meghalaya/chi_analysis\\"
    Filename = "Mega_for_chi.bil"
    HSFilename = "Mega_for_chi_hs.bil"
    
    DEMname = DataDirectory+Filename
    HSname = DataDirectory+HSFilename
    
    FigName = DataDirectory+'Image3.pdf'
    
    #FigFormat = 'svg'
    #FigFileN= 'Sorbas_chi.svg'
    #FigFileName= DataDirectory+FigFileN
    

    #LSDP.BasicDrapedPlotGridPlot(DEMname,HSname, 'gray','gray',
    #                        'Elevation in meters',(0,0),
    #                        0.4,FigName,'pdf')  
    
    EPSG_string = LSDP.GetUTMEPSG(DEMname)
    print "EPSG string is: " + EPSG_string
    
    
    #Plot the basin over the elevation    
    #tcmapcolorbarlabel = "Elevation (m)"
    #tcmap = 'jet'
    #tcmapcolorbarlabel='Chi'
    #clim_val = (50,300)
    #LSDP.BasicDensityPlotGridPlot(ThisFile,tcmap,tcmapcolorbarlabel,clim_val)

    # Now load some point data
    Filename = "Test_chi.csv"

    fname = DataDirectory+Filename
    thisPointData = LSDP.LSDMap_PointData(fname) 
    
    # convert to easting and northing
    [easting,northing] = thisPointData.GetUTMEastingNorthing(EPSG_string)
    
    print easting
    print northing

if __name__ == "__main__":
    ChiMappingToolsTest()

    