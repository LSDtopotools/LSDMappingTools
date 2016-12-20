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
    #DataDirectory = "/home/smudd/SMMDataStore/analysis_for_papers/Meghalaya/chi_analysis/"
    #DataDirectory = "T:\\analysis_for_papers\\Meghalaya/chi_analysis\\"
    DataDirectory = "C:\\Vagrantboxes\\LSDTopoTools\\Topographic_projects\\Meghalaya\\"
    Filename = "Mega_clip.bil"
    HSFilename = "Mega_clip_hs.bil"
    
    DEMname = DataDirectory+Filename
    HSname = DataDirectory+HSFilename
    
    FigName = DataDirectory+'TestChiFull2.png'
    ChiName = DataDirectory+'Mega_clip_MChiSegmented.csv'
    
    #FigFormat = 'svg'
    #FigFileN= 'Sorbas_chi.svg'
    #FigFileName= DataDirectory+FigFileN
    elevation_threshold = 1
    
    #LSDP.BasicChiPlotGridPlot(DEMname,HSname,ChiName, 'gray','gray',
    #                        '$k_{sn}$',(0,0),
    #                        0.4,FigName,'png',elevation_threshold)  
    #FigName2 = DataDirectory+'TestChannelMap.png'    
    #LSDP.BasicChannelPlotGridPlotCategories(DEMname,HSname,ChiName, 'gray','gray',
    #                        '$Channel$',(0,0),
    #                        0.4,FigName2,'png',elevation_threshold,'source_key')      
    
    #FigName3 =  DataDirectory+'ChiProfiles.png'    
    #LSDP.ChiProfiles(ChiName, FigName3,'png',elevation_threshold)     
    
    FigName4 =  DataDirectory+'ChiStackProfiles.png'  
    first_basin = 0
    last_basin = 10
    LSDP.StackedChiProfiles(ChiName, FigName4,'png',elevation_threshold,first_basin,last_basin)  
    #LSDP.BasicDrapedPlotGridPlot(DEMname,HSname, 'gray','gray',
    #                        'Elevation in meters',(0,0),
    #                        0.4,FigName,'pdf')  
    
    #EPSG_string = LSDP.GetUTMEPSG(DEMname)
    #print "EPSG string is: " + EPSG_string
    
    
    #Plot the basin over the elevation    
    #tcmapcolorbarlabel = "Elevation (m)"
    #tcmap = 'jet'
    #tcmapcolorbarlabel='Chi'
    #clim_val = (50,300)
    #LSDP.BasicDensityPlotGridPlot(ThisFile,tcmap,tcmapcolorbarlabel,clim_val)

    # Now load some point data
    #Filename = "Mega_for_chi_MChiSegmented.csv"

    #fname = DataDirectory+Filename
    #thisPointData = LSDP.LSDMap_PointData(fname) 
    
    # convert to easting and northing
    #[easting,northing] = thisPointData.GetUTMEastingNorthing(EPSG_string)
    
    #print easting
    #print northing

if __name__ == "__main__":
    ChiMappingToolsTest()

    