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
import matplotlib.pyplot as plt

def ChiMappingToolsTest():
    #DataDirectory = "/home/smudd/SMMDataStore/analysis_for_papers/Meghalaya/chi_analysis/"
    DataDirectory = "T:\\analysis_for_papers\\Meghalaya/chi_analysis\\"
    #DataDirectory = "C:\\Vagrantboxes\\LSDTopoTools\\Topographic_projects\\Meghalaya\\"
    Filename = "Mega_clip.bil"
    HSFilename = "Mega_clip_hs.bil"
    BasinFilename = "Mega_clip_AllBasins.bil"
    
    DEMname = DataDirectory+Filename
    HSname = DataDirectory+HSFilename
    Basinname = DataDirectory+BasinFilename 
    
    FigName = DataDirectory+'TestChiFull2.png'
    ChiName = DataDirectory+'Mega_clip_MChiSegmented.csv'
    BasinInfoName = DataDirectory+'Mega_clip_AllBasinsInfo.csv'
    
    #FigFormat = 'svg'
    #FigFileN= 'Sorbas_chi.svg'
    #FigFileName= DataDirectory+FigFileN
    elevation_threshold = 1
    
    #======================================
    # Uncomment this for a basic plot of the hillshade draped over the elevation, 
    # With chi on top
    #LSDP.BasicChiPlotGridPlot(DEMname,HSname,ChiName, 'gray','gray',
    #                        '$k_{sn}$',(0,0),
    #                        0.4,FigName,'png',elevation_threshold)  
    #======================================
    
    #======================================
    # Uncomment this for a plot of channels color coded by their sources
    #FigName2 = DataDirectory+'TestChannelMap.png'    
    #LSDP.BasicChannelPlotGridPlotCategories(DEMname,HSname,ChiName, 'gray','gray',
    #                        '$Channel$',(0,0),
    #                        0.4,FigName2,'png',elevation_threshold,'source_key')      
    #======================================
    
    #======================================
    # Uncomment this for a very rudimentary plot of the chi profiles
    #FigName3 =  DataDirectory+'ChiProfiles.png'    
    #LSDP.ChiProfiles(ChiName, FigName3,'png',elevation_threshold) 
    #======================================
    
    #======================================
    # Uncomment this for a stack of chi plots
    #FigName4 =  DataDirectory+'ChiStackProfiles.png'  
    #first_basin = 0
    #last_basin = 10
    #LSDP.StackedChiProfiles(ChiName, FigName4,'png',elevation_threshold,first_basin,last_basin)  
    #======================================
    
    
    #======================================
    # Uncomment this for just a plot of the hillshade draped over topography
    #FigName = DataDirectory+'Basic_plot.png'
    #LSDP.BasicDrapedPlotGridPlot(DEMname,HSname, 'gray','gray',
    #                        'Elevation in meters',(0,0),
    #                        0.4,FigName,'png') 
    #======================================

    #======================================
    # Uncomment this for a plot of the basins draped over a fancy hillshde map
    #FigName8 = DataDirectory+'BasinPlot.png'
    #LSDP.DrapedOverFancyHillshade(DEMname,HSname,Basinname, 'gray','cubehelix',
    #                        'Basin Number',(0,0),
    #                        0.4,FigName8,'png',elevation_threshold)  
    #======================================    
    
    
    #======================================    
    # Uncomment this for a plot of the basins draped over a fancy hillshde map
    # with the basins annotated onto the figure
    FigName9 = DataDirectory+'NumberedBasinPlot.png'
    LSDP.BasinsOverFancyHillshade(DEMname,HSname,Basinname, BasinInfoName, 'gray','cubehelix',
                            (0,0), 0.4 ,FigName9,'png',elevation_threshold)  
    #======================================    
    
    
    
    
    
    #======================================
    # Uncomment this for a stack of gradient profiles 
    #FigName5 =  DataDirectory+'ChiGradientProfiles.png'
    #FigName6 =  DataDirectory+'FDGradientProfiles.png'
    #first_basin = 0
    #last_basin = 5
    #LSDP.StackedProfilesGradient(ChiName,FigName5,'png',elevation_threshold,first_basin,last_basin,plt.cm.afmhot,'chi',10,'log')  
    #LSDP.StackedProfilesGradient(ChiName,FigName6,'png',elevation_threshold,first_basin,last_basin,plt.cm.afmhot,'flow_distance',100000,'log')  
    #======================================
    
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

    