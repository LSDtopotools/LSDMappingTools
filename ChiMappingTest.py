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
    #DataDirectory = "C:\\Vagrantboxes\\LSDTopoTools\\Topographic_projects\\Meghalaya\\"
    DataDirectory = "T:\\analysis_for_papers\\Meghalaya/chi_analysis\\"
    Base_file = "Mega_clip"
    
    
    #DataDirectory = "/home/smudd/LSDTopoData/India/Southern_india/"
    #Base_file = "SIndia_clip"
    
        
    bil = ".bil"
    
    #Filename = "Mega_clip.bil"
    #HSFilename = "Mega_clip_hs.bil"
    #BasinFilename = "Mega_clip_AllBasins.bil"
    
    DEMname = DataDirectory+Base_file+bil
    HSname = DataDirectory+Base_file+"_hs"+bil
    Basinname = DataDirectory+Base_file+"_AllBasins"+bil
    
    FigName = DataDirectory+'TestChiFull2.png'
    ChiName = DataDirectory+Base_file+'_MChiSegmented.csv'
    BasinInfoName = DataDirectory+Base_file+'_AllBasinsInfo.csv'
    
    #FigFormat = 'svg'
    #FigFileN= 'Sorbas_chi.svg'
    #FigFileName= DataDirectory+FigFileN
    elevation_threshold = 1
    size_format = "geomorphology"

    #======================================
    # This tests the basin sorting
    #Junction_list = [[3,0,2,1],[4],[6,7,8,9,10]]
    Junction_list = [[3,0,1,2],[4,5,6],[7,10,9,12,8,11,13]]
    
    
    #LSDP.BasinKeyToJunction(Junction_list,BasinInfoName)
    #LSDP.BasinOrderer(BasinInfoName, DEMname, "outlet_longitude",reverse=True) 
    threshold_length = 5
    #thisPointData = LSDP.LSDMap_PointData(ChiName)
    #these_source_nodes = LSDP.FindSourceInformation(thisPointData)
    #remaining_sources = LSDP.FindShortSourceChannels(these_source_nodes,threshold_length)
    #print("The remaining sources are: ")
    #print(remaining_sources)
    #print("The number of remaining sources are: "+str(len(remaining_sources)))
    #======================================
    
    
    #======================================
    # Uncomment this for a basic plot of the hillshade draped over the elevation, 
    # With chi on top, using the cubehelix colour scheme
    #FigName = DataDirectory+'Meghalaya_Ksn_plot_CubeHelix.png'
    #LSDP.BasicChiPlotGridPlot(DEMname,HSname,ChiName, 'gray','gray',
    #                        '$k_{sn}$',(0,0),
    #                        0.4,FigName,'png',elevation_threshold,size_format)  
    #======================================

    #======================================
    # Uncomment this for a basic plot of the hillshade draped over the elevation, 
    # With chi on top, using the kirby and whipple colour scheme
    #FigName = DataDirectory+'Chi_plot_KW.png'
    #LSDP.BasicChiPlotGridPlotKirby(DEMname,HSname,ChiName, 'gray','gray',
    #                        '$k_{sn}$',(0,0),
    #                        0.4,FigName,'png',elevation_threshold)  
    #======================================

    
    #======================================
    # Uncomment this for a plot of channels color coded by their sources
    FigName2 = DataDirectory+'Meghalaya_ChannelMap.png'    
    source_thinning_threshold = 10
    LSDP.BasicChannelPlotGridPlotCategories(DEMname,HSname,ChiName, 'gray','gray',
                            '$Channel$',(0,0),
                            0.4,FigName2,'png',elevation_threshold,'source_key',
                            source_thinning_threshold,
                            size_format)      
    #======================================
    
    #======================================
    # Uncomment this for a very rudimentary plot of the chi profiles
    #FigName3 =  DataDirectory+'ChiProfiles.png'
    #this_basins_list = [3,0,2,1]
    #basin_rename_list = []
    #basin_rename_list =[1,3,2,0,4,5,8,6,7,9,10]  
    
    #source_thinning_threshold = 7
    #label_sources = False
    #this_basins_list = [7,10]    
    #basin_rename_order = [7,10,9,12,8,11,13,4,2,1,0,3,5,6]
    #basin_rename_list = LSDP.BasinOrderToBasinRenameList(basin_rename_order)

    
    #LSDP.ChiProfiles(ChiName, FigName3,'png',this_basins_list,
    #                        basin_rename_list,False,elevation_threshold,source_thinning_threshold,
    #                        size_format) 
    #======================================
    
    #======================================
    # Uncomment this for a stack of chi plots
    FigName4 =  DataDirectory+'Meghalaya_SWbasins_chi.png'  
    first_basin = 0
    last_basin = 10
    chi_offset = 15
    this_basins_list = [7,10]    
    basin_rename_order = [7,10,9,12,8,11,13,4,2,1,0,3,5,6]
    basin_rename_list = LSDP.BasinOrderToBasinRenameList(basin_rename_order)
    LSDP.StackedChiProfiles(ChiName, FigName4,'png',elevation_threshold,
                            first_basin,last_basin,this_basins_list,
                            basin_rename_list,chi_offset,False,
                            source_thinning_threshold,
                            size_format)  
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
    # The basins will have no numbering
    #FigName8 = DataDirectory+'BasinPlot.png' 
    #LSDP.DrapedOverFancyHillshade(DEMname,HSname,Basinname, 'gray','cubehelix',
    #                        'Basin Number',(0,0),
    #                        0.4,FigName8,'png',elevation_threshold)  
    #======================================    
    
    
    #======================================    
    # Uncomment this for a plot of the basins draped over a fancy hillshde map
    # with the basins annotated onto the figure
    #FigName9 = DataDirectory+'MeghalayaBasins.png'
    #spread = 15
    #basin_rename_order = [7,10,9,12,8,11,13,4,2,1,0,3,5,6]
    #basin_rename_list = LSDP.BasinOrderToBasinRenameList(basin_rename_order)
    
    #print basin_rename_list
    
    #LSDP.BasinsOverFancyHillshade(DEMname,HSname,Basinname, BasinInfoName, 'gray','cubehelix',
    #                        (0,0), 0.4 ,FigName9,'png',
    #                        elevation_threshold,
    #                        Junction_list, basin_rename_list,spread,
    #                        ChiName,
    #                        False, threshold_length,size_format)  
    #======================================    

    
    
    
    #======================================
    # Uncomment this for a stack of gradient profiles 
    #FigName5 =  DataDirectory+'Meghalaya_ChiGradientProfilesNorth.png'
    #FigName6 =  DataDirectory+'Meghalaya_FDGradientProfilesNorth.png'
    
    #FigName55 =  DataDirectory+'Meghalaya_ChiGradientProfilesSouth.png'
    #FigName66 =  DataDirectory+'Meghalaya_FDGradientProfilesSouth.png'

    #FigName555 =  DataDirectory+'Meghalaya_ChiGradientProfilesSouth_Funny.png'
    
    #first_basin = 0
    #last_basin = 6
    
    #source_thinning_threshold = 4
    #label_sources = False
    #this_basins_list = [3,0,1,2]
    #basin_rename_list = [1,3,2,0,4,5,8,6,7,9,10]
    
    #basin_rename_order = [7,10,9,12,8,11,13,4,2,1,0,3,5,6]
    #basin_rename_list = LSDP.BasinOrderToBasinRenameList(basin_rename_order)
    
    #LSDP.StackedProfilesGradient(ChiName,FigName5,'png',elevation_threshold,first_basin,last_basin,this_basins_list,basin_rename_list,plt.cm.afmhot,'chi',15,'log',label_sources,source_thinning_threshold,size_format)  
    #LSDP.StackedProfilesGradient(ChiName,FigName6,'png',elevation_threshold,first_basin,last_basin,this_basins_list,basin_rename_list,plt.cm.afmhot,'flow_distance',100000,'log',label_sources,source_thinning_threshold,size_format)

    #this_basins_list = [10,9,12,8,11,13]
    #LSDP.StackedProfilesGradient(ChiName,FigName55,'png',elevation_threshold,first_basin,last_basin,this_basins_list,basin_rename_list,plt.cm.afmhot,'chi',10,'log',label_sources,source_thinning_threshold,size_format)  
    #LSDP.StackedProfilesGradient(ChiName,FigName66,'png',elevation_threshold,first_basin,last_basin,this_basins_list,basin_rename_list,plt.cm.afmhot,'flow_distance',50000,'log',label_sources,source_thinning_threshold,size_format)  

    #this_basins_list = [7,10]
    #LSDP.StackedProfilesGradient(ChiName,FigName555,'png',elevation_threshold,first_basin,last_basin,this_basins_list,basin_rename_list,plt.cm.afmhot,'chi',10,'log',label_sources,source_thinning_threshold,size_format)  
    #LSDP.StackedProfilesGradient(ChiName,FigName66,'png',elevation_threshold,first_basin,last_basin,this_basins_list,basin_rename_list,plt.cm.afmhot,'flow_distance',50000,'log',label_sources,source_thinning_threshold,size_format)  

    
    
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

    