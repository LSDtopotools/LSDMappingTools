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
    #DataDirectory = "T://analysis_for_papers//Cosmo_paper//Tibet//for_plotting//"
    #DataDirectory = "C://basin_data//CosmoPaper//DEMs//"
    DataDirectory = "M://papers//Mudd_SOS//Betics//"
    #Filename = "SanBern.bil"
    Filename = "betics_chi_TA3000_clip_UTM30.bil"
    #Filename = "SpawnedBasin_07C13-(Q8)_SH.bil"
    #Filename = "CRNvariable_long_0_0_1var128.asc"
    #DrapeFileName = "CRNvariable_long_0_0_1var128_erosion.asc"
    DrapeFileName = "SpawnedBasin_07C13-(Q8)_BASINS.bil"
    DEMName = "SpawnedBasin_07C13-(Q8).bil"
    ThisFile = DataDirectory+Filename
    #DrapeFile =DataDirectory+DrapeFileName 
    #DEMFile = DataDirectory+DEMName
    #Shield2 = DataDirectory+"SpawnedBasin_Q8_phi30.bil"
    
    FigFormat = 'svg'
    FigFileN= 'Sorbas_chi.svg'
    FigFileName= DataDirectory+FigFileN
    
    #ShieldFigName = DataDirectory+'Shielding.svg'

    
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
    
    
    tcmap = 'jet'
    drape_cmap = 'gray'
    shield_map = 'summer'
    drape_alpha = 0.4
    tcmapcolorbarlabel='Chi'
    clim_val = (50,300)
    auto_clim = (1000,4500)
    #LSDP.BasicDensityPlot(ThisFile,tcmap,tcmapcolorbarlabel,clim_val)
    #LSDP.DrapedPlot(ThisFile,DrapeFile)
    
    #LSDP.BasicDensityPlotGridPlot(ThisFile,tcmap,tcmapcolorbarlabel)

    #Plot the basin over the elevation    
    cmap_label = "Elevation (m)"
    #LSDP.BasicDrapedPlotGridPlot(DEMFile,DrapeFile,tcmap,drape_cmap, cmap_label,
    #                             auto_clim,drape_alpha,FigFileName,FigFormat)

    # Now plot the two shielding rasters
    LSDP.BasicDensityPlotGridPlot(ThisFile,shield_map,tcmapcolorbarlabel,clim_val,
                                  FigFileName,FigFormat)

    #FigFormat = 'svg'
    #SNFileN= 'Shield30.svg'
    #ShieldFigName2= DataDirectory+SNFileN
    # Now plot the two shielding rasters
    #LSDP.BasicDensityPlotGridPlot( Shield2,shield_map,tcmapcolorbarlabel,clim_val,
    #                              ShieldFigName2,FigFormat)
                                  

if __name__ == "__main__":
    #fit_weibull_from_file(sys.argv[1]) 
    TestNewMappingTools() 