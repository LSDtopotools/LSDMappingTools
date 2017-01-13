#==============================================================================
# These are some scripts for testing the functionality of LSDMappingTools
#==============================================================================
# -*- coding: utf-8 -*-
"""
Created on Tue May 05 14:08:16 2015

@author: smudd
"""
import glob as glob
import numpy as np
import LSDPlottingTools as LSDP
import LSDPlottingTools.LSDMatplotlibExtensions as mplext

import matplotlib.pyplot as plt



#fit_weibull_from_file(sys.argv[1]) 
#TestNewMappingTools2() 
#ResetErosionRaster()
#FloodThenHillshade()
#FixStupidNoData()

def MultiDrapeMaps(DataDir, ElevationRaster, DrapeRasterWild, cmap):
    """
    Plots flood extents from water depth rasters
    draped over the catchment elevation raster
    in a series of subplots
    
    Takes a wildcard for the drapes
    Expexts a fixed elevation raster, but this could be
    modified in future.
    """
    f, ax_arr = plt.subplots(2, 2, figsize=(10, 5), sharex=True, sharey=True)
    ax_arr = ax_arr.ravel()
    
    FPFiles = sorted(glob.glob(DataDirectory+DrapeRasterWild), key=str)
    n_files = len(FPFiles)
    print "Number of files = ", n_files
    
    elev_raster_file = DataDir + ElevationRaster
    
    hillshade = LSDP.Hillshade(elev_raster_file)
    #hillshade_array = LSDP.ReadRasterArrayBlocks(elev_raster_file)
    
    # now get the extent
    extent_raster = LSDP.GetRasterExtent(elev_raster_file)
    
    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # now get the tick marks    
    n_target_tics = 5
    xlocs,ylocs,new_x_labels,new_y_labels = LSDP.GetTicksForUTM(elev_raster_file,x_max,x_min,y_max,y_min,n_target_tics)  

    print "xmax: " + str(x_max)
    print "xmin: " + str(x_min)
    print "ymax: " + str(y_max)
    print "ymin: " + str(y_min)
    
    for i in range(n_files):
        
        print "The floodplain file name is: ", FPFiles[i]
        FP_raster = LSDP.ReadRasterArrayBlocks(FPFiles[i])
        #FP_raster = np.ma.masked_where(FP_raster <= 0, FP_raster)
        
        low_values_index = FP_raster < 0.0005
        FP_raster[low_values_index] = np.nan
        
        ax_arr[i].imshow(hillshade, "gray", extent=extent_raster, interpolation="nearest")
        ax_arr[i].imshow(FP_raster, cmap, extent=extent_raster, alpha=1.0, interpolation="none")
        
        
        
# truncate the colourmap    
cmap = plt.get_cmap("Blues")
trunc_cmap = mplext.truncate_colormap("Blues", 0.4, 1.0)

#DataDirectory = "/run/media/dav/SHETLAND/Analyses/Ryedale_storms_simulation/Gridded/DetachLim/"
DataDirectory = "/mnt/SCRATCH/Dev/LSDMappingTools/test_rasters/"
filename = DataDirectory + "Elevations0.asc"
drapename = DataDirectory + "WaterDepths2880.asc"

MultiDrapeMaps(DataDirectory, "Elevations0.asc", "WaterDepths*.asc", trunc_cmap)

# Create the drape array from one of the Catchment model output rasters
#drape_array = LSDP.ReadRasterArrayBlocks(drapename)
# Optional: A lot of the output rasters contain very small values for certain 
# things like water depth or elevation difference, so you can mask this below:
#low_values_index = drape_array < 0.005
#drape_array[low_values_index] = np.nan


"""
LSDP.DrapedOverHillshade(filename,drape_array,clim_val=(0,400), \
                         drape_cmap=trunc_cmap, colorbarlabel='Elevation in meters',\
                         ShowColorbar=True, ShowDrapeColorbar=True,
                         drape_cbarlabel = "Water depth (m)",
                         drape_alpha=1.0)
"""