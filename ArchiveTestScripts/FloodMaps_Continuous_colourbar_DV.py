#==============================================================================
# These are some scripts for testing the functionality of LSDMappingTools
#==============================================================================
# -*- coding: utf-8 -*-
"""
Created on Tue May 05 14:08:16 2015

@author: smudd
"""
from __future__ import print_function
import glob as glob
import os.path
import numpy as np
import LSDPlottingTools as LSDP
import lsdmatplotlibextensions as mplext

import matplotlib.pyplot as plt
from matplotlib import ticker



def findmaxval_multirasters(FileList):
    """
    Loops through a list or array of rasters (np arrays)
    and finds the maximum single value in the set of arrays.
    """
    overall_max_val = 0
    
    for i in range (len(FileList)):
        
        raster_as_array = LSDP.ReadRasterArrayBlocks(FileList[i])
        this_max_val = np.max(raster_as_array)
        
        if this_max_val > overall_max_val:
            overall_max_val = this_max_val
            print(overall_max_val)
            
    return overall_max_val

def MultiDrapeMaps(DataDir, ElevationRaster, DrapeRasterWild, cmap, drape_min_threshold=None, drape_max=None):
    """
    Plots flood extents from water depth rasters
    draped over the catchment elevation raster
    in a series of subplots
    
    Takes a wildcard for the drapes
    Expexts a fixed elevation raster, but this could be
    modified in future.
    
    Thought: consider, if plotting multiple datasets, how you
    are going to deal with min a max values in the colur range.
    imshow will automatically set vmin and vmax and stretch the colour bar 
    over this - which can be visually misleading. Ideally, you
    want to have the same colour map used for *all* subplots, and 
    this is not default behaviour.
    """
    f, ax_arr = plt.subplots(2, 2, figsize=(10, 5), sharex=True, sharey=True)
    ax_arr = ax_arr.ravel()
    
    FPFiles = sorted(glob.glob(DataDirectory+DrapeRasterWild), key=str)
    n_files = len(FPFiles)
    print("Number of files = ", n_files)
    
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

    print("xmax: " + str(x_max))
    print("xmin: " + str(x_min))
    print("ymax: " + str(y_max))
    print("ymin: " + str(y_min))
    
    """
    Find the maximum water depth in all rasters.
    You need this to normalize the colourscale accross
    all plots when teh imshow is done later.
    """

    try:
        print("Calculating max drape raster value by scanning rasters...")
        max_water_depth = findmaxval_multirasters(FPFiles)
        drape_max = max_water_depth
        
    except:
        print("Something went wrong trying to obtain the max value in \
                your drape raster file list.")
    finally:
        print("The drape(s) max value is set to: ", drape_max) 
    
    
    for i in range(n_files):
        
        print("The floodplain file name is: ", FPFiles[i])
        FP_raster = LSDP.ReadRasterArrayBlocks(FPFiles[i])
        #FP_raster = np.ma.masked_where(FP_raster <= 0, FP_raster)
        
        filename = os.path.basename(FPFiles[i])
        title = mplext.labels.make_line_label(filename)
        print(title)
        
        low_values_index = FP_raster < drape_min_threshold
        FP_raster[low_values_index] = np.nan
        
        im = ax_arr[i].imshow(hillshade, "gray", extent=extent_raster, interpolation="nearest")
        """
        Now we can set vmax to be the maximum water depth we calcualted earlier, making our separate
        subplots all have the same colourscale
        """
        im = ax_arr[i].imshow(FP_raster, cmap, extent=extent_raster, 
                                alpha=1.0, interpolation="none", 
                                vmin=drape_min_threshold, 
                                vmax=drape_max)
        ax_arr[i].set_title(title)
 
    f.subplots_adjust(right=0.8)
    cax = f.add_axes([0.9, 0.1, 0.03, 0.8])
    
    cbar = f.colorbar(im, cax=cax) 
    cbar.set_label("Water depth (m)")
    #cbar.set_ticks(np.linspace(0, 8, 8))
    #cbar = mplext.colours.colorbar_index(f, cax, 8, cmap, 
    #                                     drape_min_threshold, drape_max)
    cbar.set_label("Water depth (m)")
    #tick_locator = ticker.MaxNLocator(nbins=8)
    #cbar.locator = tick_locator
    #cbar.update_ticks()
    
    f.text(0.5, 0.04, 'Easting (m)', ha='center')
    f.text(0.04, 0.5, 'Northing (m)', va='center', rotation='vertical')
        
# truncate the colourmap    
#cmap =  plt.get_cmap("Blues")

trunc_cmap = mplext.colours.truncate_colormap("Blues", 0.4, 1.0)

#discreet_cmap = mplext.colours.discrete_colourmap(8, "Blues")
#discreet_cmap = mplext.colours.cmap_discretize(8, trunc_cmap)

#DataDirectory = "/run/media/dav/SHETLAND/Analyses/Ryedale_storms_simulation/Gridded/DetachLim/"
DataDirectory = "/mnt/SCRATCH/Dev/LSDMappingTools/test_rasters/peak_flow_rasters/"
filename = DataDirectory + "Elevations0.asc"
drapename = DataDirectory + "WaterDepths2880.asc"

MultiDrapeMaps(DataDirectory, "Elevations0.asc", "WaterDepths*.asc", 
               trunc_cmap, 
               drape_min_threshold=0.01)

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