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
import LSDPlottingTools.LSDMatplotlibExtensions as mplext


#fit_weibull_from_file(sys.argv[1]) 
#TestNewMappingTools2() 
#ResetErosionRaster()
#FloodThenHillshade()
#FixStupidNoData()

DataDirectory = "/run/media/dav/SHETLAND/Analyses/Ryedale_storms_simulation/Gridded/DetachLim/"
filename = DataDirectory + "Elevations0.asc"
drapename = DataDirectory + "WaterDepths2880.asc"

# Create the drape array from one of the Catchment model output rasters
drape_array = LSDP.ReadRasterArrayBlocks(drapename)
# Optional: A lot of the output rasters contain very small values for certain 
# things like water depth or elevation difference, so you can mask this below:
low_values_index = drape_array < 0.005
drape_array[low_values_index] = np.nan

#cmap = plt.get_cmap("Blues")

trunc_cmap = mplext.truncate_colormap("Blues", 0.4, 1.0)

LSDP.DrapedOverHillshade(filename,drape_array,clim_val=(0,400), \
                         drape_cmap=trunc_cmap, colorbarlabel='Elevation in meters',\
                         ShowColorbar=True, ShowDrapeColorbar=True,
                         drape_cbarlabel = "Water depth (m)",
                         drape_alpha=1.0)
