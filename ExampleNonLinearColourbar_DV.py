#==============================================================================
# These are some scripts for testing the functionality of LSDMappingTools
#==============================================================================
# -*- coding: utf-8 -*-
"""
Created on Tue May 05 14:08:16 2015

@author: dav
"""
import glob as glob
import os.path
import numpy as np
import LSDPlottingTools as LSDP
import matplotlib.cm as cm
#import lsdmatplotlibextensions as mplext

import matplotlib.pyplot as plt
from matplotlib import ticker

# Get favourite plotting fonts and sizes
LSDP.init_plotting_DV()

# Truncate the colour map
#trunc_cmap = LSDP.colours.truncate_colormap("Blues", 0.4, 1.0)

# Option for getting a discrete colour map
#discreet_cmap = mplext.colours.discrete_colourmap(8, "Blues")
#discreet_cmap = mplext.colours.cmap_discretize(8, trunc_cmap)

# Non-linear colourmap
levels = [-2.5, -1.5 -0.5, -0.25, 0, 0.25, 0.5, 1.5, 2.5 ]
#levels = levels[levels <= tmax]
print levels
levels.sort()
print levels

cmap_lin = cm.jet

nonlincmap = LSDP.colours.nonlinear_colourmap(cmap_lin, levels)
print type(nonlincmap)

#DataDirectory = "/run/media/dav/SHETLAND/Analyses/Ryedale_storms_simulation/Gridded/DetachLim/"
#DataDirectory = "/mnt/SCRATCH/Analyses/HydrogeomorphPaper/peak_flood_maps/boscastle/erode_diff/"
#DataDirectory = "/mnt/SCRATCH/Analyses/HydrogeomorphPaper/peak_flood_maps/boscastle/erode_diff/test_raster_diff_func/"
DataDirectory = "/run/media/dav/SHETLAND/Analyses/HydrogeomorphPaper/erode_diff/test_raster_diff_func/"

LSDP.MultiDrapeErodeDiffMaps(DataDirectory, "BoscastleElevations0.asc", "Boscastle*.bil", 
               cmap=nonlincmap, 
               drape_min_threshold= -2.5,
               drape_max_threshold= 2.5,
               cbar_label = "DEM difference (m)",
               middle_mask_range = (-0.01, 0.01) 
               )
