#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 11:47:22 2017

An example of using the LSDMapArtist to create drape plots

@author: dav
"""

import matplotlib.cm as cm
from LSDMapArtist.drapeplot import DrapeAxes
from LSDPlottingTools import colours as lsdcolours
from LSDPlottingTools import init_plotting_DV

init_plotting_DV()
#Directory = "/mnt/SCRATCH/Dev/LSDMappingTools/test_rasters/peak_flow_rasters/"
Directory = "/mnt/SCRATCH/Analyses/HydrogeomorphPaper/rainfall_maps/"
#Directory = "/mnt/SCRATCH/Analyses/HydrogeomorphPaper/erode_diff/Difference_UNIFORM_GRIDDED/"
#Directory = "/mnt/SCRATCH/Analyses/HydrogeomorphPaper/erode_diff/test_raster_diff_func/"
BackgroundRasterName = "BoscastleElevations0.asc"
#DrapeRasterName = "BoscastleElevDiff_UNIFORM_TLIM.bil"
DrapeRasterName = "rainfall_totals_boscastle_downscaled.asc"

# Standard colourmap
Colourmap = "Blues"

#Non-linear colourmap
##ColourLevels = lsdcolours.nonlinear_colourmap.create_levels(-3.0, 3.0, -0.2, 0.2, -0.5, 0.5)
##Colourmap = lsdcolours.nonlinear_colourmap("seismic", ColourLevels)

# Transformed colourmap
#c = lsdcolours.TransformedColourmap(lambda x: x/2+0.5, cm.jet)

drape_min_threshold = None
drape_max_threshold = None
colourbar_label = "Total accumulated rainfall (mm)"

#raster = BaseRaster(RasterName, DataDirectory)
dp = DrapeAxes(DrapeRasterName, BackgroundRasterName, Directory,
                      Colourmap, background_type="Hillshade", 
                      show_background_colourbar=False,
                      colourbar_label=colourbar_label, drape_alpha=0.7
                      )
                      #colourbar_norm_type="SymLogNorm")


# Customise the DrapePlot
#dp.make_drape_colourbar(cbar_label=colourbar_label)
#dp.set_fig_axis_labels()

dp.show_plot()