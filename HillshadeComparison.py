#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 11:47:22 2017

An example of using the LSDMapArtist to create drape plots

@author: dav
"""

import matplotlib.pyplot as plt

from LSDPlottingTools import colours as lsdcolours
from LSDPlottingTools import init_plotting_DV

import LSDPlottingTools.LSDMap_GDALIO as LSDMap_IO
import LSDPlottingTools.LSDMap_BasicPlotting as LSDMap_BP

init_plotting_DV()
#Directory = "/mnt/SCRATCH/Dev/LSDMappingTools/test_rasters/peak_flow_rasters/"
Directory = "/mnt/SCRATCH/Analyses/HydrogeomorphPaper/rainfall_maps/"
#Directory = "/mnt/SCRATCH/Analyses/HydrogeomorphPaper/erode_diff/Difference_UNIFORM_GRIDDED/"
#Directory = "/mnt/SCRATCH/Analyses/HydrogeomorphPaper/erode_diff/test_raster_diff_func/"
BackgroundRasterName = "BoscastleElevations0.asc"
#DrapeRasterName = "BoscastleElevDiff_UNIFORM_TLIM.bil"
DrapeRasterName = "rainfall_totals_boscastle_downscaled.asc"

raster = LSDMap_IO.ReadRasterArrayBlocks(Directory + BackgroundRasterName)

hs = LSDMap_BP.Hillshade(raster)

plt.imshow(hs, cmap="gray")