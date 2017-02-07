#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 11:47:22 2017

An example of using the LSDMapArtist to create drape plots

@author: dav
"""

import lsdmapartist as lsdmap

#Directory = "/mnt/SCRATCH/Dev/LSDMappingTools/test_rasters/peak_flow_rasters/"
Directory = "/run/media/dav/SHETLAND/Analyses/HydrogeomorphPaper/erode_diff/Difference_UNIFORM_GRIDDED/"
BackgroundRasterName = "RyedaleElevations0.asc"
DrapeRasterName = "RyedaleErodeDiff_GRID_UNI_TLIMM.bil"
#DrapeRasterName = "BoscastleErodeDiff_GRID_UNI_TLIMM.bil"
Colourmap = "coolwarm"
drape_min_threshold = None
drape_max_threshold = None
colourbar_label = "Erosion/Deposition (m)"

#raster = BaseRaster(RasterName, DataDirectory)
dp = lsdmap.DrapePlot(DrapeRasterName, BackgroundRasterName, Directory,
                      Colourmap, background_type="Terrain",
                      vmin=-3, vmax=3, middle_mask_range=(-0.02,0.02))

# Customise the DrapePlot
dp.make_drape_colourbar(cbar_label=colourbar_label)
dp.set_coordinate_labels_axis()

dp.show_plot()