#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 11:47:22 2017

@author: dav
"""

import lsdmapartist as lsdmap

#Directory = "/mnt/SCRATCH/Dev/LSDMappingTools/test_rasters/peak_flow_rasters/"
Directory = "/run/media/dav/SHETLAND/Analyses/HydrogeomorphPaper/peak_flood_maps/boscastle/peak_flood/"
BackgroundRasterName = "Elevations0.asc"
DrapeRasterName = "WaterDepths2400_GRID_HYDRO.asc"

#raster = BaseRaster(RasterName, DataDirectory)
dp = lsdmap.DrapePlot(DrapeRasterName, BackgroundRasterName, Directory, "Blues", drape_min_threshold=0.05)