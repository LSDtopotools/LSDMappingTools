#!/usr/bin/env python2
# -*- coding: utf-8 -*-


"""
Raster/DEM differencer

For when you've fucked up your C++ code and have to
do the DEM-differencing by hand in Python...


Created on Sun Jan 15 16:18:35 2017

@author: dav

"""
import LSDPlottingTools as LSDP

DataDirectory = "/mnt/SCRATCH/Analyses/HydrogeomorphPaper/peak_flood_maps/boscastle/erode_diff/test_raster_diff_func/"


LSDP.RasterDifference(DataDirectory + "BoscastleElevations0.asc", 
                      DataDirectory + "Elevations4200.asc", 
                      raster_band=1, 
                      OutFileName="BoscastleElevDiff_GRID_DLIM.bil", 
                      OutFileType="ENVI")