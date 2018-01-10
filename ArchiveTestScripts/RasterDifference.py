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

DataDirectory = "/run/media/dav/SHETLAND/Analyses/HydrogeomorphPaper/erode_diff/Difference_UNIFORM_GRIDDED/"


LSDP.RasterDifference(DataDirectory + "RyedaleElevDiff_GRIDDED_TLIM.bil", 
                      DataDirectory + "RyedaleElevDiff_UNIFORM_TLIM.bil", 
                      raster_band=1, 
                      OutFileName="RyedaleErodeDiff_GRID_UNI_TLIMM.bil", 
                      OutFileType="ENVI")

