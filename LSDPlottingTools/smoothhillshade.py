#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 11:40:41 2017

@author: dav
"""

import LSDPlottingTools.LSDMap_GDALIO as lsdio

import numpy as np
import scipy as sp


Zenith = 45
Azimuth = 315
ZFactor = 1

Directory = "/mnt/SCRATCH/Analyses/HydrogeomorphPaper/erode_diff/test_raster_diff_func/"
BackgroundRasterName = "BoscastleElevations0.asc"

File = Directory + BackgroundRasterName

RasterData = lsdio.ReadRasterArrayBlocks(File)

def Hillshade_Smooth(RasterData, altitude, azimuth, z_factor):
    """Plots a Hillshade a la LSDRaster"""
    
    zenith_rad = sp.deg2rad(altitude)
    azimuth_rad = sp.deg2rad(azimuth)
    
    
    