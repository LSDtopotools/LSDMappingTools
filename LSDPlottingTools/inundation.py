# -*- coding: utf-8 -*-
"""
inundation.py

A script that calculates the area of the model/DEM domain that is inundated by
water; either across the entire domain, the floodplain (using Fiona's floodplain
ID algorithm, or in the Channel, using the channel network and measuring along 
this line. The mean, max, and total area wil be able to be calculated.)
"""

import LSDMap_GDALIO as lsdgdal
import numpy as _np


def calculate_mean_waterdepth(raster):
    return _np.mean(raster)
    
def calcualte_max_waterdepth(raster):
    #print(_np.max(raster))
    return _np.max(raster)
    
def calculate_waterinundation_area(raster, cellsize, threshold):
    """Note: this will probably need some sort of threshold as Caesar maps 
    out very small water depths and so could give huge 'inundation' areas."""
    total_cells = _np.count_nonzero(raster > threshold)
    area = DX * DX * total_cells  # metres
    
    print("Area is: ", area, " metres square")
    

raster_file = "/run/media/dav/SHETLAND/Analyses/HydrogeomorphPaper/peak_flood_maps/boscastle/peak_flood/WaterDepths2400_GRID_HYDRO.asc"
water_raster = lsdgdal.ReadRasterArrayBlocks(raster_file)

DX = lsdgdal.GetUTMMaxMin(raster_file)[0]   # I never realised you could do this!
print(DX)

calculate_mean_waterdepth(water_raster)
calcualte_max_waterdepth(water_raster)
calculate_waterinundation_area(water_raster, DX, 0.02)