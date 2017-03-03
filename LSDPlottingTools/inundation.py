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
import matplotlib.pyplot as plt


def calculate_mean_waterdepth(raster):
    print("Mean water depth (whole raster): ", _np.mean(raster))
    return _np.mean(raster)
    
def calcualte_max_waterdepth(raster):
    #print(_np.max(raster))
    return _np.max(raster)
    
def calculate_waterinundation_area(raster, cellsize, threshold):
    """Note: this will probably need some sort of threshold as Caesar maps 
    out very small water depths and so could give huge 'inundation' areas."""
    total_cells = _np.count_nonzero(raster > threshold)
    area = DX * DX * total_cells  # metres
    
    print("Inundation area is: ", area, " metres square")
    
def floodplain_mean_depth(water_raster, floodplain_mask, threshold=0.0):
    """Calculates the mean waterdepth on the floodplain"""
    
    # I.e. mask where we are NOT in a floodplain (flagged==1)
    floodplain_waters = _np.ma.masked_where(floodplain_mask != 1, water_raster)
    plt.imshow(floodplain_waters)
    
    mean_water_depth_on_floodplain = _np.mean(floodplain_waters)
    print("Mean water depth in floodplain: ", mean_water_depth_on_floodplain)
    print(floodplain_waters)

    
    
raster_file = "/mnt/SCRATCH/Analyses/HydrogeomorphPaper/peak_flood_maps/boscastle/peak_flood/WaterDepths2400_GRID_HYDRO.asc"
#raster_file = "/run/media/dav/SHETLAND/Analyses/HydrogeomorphPaper/peak_flood_maps/boscastle/peak_flood/WaterDepths2400_GRID_HYDRO.asc"
floodplain_file = "/mnt/SCRATCH/Analyses/ChannelMaskAnalysis/floodplain_boscastle/BoscastleElevations_FP.bil"

water_raster = lsdgdal.ReadRasterArrayBlocks(raster_file)
floodplain_mask = lsdgdal.ReadRasterArrayBlocks(floodplain_file)

DX = lsdgdal.GetUTMMaxMin(raster_file)[0]   # I never realised you could do this!
print(DX)

calculate_mean_waterdepth(water_raster)
calcualte_max_waterdepth(water_raster)
calculate_waterinundation_area(water_raster, DX, 0.02)

floodplain_mean_depth(water_raster, floodplain_mask)