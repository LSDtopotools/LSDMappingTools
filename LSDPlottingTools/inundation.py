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
import glob
import re


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
    return area
    
def floodplain_mean_depth(water_raster, floodplain_mask, threshold=0.0):
    """Calculates the mean waterdepth on the floodplain"""
    
    # I.e. mask where we are NOT in a floodplain (flagged==1)
    floodplain_waters = _np.ma.masked_where(floodplain_mask != 1, water_raster)
    #plt.imshow(floodplain_waters)
    
    mean_water_depth_on_floodplain = _np.mean(floodplain_waters)
    print("Mean water depth in floodplain: ", mean_water_depth_on_floodplain)
    #print(floodplain_waters)
    return mean_water_depth_on_floodplain
    
def main_channel_mean_depth(water_raster, floodplain_mask, stream_mask, threshold=0.0):
    """Calculates the mean water depth in the floodplain channel.
       I.e. does not include channel headwaters outside the floodplain."""
    #floodplain_channel_waters = _np.ma.masked_where(_np.logical_and(floodplain_mask !=1, stream_mask >= 0), water_raster)   
    #floodplain_channel_waters = water_raster[_np.logical_and(floodplain_mask==1 , (~_np.isnan(stream_mask)))]
    #floodplain_channel_waters = water_raster[(floodplain_mask==1) & (~_np.isnan(stream_mask))]
    
    # Floodplain minus channel...
    #floodplain_channel_waters = _np.ma.masked_where( (floodplain_mask !=1 & _np.isnan(stream_mask) ), water_raster)
    
    # Channels within the flood plain
    #floodplain_channel_waters = _np.ma.masked_where( (floodplain_mask != 1 & ~_np.isnan(stream_mask) ), water_raster)
    
    # Get fourth order streams
    floodplain_channel_waters = _np.ma.masked_where( stream_mask!=4 , water_raster)
    
    mean_channel_depth = _np.mean(floodplain_channel_waters)
    print("Mean main channel depth: ", mean_channel_depth)

    #plt.imshow(floodplain_mask)
    #plt.imshow(stream_mask)
    plt.imshow(floodplain_channel_waters)
    
    return mean_channel_depth
    
def timestep_string_from_filename():
    pass

def natural_key(string_):
    """See http://www.codinghorror.com/blog/archives/001018.html"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

def simulation_inundation_timeseries(glob_wildcard, floodplain_mask, stream_mask, threshold=0):
    """Creates a timeseries of a given inundation metric. 
    
    Options should be:
        Inundation Area (Entire catchment)
        Mean Water Depth (Entire catchment)
        Mean Water Depth (Floodplain only)
        Mean Water Depth (Channel)
    """
    # Create an empty array with the correct number of columns
    data_array = _np.empty((0,5), int)
    print("Data array shape: ", data_array.shape)

    for water_raster_file in sorted(glob.glob(glob_wildcard)):
        print(water_raster_file)
        water_raster = lsdgdal.ReadRasterArrayBlocks(water_raster_file)
        
        timestep_row = [] # Empty list to store elements of the row
        cur_timestep = str("TEST") # get the current timestep by parsing the filename

        # Inundation area
        this_inundation_area = calculate_waterinundation_area(water_raster, DX, 0.02)
        this_mean_catchment_waterdepth = calculate_mean_waterdepth(water_raster)
        this_mean_floodplain_waterdepth = floodplain_mean_depth(water_raster, floodplain_mask)
        this_mean_mainchannel_waterdepth = main_channel_mean_depth(water_raster, floodplain_mask, stream_mask)
        
        # Append each value to the current row list object
        timestep_row.append([cur_timestep,
                             this_inundation_area,
                             this_mean_catchment_waterdepth,
                             this_mean_floodplain_waterdepth,
                             this_mean_mainchannel_waterdepth])
        # Convert the list into a numpy array
        timestep_row = _np.asarray(timestep_row)
        
        # Now append (stack) that row onto the bottom of the array (axis=0)
        data_array = _np.append(data_array, timestep_row, axis=0)
    
    print(data_array)
    print(data_array.shape)
    



    
"""Get your rasters into arrays"""    
water_raster_wildcard = "/run/media/dav/SHETLAND/ModelRuns/Boscastle_storms/Lumped/Hydro/WaterDepths*.asc"
water_raster_file = "/mnt/SCRATCH/Analyses/HydrogeomorphPaper/peak_flood_maps/boscastle/peak_flood/WaterDepths2400_GRID_HYDRO.asc"
#raster_file = "/run/media/dav/SHETLAND/Analyses/HydrogeomorphPaper/peak_flood_maps/boscastle/peak_flood/WaterDepths2400_GRID_HYDRO.asc"
floodplain_file = "/mnt/SCRATCH/Analyses/ChannelMaskAnalysis/floodplain_boscastle/BoscastleElevations_FP.bil"
stream_raster_file = "/mnt/SCRATCH/Analyses/ChannelMaskAnalysis/floodplain_boscastle/BoscastleElevations_SO.bil"

water_raster = lsdgdal.ReadRasterArrayBlocks(water_raster_file)
floodplain_mask = lsdgdal.ReadRasterArrayBlocks(floodplain_file)
stream_mask = lsdgdal.ReadRasterArrayBlocks(stream_raster_file)
#print(stream_mask)

DX = lsdgdal.GetUTMMaxMin(water_raster_file)[0]   # I never realised you could do this!
print(DX)

"""Calculate the depths and areas"""
#calculate_mean_waterdepth(water_raster)
#calcualte_max_waterdepth(water_raster)
#calculate_waterinundation_area(water_raster, DX, 0.02)
#floodplain_mean_depth(water_raster, floodplain_mask)
#main_channel_mean_depth(water_raster, floodplain_mask, stream_mask)

"""Make the timeseries file"""
simulation_inundation_timeseries(water_raster_wildcard, floodplain_mask, stream_mask)
