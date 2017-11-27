#=============================================================================
# These functions create figures for Basic visualization
# 
#
# It creates separate plots for each basin in the DEM.
#
# Authors:
#     Simon M. Mudd
#     Fiona J. Clubb
#=============================================================================
#=============================================================================
# IMPORT MODULES
#=============================================================================
# set backend to run on server
import matplotlib
matplotlib.use('Agg')

import numpy as np
import LSDPlottingTools as LSDP
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.ticker as ticker
import pandas as pd
from matplotlib import colors
import math
import os
import subprocess
#from shapely.geometry import Polygon
from LSDMapFigure import PlottingHelpers as Helper
from LSDMapFigure.PlottingRaster import MapFigure
from LSDMapFigure.PlottingRaster import BaseRaster


def PlotTopoRaster(DataDirectory, fname_prefix, size_format='ESURF', FigFormat='png', colors = "terrain"):
    """
    Creates a basic Terrain topographic raster. Needs the Hillshade 

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.

    Returns:
        Shaded relief plot with the basins coloured by basin ID

    Author: BG, FJC
    """
    # check if a directory exists for the chi plots. If not then make it.
    raster_directory = DataDirectory+'raster_plots/'
    if not os.path.isdir(raster_directory):
        os.makedirs(raster_directory)

    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # set figure sizes based on format
    if size_format == "geomorphology":
        fig_width_inches = 6.25
    elif size_format == "big":
        fig_width_inches = 16
    else:
        fig_width_inches = 4.92126

    # going to make the basin plots - need to have bil extensions.
    print("I'm going to make a cabic topographic plot")

    # get the rasters
    raster_ext = '.bil'
    ## Just checking if you have a PP version of it
    if os.path.isfile(DataDirectory + fname_prefix +"_PP.bil"):
        BackgroundRasterName = fname_prefix+"_PP"+raster_ext
    else:
        BackgroundRasterName = fname_prefix+raster_ext

    HillshadeName = fname_prefix+'_hs'+raster_ext
    

    # create the map figure
    MF = MapFigure(BackgroundRasterName, DataDirectory,coord_type="UTM_km", colourbar_location='None')

    # Drape the hillshade and add the color
    ## Frist plot the terrain toporaster 
    MF.add_drape_image(BackgroundRasterName,DataDirectory, # Calling the function will add a drapped raster on the top of the background on
                        colourmap = colors, # colormap used for this raster, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                        alpha=1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                        show_colourbar = True, # Well, this one is explicit I think
                        colorbarlabel = "None",
                        NFF_opti = True) 
    ## Drape the Hillshade raster
    MF.add_drape_image(HillshadeName,DataDirectory, # Calling the function will add a drapped raster on the top of the background on
                        colourmap = "gray", # colormap used for this raster, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                        alpha=0.4, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                        show_colourbar = True, # Well, this one is explicit I think
                        colorbarlabel = "None",
                        NFF_opti = True) 

    # Save the figure
    ImageName = raster_directory+fname_prefix+'_Topo.'+FigFormat
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 300)

def PlotSlopeRaster(DataDirectory, fname_prefix, size_format='ESURF', FigFormat='png'):
    """
    Creates a basic Slope Map with a [0,2] scale

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.

    Returns:
        Shaded relief plot with the basins coloured by basin ID

    Author: BG, FJC
    """
    # check if a directory exists for the chi plots. If not then make it.
    raster_directory = DataDirectory+'raster_plots/'
    if not os.path.isdir(raster_directory):
        os.makedirs(raster_directory)

    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # set figure sizes based on format
    if size_format == "geomorphology":
        fig_width_inches = 6.25
    elif size_format == "big":
        fig_width_inches = 16
    else:
        fig_width_inches = 4.92126

    # going to make the basin plots - need to have bil extensions.
    print("I'm going to make a cabic topographic plot")

    # get the rasters
    raster_ext = '.bil'
    ## Just checking if you have a PP version of it
    BackgroundRasterName = fname_prefix+"_slope"+raster_ext
    
    

    # create the map figure
    MF = MapFigure(BackgroundRasterName, DataDirectory,coord_type="UTM_km", colourbar_location='None')

    # Drape the hillshade and add the color
    ## Frist plot the black background
    MF.add_drape_image(BackgroundRasterName,DataDirectory, # Calling the function will add a drapped raster on the top of the background on
                        colourmap = "gray", # colormap used for this raster, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                        alpha=1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                        show_colourbar = True, # Well, this one is explicit I think
                        colorbarlabel = "None",
                        colour_min_max = [0,100000],
                        custom_min_max = [0,0.1],
                        NFF_opti = True) 
    ## Drape the slope raster 
    MF.add_drape_image(BackgroundRasterName,DataDirectory, # Calling the function will add a drapped raster on the top of the background on
                        colourmap = "viridis", # colormap used for this raster, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                        alpha=1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                        show_colourbar = True, # Well, this one is explicit I think
                        colour_min_max = [0,2],
                        colorbarlabel = "None",
                        NFF_opti = True) 

    # Save the figure
    ImageName = raster_directory+fname_prefix+'_Slopo.'+FigFormat
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 300)

def PlotCurveRaster(DataDirectory, fname_prefix, size_format='ESURF', FigFormat='png'):
    """
    Creates a basic Slope Map with a [0,2] scale

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.

    Returns:
        Shaded relief plot with the basins coloured by basin ID

    Author: BG, FJC
    """
    # check if a directory exists for the chi plots. If not then make it.
    raster_directory = DataDirectory+'raster_plots/'
    if not os.path.isdir(raster_directory):
        os.makedirs(raster_directory)

    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # set figure sizes based on format
    if size_format == "geomorphology":
        fig_width_inches = 6.25
    elif size_format == "big":
        fig_width_inches = 16
    else:
        fig_width_inches = 4.92126

    # going to make the basin plots - need to have bil extensions.
    print("I'm going to make a cabic topographic plot")

    # get the rasters
    raster_ext = '.bil'
    ## Just checking if you have a PP version of it
    BackgroundRasterName = fname_prefix+"_curvature"+raster_ext
    
    

    # create the map figure
    MF = MapFigure(BackgroundRasterName, DataDirectory,coord_type="UTM_km", colourbar_location='None')

    # Drape the hillshade and add the color
    ## Frist plot the black background
    MF.add_drape_image(BackgroundRasterName,DataDirectory, # Calling the function will add a drapped raster on the top of the background on
                        colourmap = "gray", # colormap used for this raster, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                        alpha=1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                        show_colourbar = True, # Well, this one is explicit I think
                        colorbarlabel = "None",
                        colour_min_max = [0,100000],
                        custom_min_max = [0,0.1],
                        NFF_opti = True) 
    ## Drape the slope raster 
    MF.add_drape_image(BackgroundRasterName,DataDirectory, # Calling the function will add a drapped raster on the top of the background on
                        colourmap = "viridis", # colormap used for this raster, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                        alpha=1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                        show_colourbar = True, # Well, this one is explicit I think
                        colour_min_max = [0,2],
                        colorbarlabel = "None",
                        NFF_opti = True) 

    # Save the figure
    ImageName = raster_directory+fname_prefix+'_Curve.'+FigFormat
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 300)


