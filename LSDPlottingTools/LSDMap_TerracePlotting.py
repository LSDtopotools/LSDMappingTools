#=============================================================================
# These functions create figures for visualising the results of the terrace
# algorithm
#
# Authors:
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
import pandas as pd
from matplotlib import colors
#from shapely.geometry import Polygon
from LSDMapFigure import PlottingHelpers as Helper
from LSDMapFigure.PlottingRaster import MapFigure
from LSDMapFigure.PlottingRaster import BaseRaster
import colours

#--------------------------------------------------------------------------#
# RASTER PLOTS
# These functions interface with the MapFigure object to produce raster plots
# of the terrace locations
#--------------------------------------------------------------------------#

def MakeTerraceIDRasterPlot(DataDirectory, fname_prefix, size_format='ESURF', FigFormat='png'):
    """
    This function makes a shaded relief plot of the DEM with the terraces coloured by their unique ID

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.

    Returns:
        Shaded relief plot with the terrace coloured by ID

    Author: FJC
    """
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

    # going to make the terrace plots - need to have bil extensions.
    print("I'm going to make the terrace raster plot. Your topographic data must be in ENVI bil format or I'll break!!")

    # get the rasters
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext
    TerraceIDName = fname_prefix+'_terrace_IDs'+raster_ext

    # create the map figure
    MF = MapFigure(HillshadeName, DataDirectory,coord_type="UTM_km", colourbar_location='bottom')
    # add the terrace drape
    terrace_cmap = plt.cm.Set2
    MF.add_drape_image(TerraceIDName, DataDirectory, colourmap = terrace_cmap, alpha=0.8)

    ImageName = DataDirectory+fname_prefix+'_terrace_IDs.'+FigFormat
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 300) # Save the figure

def MakeTerraceElevationRasterPlot(DataDirectory, fname_prefix, size_format='ESURF', FigFormat='png'):
    """
    This function makes a shaded relief plot of the DEM with the terraces coloured by their elevation
    compared to the channel

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.

    Returns:
        Shaded relief plot with the terrace coloured by ID

    Author: FJC
    """
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

    # going to make the terrace plots - need to have bil extensions.
    print("I'm going to make the terrace raster plot. Your topographic data must be in ENVI bil format or I'll break!!")

    # get the rasters
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext
    TerraceIDName = fname_prefix+'_terrace_relief_final'+raster_ext

    # create the map figure
    MF = MapFigure(HillshadeName, DataDirectory,coord_type="UTM_km", colourbar_location='bottom')
    # add the terrace drape
    terrace_cmap = plt.cm.Reds
    MF.add_drape_image(TerraceIDName, DataDirectory, colourmap = terrace_cmap, alpha=0.8)

    ImageName = DataDirectory+fname_prefix+'_terrace_elev.'+FigFormat
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 1000) # Save the figure

#--------------------------------------------------------------------------#
# XY PLOTS
# These functions make plots of the distribution of terrace elevations, etc.
# as you go upstream along the channel
#--------------------------------------------------------------------------#

def make_terrace_swath_plots(DataDirectory,fname_prefix,jn_number, size_format='ESURF', FigFormat='png'):
    """
    Make a plot of the terrace elevations compared to the channel vs the distance upstream_distance
    for a specified junction number.  The points are coloured by a unique terrace ID.

    Args:
        DataDirectory: the data directory of the csv files
        fname_prefix: the DEM string
        jn_number: the outlet junction number of the channel
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.

    Author: FJC
    """
    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        #l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        #l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
        #l_pad = -35

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[10:95,5:80])
    #colorbar axis
    #ax2 = fig.add_subplot(gs[10:95,82:85])

    # read in the csv file and get the columns
    TerraceDF = Helper.ReadTerraceCSV(DataDirectory,fname_prefix,jn_number)

    # make a random colormap
    new_cmap = colours.rand_cmap(len(TerraceDF['PatchID']),type='bright',first_color_black=False,last_color_black=False,verbose=False)

    plt.errorbar(TerraceDF['DistAlongBaseline'], TerraceDF['ChannelRelief'], yerr=TerraceDF['Relief_st_err'], fmt=None, ecolor='k')
    sc = ax.scatter(TerraceDF['DistAlongBaseline'], TerraceDF['ChannelRelief'], marker='D', c=TerraceDF['PatchID'], cmap=new_cmap, s=10, lw=1, edgecolor='k',zorder=100)

    # set labels
    ax.set_xlabel('Distance upstream along main stem (km)')
    ax.set_ylabel('Elevation compared to main stem (m)')
    ax.set_ylim(0,)
    ax.set_xlim(0,)

    # add the colorbar
    # colorbarlabel = "Terrace ID"
    # cbar = plt.colorbar(sc,cmap=new_cmap,spacing='uniform', orientation='vertical',cax=ax2)
    # cbar.set_label(colorbarlabel, fontsize=10)
    # ax2.set_ylabel(colorbarlabel, fontname='Arial', fontsize=10)

    # Save figure
    OutputFigureName = DataDirectory+fname_prefix+'_xy_plots_%s' % str(jn_number)
    plt.savefig(OutputFigureName + '.' + FigFormat, format=FigFormat, dpi=300)
    #plt.close()
    ax.cla()
    plt.clf()

def make_terrace_plots_all_files(DataDirectory,fname_prefix,size_format='ESURF',FigFormat='png'):
    """
    Generate a terrace xy plot for each river in the DEM automatically
    for the junction number.

    Args:
        DataDirectory: the data directory of the csv files
        fname_prefix: the DEM string
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.

    Author: FJC
    """
    import os

    for filename in os.listdir(DataDirectory):
        if filename.endswith(".csv"):
            if "terrace_swath_plots" in filename:
                # split the filename to get the junction number
                this_fname = (filename.split("."))[0]
                this_fname = this_fname.split("_")[-1]
                jn_number = int(this_fname)

                make_terrace_swath_plots(DataDirectory, fname_prefix, jn_number, size_format=size_format, FigFormat=FigFormat)
