"""
Evia fault plots

MDH 

"""

#import modules
#import numpy as np

#import matplotlib.ticker as ticker
#import pandas as pd
#from matplotlib import colors
#import math

#import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
import sys, os
import LSDPlottingTools as LSDP
from LSDMapFigure import PlottingHelpers as Helper
from LSDMapFigure.PlottingRaster import MapFigure
from LSDMapFigure.PlottingRaster import BaseRaster
#from shapely.geometry import Polygon

def print_welcome():

    print("\n\n=======================================================================")
    print("Hello! I'm going to plot some basins")
    print("You will need to tell me which directory to look in.")
    print("Use the -dir flag to define the working directory.")
    print("If you don't do this I will assume the data is in the same directory as this script.")
    print("=======================================================================\n\n ")

def main(argv):

    # If there are no arguments, send to the welcome screen
    if not len(sys.argv) > 1:
        full_paramfile = print_welcome()
        sys.exit()

    # Get the arguments
    import argparse
    parser = argparse.ArgumentParser()
    
    # The location of the data files
    parser.add_argument("-dir", "--base_directory", type=str, help="The base directory. If not defined, current directory.")
    parser.add_argument("-fname", "--fname_prefix", type=str, help="The prefix of your DEM WITHOUT EXTENSION!")
    parser.add_argument("-fmt", "--FigFormat", type=str, default='png', help="Set the figure format for the plots. Default is png")
    args = parser.parse_args()
    
    # get the base directory
    if args.base_directory:
        DataDirectory = args.base_directory
        # check if you remembered a / at the end of your path_name
        if not DataDirectory.endswith("/"):
            print("You forgot the '/' at the end of the directory, appending...")
            DataDirectory = this_dir+"/"
    else:
        this_dir = os.getcwd()
    
    if not args.fname_prefix:
        print("WARNING! You haven't supplied your DEM name. Please specify this with the flag '-fname'")
        sys.exit()
    else:
        fname_prefix = args.fname_prefix
        
    # set to not parallel
    parallel = False
    faults = True
    FigFormat = args.FigFormat
   
    # check if a directory exists for the chi plots. If not then make it.
    raster_directory = DataDirectory+'raster_plots/'
    if not os.path.isdir(raster_directory):
        os.makedirs(raster_directory)

    # Set up fonts for plots
    label_size = 8
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # set figure sizes based on format
    size_format = ""
    if size_format == "geomorphology":
        fig_width_inches = 6.25
    elif size_format == "big":
        fig_width_inches = 16
    else:
        fig_width_inches = 4.92126

    # get the basin IDs to make a discrete colourmap for each ID
    BasinInfoDF = Helper.ReadBasinInfoCSV(DataDirectory, fname_prefix)
    
    basin_keys = list(BasinInfoDF['basin_key'])
    basin_keys = [int(x) for x in basin_keys]

    basin_junctions = list(BasinInfoDF['outlet_junction'])
    basin_junctions = [int(x) for x in basin_junctions]

    # get a discrete colormap
    cmap = plt.cm.viridis

    # going to make the basin plots - need to have bil extensions.
    print("I'm going to make the basin plots. Your topographic data must be in ENVI bil format or I'll break!!")

    # get the rasters
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext
    BasinsName = fname_prefix+'_AllBasins'+raster_ext

    # create the map figure
    MF = MapFigure(HillshadeName, DataDirectory,coord_type="UTM_km", colourbar_location='none')

    # add the basins drape
    BasinsDict = dict(zip(basin_keys,basin_keys))
#    MF.add_basin_plot(BasinsName,fname_prefix,DataDirectory, label_basins=False,
#                      use_keys_not_junctions = True, show_colourbar = False, 
#                      value_dict = BasinsDict, discrete_cmap=True, n_colours=len(basin_keys),
#                      colorbarlabel = "Basin ID", cbar_type=int, tickspacefactor=2,
#                      colourmap = cmap, edgecolour='none', adjust_text = True, parallel=parallel)

    # add the channel network
    if not parallel:
        ChannelDF = Helper.ReadChiDataMapCSV(DataDirectory,fname_prefix)
    else:
        ChannelDF = Helper.AppendChiDataMapCSVs(DataDirectory)
    
    # remove chi no data values
    ChannelDF = ChannelDF[ChannelDF.chi != -9999]
    
    ChannelPoints = LSDP.LSDMap_PointData(ChannelDF, data_type = "pandas", PANDEX = True)
    
	# add chi map
    MF.add_point_data(ChannelPoints, column_for_plotting = "chi", column_for_scaling = "chi", colorbarlabel = "$\chi$ (m)", show_colourbar=True, this_colourmap = cmap, colourbar_location="top")
	
    # add the faults
    if faults:
        LineFileName = DataDirectory + fname_prefix + "_faults.shp"
        MF.add_line_data(LineFileName, linestyle="-", linewidth=1.5, zorder=99, legend=True, label="Fault Segments")
        
    # add the basin outlines ### need to parallelise
    if not parallel:
      Basins = LSDP.GetBasinOutlines(DataDirectory, BasinsName)
    else:
      Basins = LSDP.GetMultipleBasinOutlines(DataDirectory)
    
    # Find the relay basins and plot separately
    RelayBasinIDs = [1248, 4788, 4995, 5185, 6187, 6758, 6805]
    RelayBasins = {key:value for key, value in Basins.items() if key in RelayBasinIDs}

    # Plot all basins
    MF.plot_polygon_outlines(Basins, colour='k', linewidth=0.5, alpha = 1, legend=True, label="Catchments")
    MF.plot_polygon_outlines(RelayBasins, colour='r', linewidth=0.5, alpha = 1,legend=True, label="Relay Catchments")

    # Add the legend
    MF.add_legend()
    
    # Save the figure
    ImageName = raster_directory+fname_prefix+'_chi_map.'+FigFormat
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 300)

    # Make m/n summary plots
    
if __name__ == "__main__":
    main(sys.argv[1:])

