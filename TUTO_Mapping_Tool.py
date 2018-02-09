# -*- coding: utf-8 -*-
"""
Created.

@author: Boris
"""
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import time as clock
from matplotlib import rcParams
import matplotlib.cm as cm
import LSDPlottingTools as LSDP
from LSDPlottingTools import colours as lsdcolours
from LSDPlottingTools import init_plotting_DV
from LSDMapFigure.PlottingRaster import MapFigure
from LSDMapFigure.PlottingRaster import BaseRaster


import pandas as bamboo_bears
import sys

#rcParams['font.family'] = 'sans-serif'
#rcParams['font.sans-serif'] = ['arial']
#rcParams['font.size'] = label_size
#rcParams['lines.linewidth']  = 1.5

###### Parameters ######
Directory = "/home/s1675537/PhD/DataStoreBoris/Emma/" # reading directory (if it is on windows, the path is something like C://windows/Users/blablalba/)
wDirectory = "/home/s1675537/PhD/DataStoreBoris/Emma/" # writing directory (if it is on windows, the path is something like C://windows/Users/blablalba/)
Base_file = "Betics_UTM30clip_PP" # It will be the cabkground raster. Each other raster you want to drap on it will be cropped to its extents including nodata
csv_file = Directory + "new.csv" # Name of your point file, add a similar line with different name if you have more than one point file
DrapeRasterName = "Betics_UTM30clip_hs.bil" # if you want to drap a raster on your background one. Just add a similar line in case you want another raster to drap and so on
wname = "output" # name of your output file
dpi = 300 # Quality of your output image, don't exceed 900
fig_size_inches = 9 # Figure size in Inches



##### Now we can load and plot the data

BackgroundRasterName = Base_file + ".bil" # Ignore this line

# There are several way to load point data: Most of the time you will want to
# directly load a csv_file as follow:
thisPointData = LSDP.LSDMap_PointData(csv_file,data_type = 'csv', PANDEX = True) # Load the point file #1, add a similar line with different name if you have more than one point file.
# But you can also load your file as a pandas dataframe, that allows you to sort your dataset first using pandas:
# df = bamboo_bears.read_csv(csv_file, sep = ',')
# df = df[df['elevation'] >50] # I am selecting all the points over 50 meters
# thisPointData = LSDP.LSDMap_PointData(csv_file,data_type = 'pandas', PANDEX = True) # I am loading diretly a pandas dataframe rather than the csv_file


######## Now plotting the different layers ###########
# Note that you can add as much rasters and points data as you want by just copying-pasting the add_drape_image and add_point_data functions
# You just have to adapt the parameters and load the new data earlier in the script

plt.clf() # Ignore this line, this is to make sure you're creating a brand new figure
MF = MapFigure(BackgroundRasterName, Directory,coord_type="UTM_km", NFF_opti = True) # load the background raster

MF.add_drape_image(DrapeRasterName,Directory, # Calling the function will add a drapped raster on the top of the background one
                    colourmap = "gray", # colormap used for this raster, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                    alpha = 0.5, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                    show_colourbar = False, # Well, this one is explicit I think
                    colorbarlabel = "Colourbar", # Name of your Colourbar, it might bug though
                    NFF_opti = True,
                    colour_min_max = [])



MF.add_point_data( thisPointData, # this function plot the requested point file using the lat/long column in the csv file
                   column_for_plotting = "chi",  # Column used to color the data
                   this_colourmap = "cubehelix", # Colormap used, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                   colorbarlabel = "Colourbar", # Label
                   scale_points = False, # All the point will have the same size if False
                   column_for_scaling = "None", # If scale point True, you can scale the size of your points using one of the columns
                   scaled_data_in_log = False, # If scale point True, you can log the scaling
                   max_point_size = 5, # max size if scale point True again
                   min_point_size = 0.5, # You should be able to guess that one now
                   colour_log = False, # do you want a log scale for your colorbar ?
                   colour_manual_scale = [], #Do you want to manually limit the scale of your colorbar? if not let is false
                   manual_size = 0.5, # If none of above is choosen but you want to put another value than 0.5 to scale your point
                   alpha = 1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                   minimum_log_scale_cut_off = -10) # you probably won't need this


ImageName = wDirectory+str(int(clock.time()))+wname+".png" # Ignore this
ax_style = "Normal" # Ignore this
MF.save_fig(fig_width_inches = fig_size_inches, FigFileName = ImageName, axis_style = ax_style, Fig_dpi = dpi) # Save the figure
