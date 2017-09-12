
"""
Created.

@author: Maxime
"""
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time as clock
from matplotlib import rcParams
import matplotlib.cm as cm
from LSDMapFigure.PlottingRaster import MapFigure
from LSDMapFigure.PlottingRaster import BaseRaster

from LSDPlottingTools import colours as lsdcolours
from LSDPlottingTools import init_plotting_DV
import LSDPlottingTools as LSDP
import sys
import pandas as pd
import numpy as np
import random

#rcParams['font.family'] = 'sans-serif'
#rcParams['font.sans-serif'] = ['arial']
#rcParams['font.size'] = label_size
#rcParams['lines.linewidth']  = 1.5

###### Parameters ######
Directory = "/home/s1793919/Maxime/" # reading directory (if it is on windows, the path is something like C://windows/Users/blablalba/)
wDirectory = "/home/s1793919/Maxime/" # writing directory (if it is on windows, the path is something like C://windows/Users/blablalba/)
Base_file = "north_final" # It will be the cabkground raster. Each other raster you want to drap on it will be cropped to its extents including nodata
csv_file = Directory + "csv/csv_df6_north.txt" # Name of your point file, add a similar line with different name if you have more than one point file
DrapeRasterName = "north_final_hs.bil" # if you want to drap a raster on your background one. Just add a similar line in case you want another raster to drap and so on
wname = "output" # name of your output file
dpi = 600 # Quality of your output image, don't exceed 900
fig_size_inches = 7 # Figure size in Inches

file_name = "/home/s1793919/Maxime/csv/csv_df6_north.txt" # my file is in the same folder than the script so I can just give the name
df = pd.read_csv(file_name)


lat = df["latitude"] # storing the column latitude in a new variable called lat
longit = df["longitude"]
elev = df["elevation"]
flow_dist = df["flow distance"]
drainage = df["drainage area"]
m_chi = df["m_chi"]
b_chi = df["b_chi"]
key = df["source_key"]
basin = df["basin_key"]
seg_elev = df["segmented_elevation"]
seg_num = df["segment_number"]
Geology = df["Geology"]

lithology_name = "/home/s1793919/Maxime/lithology.csv"
df_litho = pd.read_csv(lithology_name, sep = ';', encoding = "ISO-8859-1")


Geol = df_litho["Geology"]
litho = df_litho["Lithology"]
age = df_litho["Age"]


df = df.merge(df_litho, on='Geology', how='inner')

l = list(df["Geology"].unique())
df["GeoCode"] = pd.Series(np.zeros(df.shape[0]))

for i in range(len(l)):
	df["GeoCode"][df["Geology"] == l[i]] = i

GeoCode = df["GeoCode"]





#df2 = df
#df3 = pd.DataFrame()

#i = list(df["source_key"].unique())

#for h in range (len(i)):

    #k = df2[df2["source_key"] == i[h]]
        
    #df3 = k.append(df3)
    
    #h = h+1

    
#print(df3)



#liste = ["flag", "prism", "ocean", "gist_earth", "terrain", "gist_stern", "gnuplot", "gnuplot2", "CMRmap", "cubehelix", "brg", "hsv", "gist_rainbow", "rainbow", "jet", "nipy_spectral", "gist_ncar"]
liste = ["autumn"]

for i in range(len(liste)) :
    
    color = random.choice(liste)
    

##### Now we can load and plot the data

    BackgroundRasterName = Base_file + ".bil" # Ignore this line
    
    thisPointData = LSDP.LSDMap_PointData(df, data_type = "pandas", PANDEX = True) # Load the point file #1, add a similar line with different name if you have more than one point file.
    
    plt.clf() # Ignore this line
    
    MF = MapFigure(BackgroundRasterName, Directory,coord_type="UTM_km", NFF_opti = True, colourbar_location = 'bottom') # load the background raster
    
    MF.add_drape_image(BackgroundRasterName,Directory, # Calling the function will add a drapped raster on the top of the background one
                    colourmap = "gray", # colormap used for this raster, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                    alpha = 0.5, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                    show_colourbar = False, # Well, this one is explicit I think
                    colorbarlabel = "Colourbar", # Name of your Colourbar, it might bug though
                    NFF_opti = True)
                    
    MF.add_drape_image(DrapeRasterName,Directory, # Calling the function will add a drapped raster on the top of the background one
                    colourmap = "gray", # colormap used for this raster, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                    alpha = 0.5, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                    show_colourbar = False, # Well, this one is explicit I think
                    colorbarlabel = "Colourbar", # Name of your Colourbar, it might bug though
                    NFF_opti = True)
                                                
    
    MF.add_point_data( thisPointData, # this function plot the requested point file using the lat/long column in the csv file
                       column_for_plotting = "m_chi",  # Column used to color the data
                       this_colourmap = color, # Colormap used, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                       colorbarlabel = "Colourbar", # Label
                       scale_points = True, # All the point will have the same size if False
                       column_for_scaling = "None", # If scale point True, you can scale the size of your points using one of the columns
                       scaled_data_in_log = False, # If scale point True, you can log the scaling
                       max_point_size = 5, # max size if scale point True again
                       min_point_size = 0.5, # You should be able to guess that one now
                       colour_log = False, # do you want a log scale for your colorbar ?
                       colour_manual_scale = [0,100], #Do you want to manually limit the scale of your colorbar? if not let is false
                       manual_size = 2.5, # If none of above is choosen but you want to put another value than 0.5 to scale your point
                       alpha = 1, # transparency of this specific layer, 0 sfor fully transparent (why not) and 1 for fully opaque
                       minimum_log_scale_cut_off = -10) # you probably won't need this
    
  
    plt.title("Map") 
    ImageName = wDirectory+str(int(clock.time()))+wname+ liste[i] + "new_terrain"+".png" # Ignore this
    
    ax_style = "Normal" # Ignore this
    MF.save_fig(fig_width_inches = fig_size_inches, FigFileName = ImageName, axis_style = ax_style, Fig_dpi = dpi) # Save the figure
