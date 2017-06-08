## LSDMap_KnickpointPlotting.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## These functions are tools for analysing and plotting knickpoint data
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## FJC
## 05/06/2017
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import LSDPlottingTools.LSDMap_PointTools as LSDMap_PD
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

def plot_knickpoint_elevations(PointData, DataDirectory, basin_key=0, kp_threshold=0,
                               FigFileName='Image.pdf', FigFormat='pdf', size_format='ESURF', kp_type = "diff"):
    """
    Function to create a plot of knickpoint elevation vs flow distance for each
    basin. Knickpoints are colour-coded by source node, and the marker size represents
    the magnitude of the knickpoint.

    Args:
        PointData: the LSDMap_PointData object with the knickpoint information
        DataDirectory (str): the data directory for the knickpoint file
        csv_name (str): name of the csv file with the knickpoint information
        basin_key (int): key to select the basin of interest
        kp_threshold (int): threshold knickpoint magnitude, any knickpoint below this will be removed (This option may be removed soon)
        kp_type (string): switch between diff and ratio data
        FigFileName (str): The name of the figure file
        FigFormat (str): format of output figure, can be 'pdf' (default), 'png', 'return', or 'show'
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).

    Returns:
        Plot of knickpoint elevations against flow distance

    Author: FJC
    """
    #PointData = LSDMap_PD.LSDMap_PointData(kp_csv_fname)
    # thin out small knickpoints
    KPData = PointData
    #KPData.ThinData(kp_type,kp_threshold)

    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])

    # get the data

    # Soert the dataset to the basin key
    KPData.selectValue('basin_key', value = basin_key, operator = "==")

    elevation = KPData.QueryData('elevation')

    flow_distance = KPData.QueryData('flow distance')

    magnitude = KPData.QueryData(kp_type)
    print("For the plotting, if you want to manage the scale, " +kp_type + " max is "+ str(magnitude.max()) +" and min is " + str(magnitude.min()))

    source = KPData.QueryData('source_key')




    #colour by source - this is the same as the script to colour channels over a raster,
    # (BasicChannelPlotGridPlotCategories) so that the colour scheme should match
    # make a color map of fixed colors
    NUM_COLORS = len(np.unique(source))
    this_cmap = plt.cm.Set1
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    plt.cm.ScalarMappable(norm=cNorm, cmap=this_cmap)
    channel_data = [x % NUM_COLORS for x in source]

    # normalise magnitude for plotting
    min_size = np.min(magnitude)
    max_size = np.max(magnitude)
    #normSize = [100*((x - min_size)/(max_size - min_size)) for x in magnitude]
    normSize = 100 * (magnitude - min_size)/(max_size - min_size)

    # now get the channel profiles that correspond to each knickpoint source and basin
    # PointData.ThinDataFromKey('basin_key',basin_key)
    # PointData.ThinDataSelection('source_key',maskSource)
    #
    # channel_elev = PointData.QueryData('elevation')
    # channel_elev = [float(x) for x in channel_elev]
    # channel_dist = PointData.QueryData('flow_distance')
    # channel_dist = [float(x) for x in channel_dist]

    # now plot the knickpoint elevations and flow distances
    #ax.scatter(channel_dist, channel_elev, s=0.1, c='k')
    ax.scatter(flow_distance, elevation, c = channel_data, cmap=this_cmap, s = normSize, lw=0.5,edgecolors='k',zorder=100)
    ax.set_xlabel('Flow distance (m)')
    ax.set_ylabel('Elevation (m)')

    # return or show the figure
    print("The figure format is: " + FigFormat)
    if FigFormat == 'show':
        plt.show()
    elif FigFormat == 'return':
        return fig
    else:
        save_fmt = FigFormat
        plt.savefig(DataDirectory+FigFileName,format=save_fmt,dpi=500)
        fig.clf()


def plot_diff_ratio(PointData, DataDirectory, saveName = "Basic_diff_ratio", save_fmt = ".png", size_format = "ESURF", log_data = False):
    """
    Basic plot to have a general view of the knickpoints: diff against ratio colored by elevation

    Args:
        PointData: A PointData object
        DataDirectory: Where the data is saved
        saveName: save name

    returns:
        Nothing, sorry.
    Author: BG
    """
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])

    diff = PointData.QueryData("diff")
    ratio = PointData.QueryData("ratio")

    if(log_data):
        print("I am logging the data")
        diff = np.log10(diff)
        ratio = np.log10(ratio)

    elevation =PointData.QueryData("elevation")
    ax.scatter(diff,ratio, s=0.5, lw = 0, c = elevation)
    ax.set_xlabel("Diff")
    ax.set_ylabel("Ratio")

    plt.savefig(DataDirectory+saveName+save_fmt,dpi=500)




def map_custom():
    """
    Testing function to plot custom maps before creating real function for mapping routines

    Args:
        Yes.
    returns:
        No.
    Author:
        BG
    """
    ###### Parameters ######
    Directory = "/home/s1675537/PhD/DataStoreBoris/GIS/Data/Carpathian/knickpoint/" # reading directory
    wDirectory = Directory # writing directory
    Base_file = "Buzau" # It will be the cabkground raster. Each other raster you want to drap on it will be cropped to its extents including nodata
    csv_file = Directory + "test_sign_m_chi.csv" # Name of your point file, add a similar line with different name if you have more than one point file
    DrapeRasterName = "Buzau_hs.bil" # if you want to drap a raster on your background one. Just add a similar line in case you want another raster to drap and so on
    wname = "sign_test" # name of your output file
    dpi = 500 # Quality of your output image, don't exceed 900
    fig_size_inches = 7 # Figure size in Inches

    ##### Now we can load and plot the data

    BackgroundRasterName = Base_file + ".bil" # Ignore this line
    thisPointData = LSDP.LSDMap_PointData(csv_file, PANDEX = True) # Load the point file #1, add a similar line with different name if you have more than one point file.

    plt.clf() # Ignore this line

    MF = MapFigure(BackgroundRasterName, Directory,coord_type="UTM_km") # load the background raster

    MF.add_drape_image(DrapeRasterName,Directory, # Calling the function will add a drapped raster on the top of the background one
                        colourmap = "gray", # colormap used for this raster, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                        alpha=0.5, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                        show_colourbar = False, # Well, this one is explicit I think
                        colorbarlabel = "Colourbar") # Name of your Colourbar, it might bug though



    MF.add_point_data( thisPointData, # this function plot the requested point file using the lat/long column in the csv file
                       column_for_plotting = "m_chi_sign",  # Column used to color the data
                       this_colourmap = "cubehelix", # Colormap used, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                       colorbarlabel = "Colourbar", # Label
                       scale_points = False, # All the point will have the same size if False
                       column_for_scaling = "None", # If scale point True, you can scale the size of your points using one of the columns
                       scaled_data_in_log = False, # If scale point True, you can log the scaling
                       max_point_size = 5, # max size if scale point True again
                       min_point_size = 0.5, # You should be able to guess that one now
                       coulor_log = False, # do you want a log scale for your colorbar ?
                       coulor_manual_scale = [], #Do you want to manually limit the scale of your colorbar? if not let is false
                       manual_size = 0.5, # If none of above is choosen but you want to put another value than 0.5 to scale your point
                       alpha = 1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                       minimum_log_scale_cut_off = -10) # you probably won't need this

    ImageName = wDirectory+str(int(clock.time()))+wname+".png" # Ignore this
    ax_style = "Normal" # Ignore this
    MF.save_fig(fig_width_inches = fig_size_inches, FigFileName = ImageName, axis_style = ax_style, Fig_dpi = dpi) # Save the figure

def map_knickpoint_sign(PointData, DataDirectory, Raster_base_name, HS_name = "none",Time_in_name = False):
    """
    Will create a map of the knickpoint simply colored by sign.

    Args:
        PointData (PointTools object)
        DataDirectory (str): directory where the data will be saved and loaded.
        Raster_base_name (str): Base name of your files without the .bil
        HS_name (str): name of your Hillshade file, by default baseName + _hs.bil like LSDTT create it
        Time_in_name (bool): Option to add timing info in the nae of the figure. Can be useful if you test loads of parameters and you want to be sure that your files names are different (but awful).
    returns:
        No, but creates a map named map_knickpoint_sign.png
    Author:
        BG
    """
    ###### Parameters ######
    Directory = DataDirectory # reading directory
    wDirectory = Directory # writing directory
    Base_file = Raster_base_name # It will be the cabkground raster. Each other raster you want to drap on it will be cropped to its extents including nodata

    wname = "map_knickpoint_sign" # name of your output file
    dpi = 500 # Quality of your output image, don't exceed 900
    fig_size_inches = 7 # Figure size in Inches
    if(HS_name == "none"):
        HS_name = Raster_base_name+("_hs.bil")
    DrapeRasterName = HS_name # if you want to drap a raster on your background one. Just add a similar line in case you want another raster to drap and so on

    ##### Now we can load and plot the data

    BackgroundRasterName = Base_file + ".bil" # Ignore this line
    plt.clf() # Ignore this line

    MF = MapFigure(BackgroundRasterName, Directory,coord_type="UTM_km") # load the background raster

    MF.add_drape_image(DrapeRasterName,Directory, # Calling the function will add a drapped raster on the top of the background one
                        colourmap = "gray", # colormap used for this raster, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                        alpha=0.5, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                        show_colourbar = False, # Well, this one is explicit I think
                        colorbarlabel = "Colourbar") # Name of your Colourbar, it might bug though



    MF.add_point_data( PointData, # this function plot the requested point file using the lat/long column in the csv file
                       column_for_plotting = "sign",  # Column used to color the data
                       this_colourmap = "cubehelix", # Colormap used, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                       colorbarlabel = "Colourbar", # Label
                       scale_points = False, # All the point will have the same size if False
                       column_for_scaling = "None", # If scale point True, you can scale the size of your points using one of the columns
                       scaled_data_in_log = False, # If scale point True, you can log the scaling
                       max_point_size = 5, # max size if scale point True again
                       min_point_size = 0.5, # You should be able to guess that one now
                       coulor_log = False, # do you want a log scale for your colorbar ?
                       coulor_manual_scale = [], #Do you want to manually limit the scale of your colorbar? if not let is false
                       manual_size = 0.5, # If none of above is choosen but you want to put another value than 0.5 to scale your point
                       alpha = 1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                       minimum_log_scale_cut_off = -10) # you probably won't need this
    if(Time_in_name):
        ImageName = wDirectory+str(int(clock.time()))+wname+".png" # Ignore this
    else:
        ImageName = wDirectory+wname+".png" # Ignore this
    ax_style = "Normal" # Ignore this
    MF.save_fig(fig_width_inches = fig_size_inches, FigFileName = ImageName, axis_style = ax_style, Fig_dpi = dpi) # Save the figure
















#
