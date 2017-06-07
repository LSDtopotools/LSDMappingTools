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
