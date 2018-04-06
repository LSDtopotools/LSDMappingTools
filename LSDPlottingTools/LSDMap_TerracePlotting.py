#------------------------------------------------------------------------------#
# terrace-long-profiler
# Scripts to plot the long profiles of river terraces compared to the main
# channel.
# Authors: F. Clubb
#          A. Wickert
#------------------------------------------------------------------------------#

# import modules
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from LSDPlottingTools import colours
import matplotlib.cm as cm
from matplotlib import rcParams
from matplotlib import colors as colors
from LSDPlottingTools import LSDMap_GDALIO as IO
from LSDMapFigure import PlottingHelpers as H
from shapely.geometry import shape, Polygon, Point, LineString
import fiona
import os

#---------------------------------------------------------------------------------------------#
# Set up figure
#---------------------------------------------------------------------------------------------#
def CreateFigure(FigSizeFormat="default", AspectRatio=16./9.):
    """
    This function creates a default matplotlib figure object

    Args:
        FigSizeFormat: the figure size format according to journal for which the figure is intended
            values are geomorphology,ESURF, ESPL, EPSL, JGR, big
            default is ESURF

        AspectRatio: The shape of the figure determined by the aspect ratio, default is 16./9.

    Returns:
        matplotlib figure object

    Author: FJC
    """
    # set figure sizes (in inches) based on format
    if FigSizeFormat == "geomorphology":
        FigWidth_Inches = 6.25
    elif FigSizeFormat == "big":
        FigWidth_Inches = 16
    elif FigSizeFormat == "ESURF":
        FigWidth_Inches = 4.92
    elif FigSizeFormat == "ESPL":
        FigWidth_Inches = 7.08
    elif FigSizeFormat == "EPSL":
        FigWidth_Inches = 7.48
    elif FigSizeFormat == "JGR":
        FigWidth_Inches = 6.6

    else:
        FigWidth_Inches = 4.92126

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = 10
    #rcParams['text.usetex'] = True

    Fig = plt.figure(figsize=(FigWidth_Inches,FigWidth_Inches/AspectRatio))

    return Fig
#---------------------------------------------------------------------------------------------#
# ANALYSIS FUNCTIONS
# Functions to analyse the terrace info
#---------------------------------------------------------------------------------------------#

def SelectTerracesFromShapefile(DataDirectory,shapefile_name,fname_prefix):
    """
    This function takes in a shapefile of digitised terraces and
    uses it to filter the terrace DF.  Only pixels within each
    shapefile are kept, and they are assigned a new ID based on the
    ID of the shapefile polygons.

    Args:
        DataDirectory (str): the data directory
        shapefile_name (str): the name of the shapefile
        fname_prefix (str): prefix of the DEM

    Returns: terrace df filtered by the digitised terraces

    Author: FJC
    """
    # first get the terrace df
    terrace_df = H.read_terrace_csv(DataDirectory,fname_prefix)

    # now get the shapefile with the digitised terraces
    digitised_terraces = H.read_terrace_shapefile(DataDirectory,shapefile_name)

    # for each point in the df, need to check if it is in one of the polygons. This will probably be slow.

    # set up the new terrace df
    new_df = pd.DataFrame()

    print ("Filtering points by shapefile, this might take a while...")

    for idx, row in terrace_df.iterrows():
        this_point = Point(row['X'], row['Y'])
        #print this_point
        # check if this point is in one of the polygons
        for id, polygon in digitised_terraces.items():
            if polygon.contains(this_point):
                # this point is within this terrace, keep it and assign a new ID number
                row['TerraceID'] = id
                new_df = new_df.append(row)

    OutDF_name = "_terrace_info_shapefiles.csv"
    OutDF_name = DataDirectory+fname_prefix+OutDF_name
    new_df.to_csv(OutDF_name,index=False)

    return new_df

def SelectTerracePointsFromCentrelines(DataDirectory,shapefile_name,fname_prefix, distance=2):
    """
    This function takes in a shapefile of digitised terrace centrelines and finds points within a certain distance of the line. Returns as a df.

    Args:
        DataDirectory (str): the data directory
        shapefile_name (str): the name of the shapefile
        fname_prefix (str): prefix of the DEM

    Returns: terrace df filtered by the digitised terraces

    Author: FJC
    """
    # first get the terrace df
    terrace_df = H.read_terrace_csv(DataDirectory,fname_prefix)

    # now get the shapefile with the digitised terraces
    centrelines = H.read_terrace_centrelines(DataDirectory,shapefile_name)

    # for each point in the df, need to check if it is in one of the polygons. This will probably be slow.

    # set up the new terrace df
    new_df = pd.DataFrame()

    print ("Filtering points by shapefile, this might take a while...")

    for idx, row in terrace_df.iterrows():
        this_point = Point(row['X'], row['Y'])
        #print this_point
        # check if this point is in one of the polygons
        for id, line in centrelines.items():
            if line.distance(this_point) < distance:
                # this point is within this terrace, keep it and assign a new ID number
                row['TerraceID'] = id
                new_df = new_df.append(row)

    OutDF_name = "_terrace_info_centrelines.csv"
    OutDF_name = DataDirectory+fname_prefix+OutDF_name
    new_df.to_csv(OutDF_name,index=False)

    return new_df

def filter_terraces(terrace_df,min_size=5000, max_size=1000000):
    """
    This function takes the initial terrace dataframe and sorts it to remove terraces
    that are too small.

    Args:
        terrace_df: pandas dataframe with the terrace info
        min_size (int): minimum n of pixels in each terrace
        max_size (int): max n of pixels in each terrace

    Returns:
        dataframe with filtered terrace info.

    Author: FJC
    """
    # first get the unique terrace IDs
    terraceIDs = terrace_df.TerraceID.unique()

    # loop through unique IDs and check how many rows correspond to this ID, then
    # remove any that are too small
    for i in terraceIDs:
        n_pixels = len(terrace_df[terrace_df['TerraceID'] == i])
        if n_pixels < min_size or n_pixels > max_size:
            terrace_df = terrace_df[terrace_df['TerraceID'] != i]

    return terrace_df

def write_dip_and_dipdir_to_csv(DataDirectory,fname_prefix, digitised_terraces=False, shapefile_name=None):
    """
    Wrapper for dip and dipdir function

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): name of the DEM
        digitised_terraces (bool): boolean to use digitised terrace shapefile
        shapefile_name (str): name of shapefile

    Author: FJC
    """
    # read in the terrace csv
    terraces = H.read_terrace_csv(DataDirectory,fname_prefix)
    if digitised_terraces:
        # check if you've already done the selection, if so just read in the csv
        print ("File name is", DataDirectory+fname_prefix+'_terrace_info_shapefiles.csv')
        if os.path.isfile(DataDirectory+fname_prefix+'_terrace_info_shapefiles.csv'):
            terraces = pd.read_csv(DataDirectory+fname_prefix+'_terrace_info_shapefiles.csv')
        else:
            terraces = SelectTerracesFromShapefile(DataDirectory,shapefile_name,fname_prefix)
    else:
        filter_terraces(terraces, min_size)

    # get the terrace dip and dip dirs
    terrace_dips = get_terrace_dip_and_dipdir(terraces)

    # write to csv
    terrace_dips.to_csv(DataDirectory+fname_prefix+'_Dip_DipDirection.csv')


def get_terrace_dip_and_dipdir(terrace_df):
    """
    This function takes the initial terrace dataframe and calculates the dip and
    strike of the terrace surfaces. Fits a polynomial surface to the distribution
    of terrace elevations and then gets the dip and dip directions of this surface.

    Args:
        terrace_df: pandas dataframe with the terrace info

    Returns:
        dataframe with terrace dip and dip directions

    Author: AW and FJC
    """
    from scipy import linalg
    import math
    from mpl_toolkits.mplot3d import Axes3D

    # get the unique terrace IDs
    terraceIDs = terrace_df.TerraceID.unique()
    print (terraceIDs)

    dips = []
    dip_dirs = []
    strikes = []
    XbarTerraces = []
    YbarTerraces = []

    # get the info for each terrace ID
    for terraceID in terraceIDs:
        _terrace_subset = (terrace_df.TerraceID.values == terraceID)
        _x = terrace_df['DistAlongBaseline'].values[_terrace_subset]
        _y = terrace_df['DistToBaseline'].values[_terrace_subset]
        _z = terrace_df['Elevation'].values[_terrace_subset]
        _X = terrace_df['X'].values[_terrace_subset]
        _Y = terrace_df['Y'].values[_terrace_subset]

        # fit a plane to these points
        # form: Z = C[0]*X + C[1]*Y + C[2]
        _XY = np.vstack((_X, _Y, np.ones(len(_Y)))).transpose()
        C,_,_,_ = linalg.lstsq(_XY, _z)

        # going to get the dip and dip direction using the unit normal vector
        # to the plane.
        a = -C[0]
        b = -C[1]
        c = 1
        n_vec = np.array([a,b,c])
        print (n_vec)
        # n vector projected onto the xy plane (multiply n_vec by (1,1,0))
        n_xy = np.array([a,b,0])

        #print n_vec

        # now get the dip = angle between n_vec and n_xy. angle between 2 vectors:
        # cos theta = (alpha . beta) / (|alpha| |beta|)
        # This gives the dip in radians
        dip = np.arccos((np.dot(n_vec,n_xy))/(np.linalg.norm(n_vec)*np.linalg.norm(n_xy)))
        # convert to degrees
        dip = 90 - math.degrees(dip)

        # get the dip direction = angle between n_proj and the due north vector y = (0,1,0)
        y_vec = np.array([0,1,0])
        theta = np.arccos((np.dot(n_xy,y_vec))/(np.linalg.norm(n_xy)*np.linalg.norm(y_vec)))
        theta = math.degrees(theta)
        print ("Theta", theta)

        # work out strike depending on orientation
        if a > 0 and b > 0: # x is positive so dip dir is just theta
            strike = 270 + theta
        elif a < 0:
            strike = 270 - theta
        elif a > 0 and b < 0: # x is negative so dip dir is 180+theta?
            strike = theta - 90

        # now get the dip dir using the right hand rule
        dip_dir = strike+90
        if dip_dir > 359:
            dip_dir = dip_dir - 360

        dips.append(dip)
        dip_dirs.append(dip_dir)
        strikes.append(strike)

        # get the mean x and y for plotting
        XbarTerraces.append(np.mean(_X))
        YbarTerraces.append(np.mean(_Y))

    outarray = np.vstack((XbarTerraces, YbarTerraces, dips, dip_dirs, strikes)).transpose()
    _column_names = ('X', 'Y', 'dip', 'dip_azimuth', 'strike')
    _index = np.arange(len(dips))+1
    output_pd = pd.DataFrame(data = outarray, index=_index, columns=_column_names)
    return output_pd

def get_terrace_areas(terrace_df, fname_prefix):
    """
    This function takes the initial terrace dataframe and calculates the
    area of each terrace.

    Args:
        terrace_df: pandas dataframe with the terrace info
        fname_prefix: name of the DEM (to get data res)

    Returns:
        dict where key is the terrace ID and value is the terrace area in m^2

    Author: FJC
    """
    # get unique IDs
    terraceIDs = terrace_df.terraceID.unique()

    area_dict = {}

    for terraceID in terraceIDs:
        # get the n rows with this ID
        masked_df = terrace_df[terrace_df['terraceID'] == terraceID]
        n_pixels = len(masked_df.index)

        # get the data resolution of the DEM
        Cell_area = IO.GetPixelArea(fname_prefix)
        terrace_area = n_pixels * Cell_area

        area_dict[terraceID] = terrace_area

    return area_dict

#---------------------------------------------------------------------------------------------#
# XZ PLOTS
# Functions to make XZ plots of terraces
#---------------------------------------------------------------------------------------------#

def long_profiler(DataDirectory,fname_prefix, min_size=5000, FigFormat='png', size_format='ESURF'):
    """
    This function creates a plot of the terraces with distance
    downstream along the main channel.

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the name of the DEM
        min_size (int): the minimum number of pixels for a terrace. Any smaller ones will be excluded.
        FigFormat: the format of the figure, default = png
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).

    Returns:
        terrace long profile plot

    Author: AW, FJC
    """
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

    # make the plot
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.05,right=0.85,top=1.0)
    ax = fig.add_subplot(gs[5:100,10:95])

    # read in the terrace csv
    terraces = H.read_terrace_csv(DataDirectory,fname_prefix)
    filter_terraces(terraces, min_size)

    # read in the baseline channel csv
    lp = H.read_channel_csv(DataDirectory,fname_prefix)
    lp = lp[lp['Elevation'] != -9999]

    terraceIDs = sorted(list(set(list(terraces.TerraceID))))
    xTerraces = []
    zTerraces = []
    yTerraces = []
    newIDs = []

    # loop through the terrace IDs and get the x, y, and z values
    for terraceID in terraceIDs:
        _terrace_subset = (terraces.TerraceID.values == terraceID)
        _x = terraces['DistAlongBaseline'].values[_terrace_subset]
        _y = terraces['DistToBaseline'].values[_terrace_subset]
        _z = terraces['Elevation'].values[_terrace_subset]
        _x_unique = sorted(list(set(list(_x))))
        _z_unique = []
        # Filter
        if len(_x) > 50 and len(_x_unique) > 1 and len(_x_unique) < 1000:
            #print np.max(np.diff(_x_unique))
            #if len(_x_unique) > 10 and len(_x_unique) < 1000 \
            #and :
            for _x_unique_i in _x_unique:
                #_y_unique_i = np.min(np.array(_y)[_x == _x_unique_i])
                #_z_unique.append(np.min(_z[_y == _y_unique_i]))
                _z_unique.append(np.min(_z[_x == _x_unique_i]))
            if np.mean(np.diff(_z_unique)/np.diff(_x_unique)) < 10:
                xTerraces.append(_x_unique)
                zTerraces.append(_z_unique)
                newIDs.append(terraceID)

    # get discrete colours so that each terrace is a different colour
    this_cmap = cm.rainbow
    this_cmap = colours.cmap_discretize(len(newIDs),this_cmap)
    print ("N COLOURS: ", len(newIDs))
    print (newIDs)
    colors = iter(this_cmap(np.linspace(0, 1, len(newIDs))))
    # plot the terraces
    for i in range(len(xTerraces)):
        plt.scatter(xTerraces[i], zTerraces[i], s=2, c=next(colors))

    # plot the main stem channel in black
    plt.plot(lp['DistAlongBaseline'],lp['Elevation'], c='k', lw=2)

    # add a colourbar
    cax = fig.add_axes([0.83,0.15,0.03,0.8])
    sm = plt.cm.ScalarMappable(cmap=this_cmap, norm=plt.Normalize(vmin=min(newIDs), vmax=max(newIDs)))
    sm._A = []
    cbar = plt.colorbar(sm,cmap=this_cmap,spacing='uniform',cax=cax, label='Terrace ID', orientation='vertical')
    colours.fix_colourbar_ticks(cbar,len(newIDs),cbar_type=int,min_value=min(newIDs),max_value=max(newIDs),labels=newIDs)

    # set axis params and save
    ax.set_xlabel('Distance downstream (m)')
    ax.set_ylabel('Elevation (m)')
    plt.savefig(DataDirectory+fname_prefix+'_terrace_plot.'+FigFormat,format=FigFormat,dpi=300)
    plt.clf()

def long_profiler_dist(DataDirectory,fname_prefix, min_size=5000, FigFormat='png', size_format='ESURF', digitised_terraces=False, shapefile_name=None):
    """
    Make long profile plot where terrace points are binned by
    distance along the channel
    """
    # make a figure
    fig = CreateFigure()
    ax = plt.subplot(111)

    # read in the terrace csv
    terraces = H.read_terrace_csv(DataDirectory,fname_prefix)
    if digitised_terraces:
        # check if you've already done the selection, if so just read in the csv
        if os.path.isfile(DataDirectory+fname_prefix+'_terrace_info_shapefiles.csv'):
            terraces = pd.read_csv(DataDirectory+fname_prefix+'_terrace_info_shapefiles.csv')
        else:
            terraces = SelectTerracesFromShapefile(DataDirectory,shapefile_name,fname_prefix)
    else:
        filter_terraces(terraces, min_size)

    # read in the baseline channel csv
    lp = H.read_channel_csv(DataDirectory,fname_prefix)
    lp = lp[lp['Elevation'] != -9999]

    # get the distance from outlet along the baseline for each terrace pixels
    new_terraces = terraces.merge(lp, left_on = "BaselineNode", right_on = "node")
    print (new_terraces)

    xTerraces = np.array(new_terraces['DistFromOutlet'])
    yTerraces = np.array(new_terraces['DistToBaseline'])
    zTerraces = np.array(new_terraces['Elevation_x'])

    MaximumDistance = xTerraces.max()

    # now bin by distance along the baseline
    bins = np.unique(xTerraces)
    nbins = len(np.unique(xTerraces))
    n, _ = np.histogram(xTerraces, bins=nbins)
    s_zTerraces, _ = np.histogram(xTerraces, bins=nbins, weights=zTerraces)
    s_zTerraces2, _ = np.histogram(xTerraces, bins=nbins, weights=zTerraces*zTerraces)
    mean = s_zTerraces / n
    std = np.sqrt(s_zTerraces2/n - mean*mean)

    # # invert to get distance from outlet
    # MS_DistAlongBaseline = np.array(lp['DistAlongBaseline'])[::-1]
    MS_Dist = np.array(lp['DistFromOutlet'])
    MS_Elevation = np.array(lp['Elevation'])
    Terrace_Elevation = mean

    print (MS_Dist)
    print (MS_Elevation)
    print (Terrace_Elevation)


    # plot the main stem channel in black
    plt.plot(MS_Dist/1000,MS_Elevation, c='k', lw=1)
    plt.scatter((_/1000)[:-1], Terrace_Elevation, s=2, zorder=2, c='r')

    # set axis params and save
    ax.set_xlabel('Distance from outlet (km)')
    ax.set_ylabel('Elevation (m)')
    ax.set_xlim(0,80)
    plt.tight_layout()
    plt.savefig(DataDirectory+fname_prefix+'_terrace_plot_binned.'+FigFormat,format=FigFormat,dpi=300)

    plt.clf()

def long_profiler_centrelines(DataDirectory,fname_prefix, shapefile_name, colour_by_ksn=False, ages="", FigFormat='png'):
    """
    Function takes in the csv file of terrace centreline data
    and plots as a long profile against the baseline channel.

    Author: FJC
    """
    # check if a directory exists for the chi plots. If not then make it.
    T_directory = DataDirectory+'terrace_plots/'
    if not os.path.isdir(T_directory):
        os.makedirs(T_directory)

    # make a figure
    fig = CreateFigure()
    ax = plt.subplot(111)

    # check if the csv already exists. if not then select the points from the centrelines
    csv_filename = DataDirectory+fname_prefix+'_terrace_info_centrelines.csv'
    if os.path.isfile(csv_filename):
        terrace_df = pd.read_csv(DataDirectory+fname_prefix+'_terrace_info_centrelines.csv')
    else:
        terrace_df = SelectTerracePointsFromCentrelines(DataDirectory,shapefile_name,fname_prefix, distance=5)

    # read in the baseline channel csv
    lp = H.ReadMChiSegCSV(DataDirectory,fname_prefix)
    lp = lp[lp['elevation'] != -9999]

    # get the distance from outlet along the baseline for each terrace pixels
    terrace_df = terrace_df.merge(lp, left_on = "BaselineNode", right_on = "node")

    # make the long profile plot
    terrace_ids = terrace_df['TerraceID'].unique()

    # sort out the colours. We want a different colour for each terrace...
    this_cmap = cm.viridis
    norm = colors.Normalize(vmin=0, vmax=0.01, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=this_cmap)

    terrace_gradients = []

    for id in terrace_ids:
        # mask for this ID
        this_df = terrace_df[terrace_df['TerraceID'] == id]
        this_dist = this_df['flow_distance']/1000
        this_elev = this_df['Elevation']
        # work out the number of bins. We want a spacing of ~ 10 m
        bins = np.arange(this_dist.min(), this_dist.max(), 0.05)
        n_bins = len(bins)
        if n_bins < 1: n_bins = 1
        # bin the data
        n, _ = np.histogram(this_dist, bins=n_bins)
        sy, _ = np.histogram(this_dist, bins=n_bins, weights = this_elev)
        sy2, _ = np.histogram(this_dist, bins=n_bins, weights = this_elev*this_elev)
        mean = sy/n

        # work out the gradient of each terrace ID...rise/run. Use this to colour the terrace.
        delta_x = this_df['flow_distance'].max() - this_df['flow_distance'].min()
        delta_z = this_elev.max() - this_elev.min()
        gradient = delta_z / delta_x
        color=mapper.to_rgba(gradient)
        terrace_gradients.append(gradient)

        # plot the terrace
        ax.plot((_[1:] + _[:-1])/2, mean,c=color,zorder=2)

    lp_mainstem = H.read_index_channel_csv(DataDirectory,fname_prefix)
    lp_mainstem = lp_mainstem[lp_mainstem['elevation'] != -9999]
    lp_mainstem = lp_mainstem.merge(lp, left_on="id", right_on="node")
    if colour_by_ksn == True:
        ax.scatter(lp_mainstem["flow_distance_x"]/1000, lp_mainstem["elevation_x"], c=lp_mainstem["m_chi"], cmap=cm.hot, norm=colors.Normalize(lp_mainstem["m_chi"].min(), lp_mainstem["m_chi"].max()), s=0.5, lw=0.1)
    else:
        ax.plot(lp_mainstem['flow_distance_x']/1000,lp_mainstem['elevation_x'],'k',lw=1,label='_nolegend_')

    # if present, plot the ages on the profile
    print (ages)
    if ages:
        # read in the ages csv
        ages_df = pd.read_csv(DataDirectory+ages)
        upstream_dist = list(ages_df['upstream_dist'])
        elevation = list(ages_df['elevation'])
        ax.scatter(upstream_dist, elevation, s=8, c="w", edgecolors="k", label="$^{14}$C age (cal years B.P.)")
        ax.legend(loc='upper left', fontsize=8, numpoints=1)

    # set axis params and save
    ax.set_xlabel('Flow distance (km)')
    ax.set_ylabel('Elevation (m)')
    ax.set_xlim(0,(terrace_df['flow_distance'].max()/1000))
    ax.set_ylim(0,terrace_df['elevation'].max()+10)

    # add a colourbar
    mapper.set_array(terrace_gradients)
    cbar = plt.colorbar(mapper,cmap=this_cmap,norm=norm,orientation='vertical')
    cbar.set_label('Gradient (m/m)')


    plt.tight_layout()
    #plt.show()
    plt.savefig(T_directory+fname_prefix+'_terrace_plot_centrelines.'+FigFormat,format=FigFormat,dpi=300)

#---------------------------------------------------------------------------------------------#
# Terrace plots in chi space
# FOR THE FIRST TIME! WOOOOOP
#---------------------------------------------------------------------------------------------#
def MakeTerracePlotChiSpace(DataDirectory,fname_prefix,shapefile_name, colour_by_ksn=True):
    """
    This function makes a plot of the terraces in chi-elevation space. The elevation
    of each terrace is plotted based on the chi of the nearest channel
    FJC 21/03/18
    """
    # check if a directory exists for the chi plots. If not then make it.
    T_directory = DataDirectory+'terrace_plots/'
    if not os.path.isdir(T_directory):
        os.makedirs(T_directory)

    # make a figure
    fig = CreateFigure()
    ax = plt.subplot(111)

    # check if the csv already exists. if not then select the points from the centrelines
    csv_filename = DataDirectory+fname_prefix+'_terrace_info_centrelines.csv'
    if os.path.isfile(csv_filename):
        terrace_df = pd.read_csv(DataDirectory+fname_prefix+'_terrace_info_centrelines.csv')
    else:
        terrace_df = SelectTerracePointsFromCentrelines(DataDirectory,shapefile_name,fname_prefix, distance=5)

    # read in the mchi  csv
    lp = H.ReadMChiSegCSV(DataDirectory,fname_prefix)
    lp = lp[lp['elevation'] != -9999]

    # get the distance from outlet along the baseline for each terrace pixels
    terrace_df = terrace_df.merge(lp, left_on = "BaselineNode", right_on = "node")

    # make the long profile plot
    terrace_ids = terrace_df['TerraceID'].unique()

    # sort out the colours. We want a different colour for each terrace...
    this_cmap = cm.viridis
    norm = colors.Normalize(vmin=terrace_df['m_chi'].min(), vmax=terrace_df['m_chi'].max()-10, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=this_cmap)

    terrace_gradients = []

    for id in terrace_ids:
        # mask for this ID
        this_df = terrace_df[terrace_df['TerraceID'] == id]
        this_chi = this_df['chi']
        this_elev = this_df['Elevation']
        # work out the number of bins. We want a spacing of ~ 10 m
        bins = np.arange(this_chi.min(), this_chi.max(), 0.01)
        n_bins = len(bins)
        if n_bins < 1: n_bins = 1
        # bin the data
        n, _ = np.histogram(this_chi, bins=n_bins)
        sy, _ = np.histogram(this_chi, bins=n_bins, weights = this_elev)
        sy2, _ = np.histogram(this_chi, bins=n_bins, weights = this_elev*this_elev)
        mean = sy/n

        # work out the gradient of each terrace ID...rise/run. Use this to colour the terrace.
        delta_x = this_df['chi'].max() - this_df['chi'].min()
        delta_z = this_elev.max() - this_elev.min()
        gradient = delta_z / delta_x
        color=mapper.to_rgba(gradient)
        terrace_gradients.append(gradient)

        # plot the terrace
        ax.plot((_[1:] + _[:-1])/2, mean,c=color,zorder=2)

    lp_mainstem = H.read_index_channel_csv(DataDirectory,fname_prefix)
    lp_mainstem = lp_mainstem[lp_mainstem['elevation'] != -9999]
    lp_mainstem = lp_mainstem.merge(lp, left_on="id", right_on="node")
    if colour_by_ksn == True:
        ax.scatter(lp_mainstem["chi"], lp_mainstem["elevation_y"], c=lp_mainstem["m_chi"], cmap=cm.viridis, norm=colors.Normalize(lp_mainstem["m_chi"].min()-10, lp_mainstem["m_chi"].max()), s=0.5, lw=0.1)
    else:
        ax.plot(lp_mainstem['chi'],lp_mainstem['elevation_y'],'k',lw=1)

    # set axis params and save
    ax.set_xlabel('$\chi$')
    ax.set_ylabel('Elevation (m)')
    ax.set_xlim(0,(terrace_df['chi'].max()))
    ax.set_ylim(0,terrace_df['Elevation'].max()+10)

    # add a colourbar
    mapper.set_array(terrace_gradients)
    cbar = plt.colorbar(mapper,cmap=this_cmap,norm=norm,orientation='vertical')
    cbar.set_label('$k_{sn}$')


    plt.tight_layout()
    #plt.show()
    plt.savefig(T_directory+fname_prefix+'_terrace_plot_chi.png',format='png',dpi=300)

def MakeTerraceHeatMap(DataDirectory,fname_prefix, mchi_fname, prec=100, bw_method=0.03, FigFormat='png', ages=""):
    """
    Function to make a heat map of the terrace pixels using Gaussian KDE.
    see https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.gaussian_kde.html
    for more details.

    Args:
        DataDirectory(str): the data directory
        fname_prefix(str): prefix of your DEM
        prec(int): the resolution for the KDE. Increase this to get a finer resolution, decrease for coarser.
        bw_method: the method for determining the bandwidth of the KDE.  This is apparently quite sensitive to this.
        Can either be "scott", "silverman" (where the bandwidth will be determined automatically), or a scalar. Default = 0.03
        FigFormat(str): figure format, default = png
        ages (str): Can pass in the name of a csv file with terrace ages which will be plotted on the profile. Must be in the same directory

    FJC 26/03/18
    """
    import scipy.stats as st

    # check if a directory exists for the chi plots. If not then make it.
    T_directory = DataDirectory+'terrace_plots/'
    if not os.path.isdir(T_directory):
        os.makedirs(T_directory)

    # make a figure
    fig = CreateFigure()
    ax = plt.subplot(111)

    # read in the terrace DataFrame
    terrace_df = H.read_terrace_csv(DataDirectory,fname_prefix)
    terrace_df = terrace_df[terrace_df['BaselineNode'] != -9999]

    # read in the mchi csv
    lp = H.ReadMChiSegCSV(DataDirectory,mchi_fname)
    lp = lp[lp['elevation'] != -9999]

    # get the distance from outlet along the baseline for each terrace pixels
    terrace_df = terrace_df.merge(lp, left_on = "BaselineNode", right_on = "node")
    flow_dist = terrace_df['flow_distance']/1000
    print(terrace_df)

	## Getting the extent of our dataset
    xmin = 0
    xmax = flow_dist.max()
    ymin = 0
    ymax = terrace_df["Elevation"].max()

    ## formatting the data in a meshgrid
    X,Y = np.meshgrid(np.linspace(0,xmax,num = prec),np.linspace(0,ymax, num = prec))
    positions = np.vstack([X.ravel(), Y.ravel()[::-1]]) # inverted Y to get the axis in the bottom left
    values = np.vstack([flow_dist, terrace_df['Elevation']])
    if len(values) == 0:
        print("You don't have any terraces, I'm going to quit now.")
    else:
        # get the kernel density estimation
        KDE = st.gaussian_kde(values, bw_method = bw_method)
        Z = np.reshape(KDE(positions).T,X.shape)

        # plot the density on the profile
        cmap = cm.gist_heat_r
        cmap.set_bad(alpha=0)
        cb = ax.imshow(Z, interpolation = "None",  extent=[xmin, xmax, ymin, ymax], cmap=cmap, aspect = "auto")

        # plot the main stem channel
        lp_mainstem = H.read_index_channel_csv(DataDirectory,fname_prefix)
        lp_mainstem = lp_mainstem[lp_mainstem['elevation'] != -9999]
        lp_mainstem = lp_mainstem.merge(lp, left_on="id", right_on="node")
        lp_flow_dist = lp_mainstem['flow_distance_y']/1000
        ax.plot(lp_flow_dist,lp_mainstem['elevation_y'],'k',lw=1, label='_nolegend_')

        # if present, plot the ages on the profile
        if ages:
            # read in the ages csv
            ages_df = pd.read_csv(DataDirectory+ages)
            upstream_dist = list(ages_df['upstream_dist'])
            elevation = list(ages_df['elevation'])
            ax.scatter(upstream_dist, elevation, s=8, c="w", edgecolors="k", label="$^{14}$C age (cal years B.P.)")
            ax.legend(loc='upper left', fontsize=8, numpoints=1)

        # set some plot lims
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        ax.set_xlabel('Flow distance (km)')
        ax.set_ylabel('Elevation (m)')

        # add a colourbar
        cbar = plt.colorbar(cb,cmap=cmap,orientation='vertical')
        cbar.set_label('Density')

        # save the figure
        plt.tight_layout()
        plt.savefig(T_directory+fname_prefix+'_terrace_plot_heat_map.png',format=FigFormat,dpi=300)
        plt.clf()

def MakeTerraceHeatMapNormalised(DataDirectory,fname_prefix, mchi_fname, prec=100, bw_method=0.03, FigFormat='png', ages=""):
    """
    Function to make a heat map of the terrace pixels using Gaussian KDE. Pixels are normalised based on
    elevation of closest channel pixel.
    see https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.gaussian_kde.html
    for more details.

    Args:
        DataDirectory(str): the data directory
        fname_prefix(str): prefix of your DEM
        mchi_fname(str): if you specified a junction then this will be different from the junction fname
        prec(int): the resolution for the KDE. Increase this to get a finer resolution, decrease for coarser.
        bw_method: the method for determining the bandwidth of the KDE.  This is apparently quite sensitive to this.
        Can either be "scott", "silverman" (where the bandwidth will be determined automatically), or a scalar. Default = 0.03
        FigFormat(str): figure format, default = png
        ages (str): Can pass in the name of a csv file with terrace ages which will be plotted on the profile. Must be in the same directory

    FJC 26/03/18
    """
    import scipy.stats as st

    # check if a directory exists for the chi plots. If not then make it.
    T_directory = DataDirectory+'terrace_plots/'
    if not os.path.isdir(T_directory):
        os.makedirs(T_directory)

    # make a figure
    fig = CreateFigure()
    ax = plt.subplot(111)
    #ax1 = plt.subplot(212)

    # read in the terrace DataFrame
    terrace_df = H.read_terrace_csv(DataDirectory,fname_prefix)
    terrace_df = terrace_df[terrace_df['BaselineNode'] != -9999]

    # read in the mchi  csv
    lp = H.ReadMChiSegCSV(DataDirectory,mchi_fname)
    lp = lp[lp['elevation'] != -9999]

    # get the distance from outlet along the baseline for each terrace pixels
    terrace_df = terrace_df.merge(lp, left_on = "BaselineNode", right_on = "node")
    flow_dist = terrace_df['flow_distance']/1000

	## Getting the extent of our dataset
    xmin = 0
    xmax = flow_dist.max()
    ymin = 0
    ymax = terrace_df["ChannelRelief"].max()

    ## formatting the data in a meshgrid
    X,Y = np.meshgrid(np.linspace(0,xmax,num = prec),np.linspace(0,ymax, num = prec))
    positions = np.vstack([X.ravel(), Y.ravel()[::-1]]) # inverted Y to get the axis in the bottom left
    values = np.vstack([flow_dist, terrace_df["ChannelRelief"]])
    KDE = st.gaussian_kde(values, bw_method = bw_method)
    Z = np.reshape(KDE(positions).T,X.shape)
    #Z = np.ma.masked_where(Z < 0.00000000001, Z)

    # try a 2d hist
    # h, fd_bins, elev_bins = np.histogram2d(flow_dist, terrace_df['Elevation'], bins=500)
    # h = h.T
    # h = np.ma.masked_where(h == 0, h)
    # X,Y = np.meshgrid(fd_bins, elev_bins)
    # ax.pcolormesh(X,Y,h, cmap="seismic")

    #
    cmap = cm.gist_heat_r
    cmap.set_bad(alpha=0)
    #norm=colors.LogNorm(vmin=0, vmax=Z.max(),cmap=cmap)
    cb = ax.imshow(Z, interpolation = "None",  extent=[xmin, xmax, ymin, ymax], cmap=cmap, aspect = "auto")
    #ax.pcolormesh(X,Y,Z, cmap="seismic")

    # set some plot lims
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    ax.set_xlabel('Flow distance (km)')
    ax.set_ylabel('Elevation above channel (m)')

    # add a colourbar
    cbar = plt.colorbar(cb,cmap=cmap,orientation='vertical')
    cbar.set_label('Density')

    plt.tight_layout()
    plt.savefig(T_directory+fname_prefix+'_terrace_plot_heat_map_norm.png',format=FigFormat,dpi=300)
    plt.clf()



#---------------------------------------------------------------------------------------------#
# RASTER PLOTS
# Functions to make raster plots of terrace attributes
#---------------------------------------------------------------------------------------------#

def MakeRasterPlotTerraceIDs(DataDirectory,fname_prefix, FigFormat='png', size_format='ESURF'):
    """
    This function makes a hillshade of the DEM with the terraces
    plotted onto it coloured by their ID

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the name of the DEM without extension.
        FigFormat (str): the figure format, default='png'
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).

    Returns:
        Raster plot of terrace IDs

    Author: FJC

    """
    from LSDMapFigure.PlottingRaster import BaseRaster
    from LSDMapFigure.PlottingRaster import MapFigure

    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        #fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        fig_width_inches=6.25
        #l_pad = -40
    elif size_format == "big":
        #fig = plt.figure(1, facecolor='white',figsize=(16,9))
        fig_width_inches=16
        #l_pad = -50
    else:
        fig_width_inches = 4.92126
        #fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
        #l_pad = -35

    # going to make the terrace plots - need to have bil extensions.
    print("I'm going to make the terrace raster plot. Your topographic data must be in ENVI bil format or I'll break!!")

    # get the rasters
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext
    TerraceIDName = fname_prefix+'_terrace_IDs'+raster_ext

    # get the terrace csv
    terraces = H.read_terrace_csv(DataDirectory,fname_prefix)
    terraceIDs = sorted(list(set(list(terraces.TerraceID))))
    xTerraces = []
    zTerraces = []
    yTerraces = []
    newIDs = []

    # loop through the terrace IDs and get the x, y, and z values
    for terraceID in terraceIDs:
        _terrace_subset = (terraces.TerraceID.values == terraceID)
        _x = terraces['DistAlongBaseline'].values[_terrace_subset]
        _y = terraces['DistToBaseline'].values[_terrace_subset]
        _z = terraces['Elevation'].values[_terrace_subset]
        _x_unique = sorted(list(set(list(_x))))
        _z_unique = []
        # Filter
        if len(_x) > 50 and len(_x_unique) > 1 and len(_x_unique) < 1000:
            #print np.max(np.diff(_x_unique))
            #if len(_x_unique) > 10 and len(_x_unique) < 1000 \
            #and :
            for _x_unique_i in _x_unique:
                #_y_unique_i = np.min(np.array(_y)[_x == _x_unique_i])
                #_z_unique.append(np.min(_z[_y == _y_unique_i]))
                _z_unique.append(np.min(_z[_x == _x_unique_i]))
            if np.mean(np.diff(_z_unique)/np.diff(_x_unique)) < 10:
                xTerraces.append(_x_unique)
                zTerraces.append(_z_unique)
                newIDs.append(terraceID)

    n_colours=len(newIDs)

    # create the map figure
    MF = MapFigure(HillshadeName, DataDirectory, coord_type='UTM_km', colourbar_location='right')
    # add the terrace drape
    terrace_cmap = plt.cm.rainbow
    #terrace_cmap = colours.cmap_discretize(n_colours,terrace_cmap)
    MF.add_drape_image(TerraceIDName, DataDirectory, colourmap = terrace_cmap, discrete_cmap=True, cbar_type=int, n_colours=n_colours, colorbarlabel="Terrace ID", alpha=0.8)

    ImageName = DataDirectory+fname_prefix+'_terrace_IDs_raster_plot.'+FigFormat
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 300) # Save the figure

def MakeRasterPlotTerraceElev(DataDirectory,fname_prefix, FigFormat='png', size_format='ESURF'):
    """
    This function makes a hillshade of the DEM with the terraces
    plotted onto it coloured by their elevation

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the name of the DEM without extension.
        FigFormat (str): the figure format, default='png'
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).

    Returns:
        Raster plot of terrace IDs

    Author: FJC

    """
    from LSDMapFigure.PlottingRaster import BaseRaster
    from LSDMapFigure.PlottingRaster import MapFigure

    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        #fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        fig_width_inches=6.25
        #l_pad = -40
    elif size_format == "big":
        #fig = plt.figure(1, facecolor='white',figsize=(16,9))
        fig_width_inches=16
        #l_pad = -50
    else:
        fig_width_inches = 4.92126
        #fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
        #l_pad = -35

    # going to make the terrace plots - need to have bil extensions.
    print("I'm going to make a raster plot of terrace elevations. Your topographic data must be in ENVI bil format or I'll break!!")

    # get the rasters
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext
    TerraceElevName = fname_prefix+'_terrace_relief_final'+raster_ext

    # get the terrace csv
    terraces = H.read_terrace_csv(DataDirectory,fname_prefix)
    terraceIDs = sorted(list(set(list(terraces.TerraceID))))
    xTerraces = []
    zTerraces = []
    yTerraces = []
    newIDs = []

    # loop through the terrace IDs and get the x, y, and z values
    for terraceID in terraceIDs:
        _terrace_subset = (terraces.TerraceID.values == terraceID)
        _x = terraces['DistAlongBaseline'].values[_terrace_subset]
        _y = terraces['DistToBaseline'].values[_terrace_subset]
        _z = terraces['Elevation'].values[_terrace_subset]
        _x_unique = sorted(list(set(list(_x))))
        _z_unique = []
        # Filter
        if len(_x) > 50 and len(_x_unique) > 1 and len(_x_unique) < 1000:
            #print np.max(np.diff(_x_unique))
            #if len(_x_unique) > 10 and len(_x_unique) < 1000 \
            #and :
            for _x_unique_i in _x_unique:
                #_y_unique_i = np.min(np.array(_y)[_x == _x_unique_i])
                #_z_unique.append(np.min(_z[_y == _y_unique_i]))
                _z_unique.append(np.min(_z[_x == _x_unique_i]))
            if np.mean(np.diff(_z_unique)/np.diff(_x_unique)) < 10:
                xTerraces.append(_x_unique)
                zTerraces.append(_z_unique)
                newIDs.append(terraceID)

    n_colours=len(newIDs)

    # create the map figure
    MF = MapFigure(HillshadeName, DataDirectory, coord_type='UTM_km', colourbar_location='right')
    # add the terrace drape
    terrace_cmap = plt.cm.Reds
    #terrace_cmap = colours.cmap_discretize(n_colours,terrace_cmap)
    MF.add_drape_image(TerraceElevName, DataDirectory, colourmap = terrace_cmap, colorbarlabel="Elevation above channel (m)", alpha=0.8)

    ImageName = DataDirectory+fname_prefix+'_terrace_elev_raster_plot.'+FigFormat
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 300) # Save the figure

def MakeRasterPlotTerraceDips(DataDirectory,fname_prefix,min_size=5000,FigFormat='png',size_format='ESURF'):
    """
    This function makes a raster plot of terrace locations with arrows showing the terrace
    dip and dip directions.
    Dip and dip direction are calculated by fitting a plane to each terrace using least-squares
    regression.

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the name of the DEM without extension.
        min_size (int): minimum number of pixels for a terrace, smaller ones will be removed
        FigFormat (str): the figure format, default='png'
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).

    Returns:
        plot of terrace locations and dip/dip directions

    Author: FJC

    """
    from LSDMapFigure.PlottingRaster import BaseRaster
    from LSDMapFigure.PlottingRaster import MapFigure

    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        #fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        fig_width_inches=6.25
        #l_pad = -40
    elif size_format == "big":
        #fig = plt.figure(1, facecolor='white',figsize=(16,9))
        fig_width_inches=16
        #l_pad = -50
    else:
        fig_width_inches = 4.92126
        #fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
        #l_pad = -35

    # going to make the terrace plots - need to have bil extensions.
    print("I'm going to make a raster plot of terrace elevations. Your topographic data must be in ENVI bil format or I'll break!!")

    # get the rasters
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext
    TerraceElevName = fname_prefix+'_terrace_relief_final'+raster_ext

    # get the terrace csv
    terraces = H.read_terrace_csv(DataDirectory,fname_prefix)
    filter_terraces(terraces)

    # get the terrace IDs
    terraceIDs = terraces.TerraceID.unique()
    n_colours=len(terraceIDs)

    # get the terrace dip and dip dirs
    terrace_dips = get_terrace_dip_and_dipdir(terraces)

    # create the map figure
    MF = MapFigure(HillshadeName, DataDirectory, coord_type='UTM_km', colourbar_location='right')
    # add the terrace drape
    terrace_cmap = plt.cm.Reds
    #terrace_cmap = colours.cmap_discretize(n_colours,terrace_cmap)
    MF.add_drape_image(TerraceElevName, DataDirectory, colourmap = terrace_cmap, colorbarlabel="Elevation above channel (m)", alpha=0.8)

    # add arrows oriented in the direction of dip. We might want to colour these by the dip angle?
    # MF.add_arrows_from_points(terrace_dips,azimuth_header='dip_azimuth', arrow_length=100)
    MF.add_strike_and_dip_symbols(terrace_dips,symbol_length=100,linewidth=0.5)


    ImageName = DataDirectory+fname_prefix+'_terrace_dips_raster_plot.'+FigFormat
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 300) # Save the figure
