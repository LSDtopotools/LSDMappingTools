# -*- coding: utf-8 -*-
"""
This file contains are a series of helper functions that interact with
LSDTopoTools spatial data.

Written by Simon M. Mudd and Fiona J Clubb at the University of Edinburgh
Released under GPL3

22/06/2017
"""

import os
import pandas as pd
import fiona
from shapely.geometry import shape, Polygon, Point, LineString

#==============================================================================
def MapFigureSizer(figure_width_inches,aspect_ratio, cbar_loc = "None", title = "None",
                   cbar_width = 0.2,
                   cbar_text_width = 0.4,
                   cbar_padding = 0.1,
                   cbar_fraction = 1,
                   whitespace_padding = 0.2,
                   map_text_width = 0.65,
                   map_text_height = 0.45,
                   title_height=0.2):
    """This function takes a string size argument and calculates the size of the
    various components of a plot based on a map image making up the centre of the
    figure.

    We use inches because bloody yanks wrote matplotlib and figures in matplotlib use inches.
    Luckily we do not have to calculate bushels or furlongs of pixels somewhere,
    and the inches stupidity is somewhat limited.

    Args:
        figure_width_inches (flt): The figure width in inches
        aspect_ratio (flt): The width to height ratio of the data
        cbar_loc (string): the location of the colourbar, either "left", "right", "top",
        "bottom", or "none"
        title (bool): if true then adjust the height of the figure to make space for a title
        cbar_width (flt): the width of the colorbar
        text_padding (list): the padding around the map from the whitespace and the axis/tick labels

    Author: SMM
    """


    # The text padding list holds the height and width, in inches of
    # [0] = left/right tick marks+tick labels
    # [1] = left/right text (i.e., and exis label)
    # [2] = top/bottom tick marks+tick labels
    # [3] = top/bottom text (i.e., and exis label)
    # [4] = whitespace at the edge of the figure + between map and colourbar
    # [5] = colorbar text (e.g. tick labels+axis label)


    # This gets returned, we use it to make the figure.

    #======================================================

    # Now we need to figure out where the axis are. Sadly this requires
    # a load of tedious conditional statments about the location of the axes

    #Will's changes:
    #1/ Make space for the colourbar label on the left
    #2/ Double the cbar_padding to leave space for the colourbar values on the right.
    # NB: this should later be a function of font size
    #3/ Changed rotation of colourbar text to 90 and the labelpad to -75 in PlottingRaster.py

    if cbar_loc == "left":
        cbar_left_inches = whitespace_padding + cbar_text_width
        #cbar_left_inches = whitespace_padding
        map_left_inches = cbar_left_inches+cbar_width+map_text_width + 2*cbar_padding
        #map_left_inches = cbar_left_inches+cbar_width+cbar_text_width+map_text_width+cbar_padding
        map_width_inches = figure_width_inches-map_left_inches-whitespace_padding
        map_height_inches = map_width_inches/aspect_ratio



        map_bottom_inches = whitespace_padding+map_text_height
        cbar_bottom_inches = map_bottom_inches
        figure_height_inches = map_bottom_inches+map_height_inches+whitespace_padding
        if title != "None":
            # add some space for a title if needed. At the moment this is hard coded but
            # we might want to change this for the font size.
            figure_height_inches = figure_height_inches+title_height

        fig_size_inches = [figure_width_inches,figure_height_inches]

        print("cbar_left: "+str(cbar_left_inches)+" map left: "+str(map_left_inches))
        print("cbar_bottom: "+str(cbar_bottom_inches)+ " map bottom: "+str(map_bottom_inches))

        map_axes = [map_left_inches/figure_width_inches,
                    map_bottom_inches/figure_height_inches,
                    map_width_inches/figure_width_inches,
                    map_height_inches/figure_height_inches]
        cbar_axes = [cbar_left_inches/figure_width_inches,
                    map_bottom_inches/figure_height_inches,
                    cbar_width/figure_width_inches,
                    cbar_fraction*(map_height_inches/figure_height_inches)]

    elif cbar_loc == "right":

        map_left_inches = whitespace_padding+map_text_width
        cbar_left_inches= figure_width_inches-whitespace_padding-cbar_width-cbar_text_width
        map_right_inches = cbar_left_inches-cbar_padding

        map_width_inches = map_right_inches-map_left_inches
        map_height_inches = map_width_inches/aspect_ratio

        map_bottom_inches = whitespace_padding+map_text_height
        cbar_bottom_inches = map_bottom_inches
        figure_height_inches = map_bottom_inches+map_height_inches+whitespace_padding
        if title != "None":
            # add some space for a title if needed. At the moment this is hard coded but
            # we might want to change this for the font size.
            figure_height_inches = figure_height_inches+title_height

        fig_size_inches = [figure_width_inches,figure_height_inches]

        print("cbar_left: "+str(cbar_left_inches)+" map left: "+str(map_left_inches))
        print("cbar_bottom: "+str(cbar_bottom_inches)+ " map bottom: "+str(map_bottom_inches))


        map_axes = [map_left_inches/figure_width_inches,
                    map_bottom_inches/figure_height_inches,
                    map_width_inches/figure_width_inches,
                    map_height_inches/figure_height_inches]
        cbar_axes = [cbar_left_inches/figure_width_inches,
                    map_bottom_inches/figure_height_inches,
                    cbar_width/figure_width_inches,
                    cbar_fraction*(map_height_inches/figure_height_inches)]

    elif cbar_loc == "top":
        print("I am placing the colourbar on the top")

        map_left_inches = whitespace_padding+map_text_width
        map_right_inches = figure_width_inches-whitespace_padding
        map_width_inches = map_right_inches-map_left_inches
        map_height_inches = map_width_inches/aspect_ratio

        cbar_left_inches= map_left_inches

        map_bottom_inches = whitespace_padding+map_text_height
        cbar_bottom_inches = map_bottom_inches+map_height_inches+cbar_padding+cbar_text_width

        figure_height_inches = cbar_bottom_inches+cbar_width+whitespace_padding
        if title != "None":
            # add some space for a title if needed. At the moment this is hard coded but
            # we might want to change this for the font size.
            title_height = 0.5
            figure_height_inches = figure_height_inches+title_height

        fig_size_inches = [figure_width_inches,figure_height_inches]

        print("cbar_left: "+str(cbar_left_inches)+" map left: "+str(map_left_inches))
        print("cbar_bottom: "+str(cbar_bottom_inches)+ " map bottom: "+str(map_bottom_inches))


        map_axes = [map_left_inches/figure_width_inches,
                    map_bottom_inches/figure_height_inches,
                    map_width_inches/figure_width_inches,
                    map_height_inches/figure_height_inches]
        cbar_axes = [cbar_left_inches/figure_width_inches,
                    cbar_bottom_inches/figure_height_inches,
                    cbar_fraction*(map_width_inches/figure_width_inches),
                    cbar_width/figure_height_inches]

    elif cbar_loc == "bottom":
        print("I am placing the colourbar on the bottom")

        map_left_inches = whitespace_padding+map_text_width
        map_right_inches = figure_width_inches-whitespace_padding
        map_width_inches = map_right_inches-map_left_inches
        map_height_inches = map_width_inches/aspect_ratio

        cbar_left_inches= map_left_inches

        cbar_bottom_inches = whitespace_padding+cbar_text_width
        map_bottom_inches = cbar_bottom_inches+cbar_width+cbar_padding+map_text_height

        whitespace_padding+map_text_height

        figure_height_inches = map_bottom_inches+map_height_inches+whitespace_padding
        if title != "None":
            # add some space for a title if needed. At the moment this is hard coded but
            # we might want to change this for the font size.
            figure_height_inches = figure_height_inches+title_height

        fig_size_inches = [figure_width_inches,figure_height_inches]

        print("cbar_left: "+str(cbar_left_inches)+" map left: "+str(map_left_inches))
        print("cbar_bottom: "+str(cbar_bottom_inches)+ " map bottom: "+str(map_bottom_inches))


        map_axes = [map_left_inches/figure_width_inches,
                    map_bottom_inches/figure_height_inches,
                    map_width_inches/figure_width_inches,
                    map_height_inches/figure_height_inches]
        cbar_axes = [cbar_left_inches/figure_width_inches,
                    cbar_bottom_inches/figure_height_inches,
                    cbar_fraction*(map_width_inches/figure_width_inches),
                    cbar_width/figure_height_inches]

    else:
        print("No colourbar")

        map_left_inches = whitespace_padding+map_text_width
        map_right_inches = figure_width_inches-whitespace_padding
        map_width_inches = map_right_inches-map_left_inches
        map_height_inches = map_width_inches/aspect_ratio

        map_bottom_inches = whitespace_padding+map_text_height

        figure_height_inches = map_bottom_inches+map_height_inches+whitespace_padding
        if title != "None":
            # add some space for a title if needed. At the moment this is hard coded but
            # we might want to change this for the font size.
            figure_height_inches = figure_height_inches+title_height

        fig_size_inches = [figure_width_inches,figure_height_inches]

        map_axes = [map_left_inches/figure_width_inches,
                    map_bottom_inches/figure_height_inches,
                    map_width_inches/figure_width_inches,
                    map_height_inches/figure_height_inches]
        cbar_axes = None


    print("The figure size is: ")
    print(fig_size_inches)
    print("Map axes are:")
    print(map_axes)
    print("cbar_axes are:")
    print(cbar_axes)
    return fig_size_inches, map_axes, cbar_axes


#==============================================================================
# Formats ticks for an imshow plot in UTM
# Filename is the name of the file with full path
# x_max, x_min, y_max, y_min are the extent of the plotting area (NOT the DEM)
# n_target ticks are the number of ticks for plotting
#------------------------------------------------------------------------------
def GetTicksForUTMNoInversion(FileName,x_max,x_min,y_max,y_min,n_target_tics):
    """This fuction is used to set tick locations for UTM maps. It tries to optimise the spacing of these ticks.

    Args:
        x_min (float): The minimum value on the x axis (in metres).
        x_max (float): The maximum value on the x axis (in metres).
        y_min (float): The minimum value on the y axis (in metres).
        y_max (float): The maximum value on the y axis (in metres).
        n_target_ticks (int): The number of ticks you want on the axis (this is optimised so you may not get exactly this number)

    Returns:
        new_xlocs (float list): List of locations of the ticks in metres.
        new_x_labels (str list): List of strings for ticks, will be location in kilometres.
        new_ylocs (float list): List of locations of the ticks in metres.
        new_y_labels (str list): List of strings for ticks, will be location in kilometres.

    Author: SMM
    """

    CellSize,XMin,XMax,YMin,YMax = LSDMap_IO.GetUTMMaxMin(FileName)
    NDV, xsize, ysize, GeoT, Projection, DataType = LSDMap_IO.GetGeoInfo(FileName)

    #print("Getting ticks. YMin: "+str(YMin)+" and YMax: "+str(YMax))

    xmax_UTM = XMax
    xmin_UTM = XMin

    ymax_UTM = YMax
    ymin_UTM = YMin

    dy_fig = ymax_UTM-ymin_UTM
    dx_fig = xmax_UTM-xmin_UTM

    dx_spacing = dx_fig/n_target_tics
    dy_spacing = dy_fig/n_target_tics

    if (dx_spacing>dy_spacing):
        dy_spacing = dx_spacing

    str_dy = str(dy_spacing)
    str_dy = str_dy.split('.')[0]
    n_digits = str_dy.__len__()
    nd = int(n_digits)

    first_digit = float(str_dy[0])

    dy_spacing_rounded = first_digit*pow(10,(nd-1))

    str_xmin = str(xmin_UTM)
    str_ymin = str(ymin_UTM)
    str_xmin = str_xmin.split('.')[0]
    str_ymin = str_ymin.split('.')[0]
    xmin_UTM = float(str_xmin)
    ymin_UTM = float(str_ymin)

    print("minimum values are x: "+str(xmin_UTM)+ " and y: "+str(ymin_UTM))

    n_digx = str_xmin.__len__()
    n_digy = str_ymin.__len__()


    if (n_digx-nd+1) >= 1:
        front_x = str_xmin[:(n_digx-nd+1)]
    else:
        front_x = str_xmin

    if (n_digy-nd+1) >= 1:
        front_y = str_ymin[:(n_digy-nd+1)]
    else:
        front_y = str_ymin


    round_xmin = float(front_x)*pow(10,nd-1)
    round_ymin = float(front_y)*pow(10,nd-1)

    #print("UTM y in: "+str(ymin_UTM)+" and rounded min: "+str(round_ymin))


    # now we need to figure out where the xllocs and ylocs are
    xUTMlocs = np.zeros(2*n_target_tics)
    yUTMlocs = np.zeros(2*n_target_tics)
    xlocs = np.zeros(2*n_target_tics)
    ylocs = np.zeros(2*n_target_tics)

    new_x_labels = []
    new_y_labels = []

    for i in range(0,2*n_target_tics):

        # Note we use dy spacing here in both x and y directions since we want
        # Ticks spaced the same in each direction
        xUTMlocs[i] = round_xmin+(i)*dy_spacing_rounded
        yUTMlocs[i] = round_ymin+(i)*dy_spacing_rounded
        xlocs[i] = xUTMlocs[i]
        ylocs[i] = yUTMlocs[i]

        new_x_labels.append( str(xUTMlocs[i]).split(".")[0] )
        new_y_labels.append( str(yUTMlocs[i]).split(".")[0] )


    new_xlocs = []
    new_xUTMlocs = []
    x_labels = []

    # Now loop through these to get rid of those not in range
    for index,xloc in enumerate(xlocs):
        if (xloc < XMax and xloc > XMin):
            new_xlocs.append(xloc)
            new_xUTMlocs.append(xUTMlocs[index])
            x_labels.append(new_x_labels[index])

    new_ylocs = []
    new_yUTMlocs = []
    y_labels = []



    for index,yloc in enumerate(ylocs):
        if (yloc < YMax and yloc > YMin):
            new_ylocs.append(yloc)
            new_yUTMlocs.append(yUTMlocs[index])
            y_labels.append(new_y_labels[index])


    #return xlocs,ylocs,new_x_labels,new_y_labels
    return new_xlocs,new_ylocs,x_labels,y_labels
#==============================================================================


#=============================================================================
# CSV READERS
# Read in the csv files to pandas dataframes
#=============================================================================
def ReadBaselevelKeysCSV(DataDirectory, fname_prefix):
    """
    This function reads in the file with the suffix '__BaselevelKeys.csv'
    to a pandas dataframe

    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix

    Returns:
        pandas dataframe with the csv file

    Author: FJC
    """
    # get the csv filename
    baselevel_suffix = '_BaselevelKeys.csv'
    fname = fname_prefix+baselevel_suffix
    # read in the dataframe using pandas
    df = pd.read_csv(DataDirectory+fname)

    return df

def ReadSourceKeysCSV(DataDirectory, fname_prefix):
    """
    This function reads in the file with the suffix '_SourceKeys.csv'
    to a pandas dataframe

    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix

    Returns:
        pandas dataframe with the csv file

    Author: FJC
    """
    
    # get the csv filename
    source_keys_suffix = '_SourceKeys.csv'
    fname = fname_prefix+source_keys_suffix
    # read in the dataframe using pandas
    df = pd.read_csv(DataDirectory+fname)

    return df

def ReadBasinInfoCSV(DataDirectory, fname_prefix):
    """
    This function reads in the file with the suffix '_AllBasinsInfo.csv'
    to a pandas dataframe

    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix

    Returns:
        pandas dataframe with the csv file

    Author: FJC
    """
    
    # get the csv filename
    basin_suffix = '_AllBasinsInfo.csv'
    fname = fname_prefix+basin_suffix
    # read in the dataframe using pandas
    df = pd.read_csv(DataDirectory+fname)

    return df

def ReadFullStatsCSV(DataDirectory, fname_prefix, m_over_n):
    """
    This function reads in the file with the suffix '_fullstats.csv'
    to a pandas dataframe. Must specify the m/n value as an argument

    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix
        m_over_n: the m/n value

    Returns:
        pandas dataframe with the csv file

    Author: FJC
    """
    
    # get the csv filename
    fullstats_suffix = '_movernstats_%s_fullstats.csv' % m_over_n
    fname = fname_prefix+fullstats_suffix
    # read in the dataframe using pandas
    df = pd.read_csv(DataDirectory+fname)

    return df

def ReadChiProfileCSV(DataDirectory, fname_prefix):
    """
    This function reads in the file with the suffix '_movern.csv', which
    contains the data for the full chi profiles, to a pandas dataframe.

    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix

    Returns:
        pandas dataframe with the csv file

    Author: FJC
    """
    
    import os
    
    # get the csv filename
    profile_suffix = '_burned_movern.csv'
    fname = fname_prefix+profile_suffix

    # check if this exists. if not then use the burned one.
    if not os.path.isfile(DataDirectory+fname):
        print ("Reading the csv...")
        fname = fname_prefix+'_movern.csv'
    else:
        print("Reading the burned csv...")

    # read in the dataframe using pandas
    df = pd.read_csv(DataDirectory+fname)
    return df

def ReadBasinStatsCSV(DataDirectory, fname_prefix):
    """
    This function reads in the file with the suffix '_disorder_basinstats.csv'
    to a pandas dataframe

    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix

    Returns:
        pandas dataframe with the csv file

    Author: FJC
    """
    
    # get the csv filename
    basin_stats_suffix = '_disorder_basinstats.csv'
    fname = fname_prefix+basin_stats_suffix
    # read in the dataframe using pandas
    df = pd.read_csv(DataDirectory+fname)

    return df

def AppendBasinStatsCSVs(DataDirectory, FilenamePrefix):
    """
    This function reads in the files with the prefic "basin"
    and the suffix '_movernstats_basinstats.csv'
    to a pandas dataframe, for use with parallel analysis

    Args:
        DataDirectory: the data directory

    Returns:
        pandas dataframe with the csv file

    Author: FJC, MDH
    """

    # get the csv filename
    basin_stats_suffix = '_disorder_basinstats.csv'

    # decclare empty data frame
    MasterDF = pd.DataFrame()

    # get the list of basins as a dict
    basin_dict = MapBasinsToKeysFromJunctionList(DataDirectory,FilenamePrefix)

    # loop through and get each basin csv
    for outlet_jn, basin_key in basin_dict.iteritems():
        print ((outlet_jn, basin_key))
        this_fname = "basin"+str(outlet_jn)+basin_stats_suffix
        df = pd.read_csv(DataDirectory+this_fname)
        df = df[df['basin_key'] == 0]
        df['basin_key'] = basin_key
        MasterDF = MasterDF.append(df, ignore_index = True)

    return MasterDF

def ReadBasinStatsPointCSV(DataDirectory, fname_prefix):
    """
    This function reads in the file with the suffix '_point_movernstats_basinstats.csv'
    to a pandas dataframe

    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix

    Returns:
        pandas dataframe with the csv file

    Author: FJC
    """
    # get the csv filename
    csv_suffix = '_point_movernstats_basinstats.csv'
    fname = fname_prefix+csv_suffix
    # read in the dataframe using pandas
    df = pd.read_csv(DataDirectory+fname)

    return df

def ReadChainCSV(DataDirectory, fname_prefix, basin_key):
    """
    This function reads in the file with the suffix '_BasinX_chain.csv'
    to a pandas dataframe, where X is the basin key

    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix
        basin_key: the basin key

    Returns:
        pandas dataframe with the csv file

    Author: FJC
    """
    # get the csv filename
    chain_suffix = '_Basin%s_chain.csv' %str(basin_key)
    fname = fname_prefix+chain_suffix
    # read in the dataframe using pandas
    df = pd.read_csv(DataDirectory+fname)

    return df

def ReadMCPointsCSV(DataDirectory, fname_prefix):
    """
    This function reads in the file with the suffix
    '_MCpoint__points_MC_basinstats.csv'

    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix

    Returns:
        pandas dataframe with the csv file

    Author: FJC
    """
    # get the csv filename
    mc_points_suffix = '_MCpoint_points_MC_basinstats.csv'
    fname = fname_prefix+mc_points_suffix
    # read in the dataframe using pandas
    df = pd.read_csv(DataDirectory+fname)

    return df


def ReadMChiSegCSV(DataDirectory, fname_prefix, type = "Normal"):
    """
    This function reads in the file with the suffix
    '_MChiSegmented.csv'
    This file holds the MCHI segmented data

    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix

    Returns:
        pandas dataframe with the csv file

    Author: BG
    """
    # get the csv filename depending of what you need

    # The "knickpoint" type is a special M_Chi file generated with the exact degmented elevation
    if(type == "Normal"):
        suffix = '_MChiSegmented.csv'
    elif(type == "knickpoint"):
        suffix = "_ksnkp_mchi.csv"
    fname = fname_prefix+suffix
    # read in the dataframe using pandas
    df = pd.read_csv(DataDirectory+fname)

    # Getting rid of NoData

    df = df[df["chi"] >= 0]

    return df

def ReadDisorderCSV(DataDirectory, fname_prefix):
    """
    Function to read in the CSV from the chi disorder
    analysis
    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix

    Returns:
        pandas dataframe with the csv file

    Author: FJC
    """
    # get the csv filename
    mc_points_suffix = '_disorder_movernstats_disorder_basinstats.csv'
    fname = fname_prefix+mc_points_suffix
    # read in the dataframe using pandas
    df = pd.read_csv(DataDirectory+fname)

    return df

def ReadDisorderUncertCSV(DataDirectory, fname_prefix):
    """
    Function to read in the CSV from the chi disorder
    analysis
    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix

    Returns:
        pandas dataframe with the csv file

    Author: FJC
    """
    # get the csv filename
    mc_points_suffix = '_fullstats_disorder_uncert.csv'
    fname = fname_prefix+mc_points_suffix
    # read in the dataframe using pandas
    df = pd.read_csv(DataDirectory+fname)

    return df

def AppendDisorderCSV(DataDirectory, FilenamePrefix):
    """
    This function reads in a series of disorder csvs and appends them together
    into one function for plotting

    Args:
        DataDirectory (str): the data DataDirectory

    Returns:
        pandas dataframe with the appended fullstats csvs

    Author: FJC
    """
    # get the csv filename
    csv_suffix =  '_fullstats_disorder_uncert.csv'

    MasterDF = pd.DataFrame()
    basin_dict = MapBasinsToKeysFromJunctionList(DataDirectory, FilenamePrefix)

    # loop through and get each basin csv
    for outlet_jn, basin_key in basin_dict.iteritems():
        this_fname = "basin"+str(outlet_jn)+csv_suffix
        # append to master DF and change the basin key and the junction
        df = pd.read_csv(DataDirectory+this_fname)
        df = df[df['basin_key'] == 0]
        df['basin_key'] = basin_key
        MasterDF = MasterDF.append(df, ignore_index = True)

    return MasterDF

def readSKKPstats(DataDirectory, fname_prefix):
    """
    This function reads in the file with the suffix
    '_KsnKn.csv'
    This file holds the MCHI segmented data

    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix

    Returns:
        pandas dataframe with the csv file

    Author: BG
    """
    # get the csv filename
    suffix = '_ksnkp_SK.csv'
    fname = fname_prefix+suffix
    # read in the dataframe using pandas
    df = pd.read_csv(DataDirectory+fname)

    return df

def ReadKnickpointCSV(DataDirectory, fname_prefix, ftype = "normal"):
    """
    This function reads in the file with the suffix
    '_KsnKn.csv'
    This file holds the MCHI segmented data

    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix

    Returns:
        pandas dataframe with the csv file

    Author: BG
    """
    # get the csv filename
    if(ftype == "raw"):
        suffix = '_ksnkp_raw.csv'
    else:
        suffix = '_ksnkp.csv'
    fname = fname_prefix+suffix
    # read in the dataframe using pandas
    df = pd.read_csv(DataDirectory+fname)

    return df

def ReadKnickzoneCSV(DataDirectory, fname_prefix):
    """
    This function reads in the file with the suffix
    '_KsnKn.csv'
    This file holds the MCHI segmented data

    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix

    Returns:
        pandas dataframe with the csv file

    Author: BG
    """
    # get the csv filename
    suffix = '_KsnKz.csv'
    fname = fname_prefix+suffix
    # read in the dataframe using pandas
    df = pd.read_csv(DataDirectory+fname)

    return df

def ReadChiResidualsCSVs(DataDirectory, fname_prefix):
    """
    This function reads in the 3 CSV files for the residuals analysis
    They have the format:
        "_residual_movernstats_movern_residuals_median.csv"
        "_residual_movernstats_movern_residuals_Q1.csv"
        "_residual_movernstats_movern_residuals_Q3.csv"

    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix

    Returns:
        list of pandas dataframes with the csv files. List format is:
        0 = medians
        1 = 1st quartile
        2 = 3rd quartile

    Author: FJC
    """
    # get the csv filename
    fnames = ["_residual_movernstats_movern_residuals_median.csv","_residual_movernstats_movern_residuals_Q1.csv","_residual_movernstats_movern_residuals_Q3.csv"]
    #print "The fnames are: ", fnames
    dfs = []
    for f in fnames:
        fname = fname_prefix+f
        dfs.append(pd.read_csv(DataDirectory+fname))

    return dfs

def ReadRawSAData(DataDirectory, fname_prefix):
    """
    This function reads in the raw SA data to a pandas dataframe

    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix

    Returns:
        pandas dataframe with the raw SA data

    Author: FJC
    """
    
    # get the csv filename
    fname_suffix = "_SAvertical.csv"
    fname = fname_prefix+fname_suffix
    df = pd.read_csv(DataDirectory+fname)

    return df

def AppendRawSAData(DataDirectory, FilenamePrefix):
    """
    This function reads in the raw SA data to a pandas dataframe
    from multiple CSV files with the filename prefix "basin"
    For use with parallelised methods

    Args:
        DataDirectory: the data directory

    Returns:
        pandas dataframe with the raw SA data

    Author: MDH, FJC
    """
    
    # get the csv filename
    csv_suffix = "_SAvertical.csv"

    MasterDF = pd.DataFrame()
    basin_dict = MapBasinsToKeysFromJunctionList(DataDirectory, FilenamePrefix)

    # loop through and get each basin csv
    for outlet_jn, basin_key in basin_dict.iteritems():
        this_fname = "basin"+str(outlet_jn)+csv_suffix
        df = pd.read_csv(DataDirectory+this_fname)
        df = df[df['basin_key'] == 0]
        df['basin_key'] = basin_key
        MasterDF = MasterDF.append(df, ignore_index = True)

    return MasterDF

def ReadSegmentedSAData(DataDirectory, fname_prefix):
    """
    This function reads in the segmented SA data to a pandas dataframe

    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix

    Returns:
        pandas dataframe with the segmented SA data

    Author: FJC
    """
    
    # get the csv filename
    fname_suffix = "_SAsegmented.csv"
    fname = fname_prefix+fname_suffix
    df = pd.read_csv(DataDirectory+fname)

    return df

def ReadBinnedSAData(DataDirectory, fname_prefix):
    """
    This function reads in the binned SA data to a pandas dataframe
    csv with the suffix "_SAbinned.csv"

    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix

    Returns:
        pandas dataframe with the segmented SA data

    Author: FJC
    """
    
    # get the csv filename
    fname_suffix = "_SAbinned.csv"
    fname = fname_prefix+fname_suffix
    df = pd.read_csv(DataDirectory+fname)

    return df


def ReadMOverNSummaryCSV(DataDirectory, fname_prefix):
    """
    This function reads in the summary csv with the best fit movern info
    to a pandas dataframe

    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix

    Returns:
        pandas dataframe with the segmented SA data

    Author: FJC
    """
    # get the csv filename
    fname_suffix = "_movern_summary.csv"
    fname = fname_prefix+fname_suffix
    df = pd.read_csv(DataDirectory+fname)

    return df

def ReadChannelNetworkCSV(DataDirectory, fname_prefix):
    """
    This function reads in the channel network csv to a df

    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix

    Returns:
        pandas dataframe with the channel network

    Author: FJC
    """
    # get the csv filename
    fname_suffix = "_CN.csv"
    fname = fname_prefix+fname_suffix
    df = pd.read_csv(DataDirectory+fname)

    return df

def ReadChiDataMapCSV(DataDirectory, fname_prefix):
    """
    This function reads in the chi data map csv to a df

    Args:
        DataDirectory: the data directory
        fname_prefix: the file name prefix

    Returns:
        pandas dataframe with the chi map

    Author: FJC
    """
    # get the csv filename
    fname_suffix = "_chi_data_map.csv"
    fname = fname_prefix+fname_suffix
    df = pd.read_csv(DataDirectory+fname)

    return df


#--------------------------------------------------------------------------------#
# Terraces
#--------------------------------------------------------------------------------#

def read_terrace_csv(DataDirectory,fname_prefix):
    """
    This function reads in the csv file with the extension "_terrace_info.csv"
    and returns it as a pandas dataframe

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the name of the DEM

    Returns:
        pandas dataframe with the terrace info

    Author: FJC
    """
    csv_suffix = '_terrace_info.csv'
    fname = DataDirectory+fname_prefix+csv_suffix

    df = pd.read_csv(fname)

    return df

def read_channel_csv(DataDirectory,fname_prefix):
    """
    This function reads in the csv file with the extension "_baseline_channel_info.csv"
    and returns it as a pandas dataframe

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the name of the DEM

    Returns:
        pandas dataframe with the channel info

    Author: FJC
    """
    csv_suffix = '_baseline_channel_info.csv'
    fname = DataDirectory+fname_prefix+csv_suffix

    df = pd.read_csv(fname)

    return df

def read_index_channel_csv(DataDirectory,fname_prefix):
    """
    This function reads in the csv file with the extension "_index_chan.csv"
    and returns it as a pandas dataframe

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the name of the DEM

    Returns:
        pandas dataframe with the channel info

    Author: FJC
    """
    csv_suffix = '_index_chan.csv'
    fname = DataDirectory+fname_prefix+csv_suffix

    df = pd.read_csv(fname)

    return df

def read_terrace_shapefile(DataDirectory, shapefile_name):
    """
    This function reads in a shapefile of digitised terraces
    using shapely and fiona

    Args:
        DataDirectory (str): the data directory
        shapefile_name (str): the name of the shapefile

    Returns: shapely polygons with terraces

    Author: FJC
    """
    Polygons = {}
    with fiona.open(DataDirectory+shapefile_name, 'r') as input:
        for f in input:
            this_shape = Polygon(shape(f['geometry']))
            this_id = f['properties']['id']
            Polygons[this_id] = this_shape

    return Polygons

def read_terrace_centrelines(DataDirectory, shapefile_name):
    """
    This function reads in a shapefile of terrace centrelines
    using shapely and fiona

    Args:
        DataDirectory (str): the data directory
        shapefile_name (str): the name of the shapefile

    Returns: shapely polygons with terraces

    Author: FJC
    """
    Lines = {}
    with fiona.open(DataDirectory+shapefile_name, 'r') as input:
        for f in input:
            this_line = LineString(shape(f['geometry']))
            this_id = f['properties']['id']
            Lines[this_id] = this_line
    return Lines

def ReadModelCSV(DataDirectory, Base_file):
    """
    This function reads in the csv file from the model run to a pandas dataframe

    Args:
        DataDirectory (str): the data directory
        Base_file (str): the base file prefix

    Returns:
        pandas dataframe with the csv file info

    Author: FJC
    """
    # get the csv filename
    csv_suffix = '_model_info.csv'

    fname = Base_file+csv_suffix
    # read in the dataframe using pandas
    df = pd.read_csv(DataDirectory+fname)
    return df

#-----------------------------------------------------------------------------#
# Drainage capture metrics
#-----------------------------------------------------------------------------#
def ReadPerimeterCSV(DataDirectory, fname_prefix):
    """
    This function reads in the csv file with the perimeter info

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the base file prefix

    Returns:
        pandas dataframe with the csv file info

    Author: FJC
    """
    csv_suffix = '_Perimeters.csv'
    df = pd.read_csv(DataDirectory+fname_prefix+csv_suffix)
    return df

#-----------------------------------------------------------------------------#
# Functions for appending csvs together for parallel basin running
# FJC 19/10/17
#-----------------------------------------------------------------------------#

def AppendBasinCSVs(DataDirectory, FilenamePrefix):
    """
    This function reads in a series of basin csv files and appends them together
    into one function for plotting

    Args:
        DataDirectory (str): the data DataDirectory

    Returns:
        pandas dataframe with the appended basin csvs

    Author: FJC
    """

    # get the csv filename
    csv_suffix = "_movernstats_basinstats.csv"

    MasterDF = pd.DataFrame()
    basin_dict = MapBasinsToKeysFromJunctionList(DataDirectory,FilenamePrefix)

    # loop through and get each basin csv
    for outlet_jn, basin_key in basin_dict.iteritems():
        this_fname = "basin"+str(outlet_jn)+csv_suffix
        # append to master DF and change the basin key and the junction
        df = pd.read_csv(DataDirectory+this_fname)
        df.ix[0, 'basin_key'] = basin_key
        df.ix[0, 'outlet_jn'] = outlet_jn
        MasterDF = MasterDF.append(df.iloc[0], ignore_index = True)

    return MasterDF

def AppendFullStatsCSVs(DataDirectory, m_over_n, FilenamePrefix):
    """
    This function reads in a series of full stats csvs and appends them together
    into one function for plotting

    Args:
        DataDirectory (str): the data DataDirectory
        m_over_n (float): the m/n value

    Returns:
        pandas dataframe with the appended fullstats csvs

    Author: FJC
    """
    # get the csv filename
    csv_suffix =  '_movernstats_%s_fullstats.csv' % m_over_n

    MasterDF = pd.DataFrame()
    basin_dict = MapBasinsToKeysFromJunctionList(DataDirectory, FilenamePrefix)

    # loop through and get each basin csv
    for outlet_jn, basin_key in basin_dict.iteritems():
        this_fname = "basin"+str(outlet_jn)+csv_suffix
        # append to master DF and change the basin key and the junction
        df = pd.read_csv(DataDirectory+this_fname)
        df = df[df['basin_key'] == 0]
        df['basin_key'] = basin_key
        MasterDF = MasterDF.append(df, ignore_index = True)

    return MasterDF

def ReadMovernCSV(DataDirectory, fname_prefix):
    """
    This function reads in a the movern csv with the suffix "_movern"

    Args:
        DataDirectory (str): the data DataDirectory
        fname_prefix
    Returns:
        pandas dataframe with the appended movern csvs

    Author: MDH
    """

    # get the csv filename
    csv_suffix = '_movern.csv'

    df = pd.read_csv(DataDirectory+fname_prefix+csv_suffix)

    return df

def AppendMovernCSV(DataDirectory, FilenamePrefix):
    """
    This function reads in a series of csvs with the suffix "_movern"
    and appends them together into one function for plotting

    Args:
        DataDirectory (str): the data DataDirectory

    Returns:
        pandas dataframe with the appended movern csvs

    Author: FJC
    """

    # get the csv filename
    csv_suffix =  '_movern.csv'

    MasterDF = pd.DataFrame()
    basin_dict = MapBasinsToKeysFromJunctionList(DataDirectory,FilenamePrefix)

    # loop through and get each basin csv
    for outlet_jn, basin_key in basin_dict.iteritems():
        this_fname = "basin"+str(outlet_jn)+csv_suffix
        # append to master DF and change the basin key and the junction
        df = pd.read_csv(DataDirectory+this_fname)
        df = df[df['basin_key'] == 0]
        df['basin_key'] = basin_key
        MasterDF = MasterDF.append(df, ignore_index = True)

    return MasterDF

def AppendBasinInfoCSVs(DataDirectory, FilenamePrefix):
    """
    This function reads in a series of csvs with the suffix "_AllBasinsInfo"
    and appends them together into one function for plotting

    Args:
        DataDirectory (str): the data DataDirectory

    Returns:
        pandas dataframe with the appended movern csvs

    Author: FJC
    """

    # get the csv filename
    csv_suffix =  '_AllBasinsInfo.csv'

    MasterDF = pd.DataFrame()
    basin_dict = MapBasinsToKeysFromJunctionList(DataDirectory, FilenamePrefix)

    # loop through and get each basin csv
    for outlet_jn, basin_key in basin_dict.iteritems():
        this_fname = "basin"+str(outlet_jn)+csv_suffix
        # append to master DF and change the basin key and the junction
        df = pd.read_csv(DataDirectory+this_fname)
        df = df[df['basin_key'] == 0]
        df['basin_key'] = basin_key
        df['outlet_junction'] = outlet_jn
        MasterDF = MasterDF.append(df, ignore_index = True)

    return MasterDF

def AppendChiDataMapCSVs(DataDirectory, FilenamePrefix):
    """
    This function reads in a series of csvs with the suffix "_chi_data_map"
    and appends them together into one function for plotting

    Args:
        DataDirectory (str): the data DataDirectory

    Returns:
        pandas dataframe with the appended csvs

    Author: FJC
    """

    # get the csv filename
    csv_suffix =  '_chi_data_map.csv'

    MasterDF = pd.DataFrame()
    basin_dict = MapBasinsToKeysFromJunctionList(DataDirectory, FilenamePrefix)

    # loop through and get each basin csv
    for outlet_jn, basin_key in basin_dict.iteritems():
        this_fname = "basin"+str(outlet_jn)+csv_suffix
        # append to master DF and change the basin key and the junction
        df = pd.read_csv(DataDirectory+this_fname)
        df = df[df['basin_key'] == 0]
        df['basin_key'] = basin_key
        MasterDF = MasterDF.append(df, ignore_index = True)

    return MasterDF

def AppendSABinnedCSVs(DataDirectory, fname_prefix):
    """
    This function reads in a series of csvs with the suffix "_SAbinned"
    and appends them together into one function for plotting

    Args:
        DataDirectory (str): the data directory

    Returns:
        pandas dataframe with the appended csvs

    Author: FJC
    """

    # get the csv filename
    csv_suffix =  '_SAbinned.csv'

    MasterDF = pd.DataFrame()
    basin_dict = MapBasinsToKeysFromJunctionList(DataDirectory, fname_prefix)

    # loop through and get each basin csv
    for outlet_jn, basin_key in basin_dict.iteritems():
        this_fname = "basin"+str(outlet_jn)+csv_suffix
        # append to master DF and change the basin key and the junction
        df = pd.read_csv(DataDirectory+this_fname)
        df = df[df['basin_key'] == 0]
        df['basin_key'] = basin_key
        MasterDF = MasterDF.append(df, ignore_index = True)

    # write to a new csv
    MasterDF.to_csv(DataDirectory+fname_prefix+csv_suffix)

    return MasterDF

def AppendSASegmentedCSVs(DataDirectory, fname_prefix):
    """
    This function reads in a series of csvs with the suffix "_SAsegmented"
    and appends them together into one function for plotting

    Args:
        DataDirectory (str): the data directory

    Returns:
        pandas dataframe with the appended csvs

    Author: FJC
    """

    # get the csv filename
    csv_suffix =  '_SAsegmented.csv'

    MasterDF = pd.DataFrame()
    basin_dict = MapBasinsToKeysFromJunctionList(DataDirectory, fname_prefix)

    # loop through and get each basin csv
    for outlet_jn, basin_key in basin_dict.iteritems():
        this_fname = "basin"+str(outlet_jn)+csv_suffix
        # append to master DF and change the basin key and the junction
        df = pd.read_csv(DataDirectory+this_fname)
        df = df[df['basin_key'] == 0]
        df['basin_key'] = basin_key
        MasterDF = MasterDF.append(df, ignore_index = True)

    # write to a new csv
    MasterDF.to_csv(DataDirectory+fname_prefix+csv_suffix)

    return MasterDF

def AppendSAVerticalCSVs(DataDirectory, fname_prefix):
    """
    This function reads in a series of csvs with the suffix "_SAvertical"
    and appends them together into one function for plotting

    Args:
        DataDirectory (str): the data directory

    Returns:
        pandas dataframe with the appended csvs

    Author: FJC
    """

    # get the csv filename
    csv_suffix =  '_SAvertical.csv'

    MasterDF = pd.DataFrame()
    basin_dict = MapBasinsToKeysFromJunctionList(DataDirectory, fname_prefix)

    # loop through and get each basin csv
    for outlet_jn, basin_key in basin_dict.iteritems():
        this_fname = "basin"+str(outlet_jn)+csv_suffix
        # append to master DF and change the basin key and the junction
        df = pd.read_csv(DataDirectory+this_fname)
        df = df[df['basin_key'] == 0]
        df['basin_key'] = basin_key
        MasterDF = MasterDF.append(df, ignore_index = True)

    # write to a new csv
    MasterDF.to_csv(DataDirectory+fname_prefix+csv_suffix)

    return MasterDF

def AppendBasinPointCSVs(DataDirectory, FilenamePrefix):
    """
    This function reads in a series of csvs with the suffix
     "_MCpoint_points_MC_basinstats" and appends them together
     into one function for plotting

    Args:
        DataDirectory (str): the data directory

    Returns:
        pandas dataframe with the appended csvs

    Author: FJC
    """

    # get the csv filename
    csv_suffix =  '_MCpoint_points_MC_basinstats.csv'

    MasterDF = pd.DataFrame()
    basin_dict = MapBasinsToKeysFromJunctionList(DataDirectory, FilenamePrefix)

    # loop through and get each basin csv
    for outlet_jn, basin_key in basin_dict.iteritems():
        this_fname = "basin"+str(outlet_jn)+csv_suffix
        # append to master DF and change the basin key and the junction
        df = pd.read_csv(DataDirectory+this_fname)
        df = df[df['basin_key'] == 0]
        df['basin_key'] = basin_key
        df['outlet_jn'] = outlet_jn
        MasterDF = MasterDF.append(df, ignore_index = True)

    return MasterDF

def AppendChiResidualsCSVs(DataDirectory, FilenamePrefix):
    """
    This function reads in a series of 3 csvs with the residuals data
     and appends them together into one function for plotting

    Args:
        DataDirectory (str): the data directory

    Returns:
        pandas dataframe with the appended csvs

    Author: FJC
    """

    # get the csv filename
    fnames = ["_residual_movernstats_movern_residuals_median.csv","_residual_movernstats_movern_residuals_Q1.csv","_residual_movernstats_movern_residuals_Q3.csv"]

    MasterDFs = []

    MasterDF = pd.DataFrame()
    basin_dict = MapBasinsToKeysFromJunctionList(DataDirectory, FilenamePrefix)

    for f in fnames:
        # loop through and get each basin csv
        for outlet_jn, basin_key in basin_dict.iteritems():
            this_fname = "basin"+str(outlet_jn)+f
            # append to master DF and change the basin key and the junction
            df = pd.read_csv(DataDirectory+this_fname)
            df = df[df['basin_key'] == 0]
            df['basin_key'] = basin_key
            df['outlet_jn'] = outlet_jn
            MasterDF = MasterDF.append(df, ignore_index = True)
        MasterDFs.append(MasterDF)

    return MasterDFs

def MapBasinsToKeys(DataDirectory):
    """
    This function reads in all of the basin files in a directory
    and assigns a basin key to each basin based on the outlet junction.
    This is used when appending these all together so that each basin still
    has a unique ID

    Args:
        DataDirectory (str): the data directory

    Returns:
        dictionary where key is the outlet junction and value is the assigned
        basin ID

    Author: FJC
    """
    import os

    # get the csv filename
    csv_suffix = "_movernstats_basinstats.csv"

    basin_dict = {}
    key = 0

    # loop through the directory and get each basin stats file
    for fname in os.listdir(DataDirectory):
        if fname.startswith("basin"):
            if fname.endswith(csv_suffix):
                # get this basin junction and map to a key
                fname = fname.split("/")[-1]
                fname = fname.split("_")[0]
                fname = fname.split("n")[-1] #stupid way of getting just the basin junction number
                basin_dict[fname] = key
                key+=1

    return basin_dict

def MapBasinsToKeysFromJunctionList(DataDirectory,FilenamePrefix):
    """
    Function to map the basins to keys in an order specified by then
    original junction list.

    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): prefix of the DEM, should be the same as the
        junctions.list file.

    Returns:
        dictionary where key is the outlet junction and value is the assigned
        basin ID

    Author: FJC
    """
    import csv

    basin_dict = {}
    key = 0
    # read in the space delimited junctons.list file
    reader = csv.reader(open(DataDirectory+FilenamePrefix+'_junctions.list'), delimiter=" ")

    for row in reader:
        for jn in row:
            basin_dict[jn] = key
            key+=1

    return basin_dict
