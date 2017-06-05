## analyse_knickpoints.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## This scripts takes in the Mchi files with the knickpoint information
## and creates various plots.  User specifies a threshold knickpoint magnitude
## (difference in MChi between the upstream and downstream segments)
## MChi file is read in using pandas
##
## I'm writing this at midnight so here is a panda:
##             ,,,         ,,,
##          ;"   ^;     ;'   ",
##          ;    s$$$$$$$s     ;
##          ,  ss$$$$$$$$$$s  ,'
##          ;s$$$$$$$$$$$$$$$
##          $$$$$$$$$$$$$$$$$$
##         $$$$P""Y$$$Y""W$$$$$
##         $$$$  p"$$$"q  $$$$$
##         $$$$  .$$$$$.  $$$$
##           $$DcaU$$$$$$$$$$
##            "Y$$$"*"$$$Y"
##               "$b.$$"
##
##
## Knickpoint plots:
## 1. For each basin, it looks at the relationship between the flow distance
## (distance from the basin outlet) and the elevation of the knickpoints
## 2. Knickpoint magnitude vs elevation
##
## IDEAS TO ADD TO PLOTS:
## - Colour points by lithology (import USGS lithologies). Theory predicts that
## - the knickpoints of all the tributaries should be at the same elevation. How does
## lithology affect this?
## - Normalise flow distance by channel length/drainage area (chi coordinate)
## - Different marker for convex/concave knickpoints
## - Plots for along the mountain front: Distance from N - S along the Sierras. How
## does the elevation of the knickpoints vary as you move from N - S?
##
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## Authors: FJC, BG
## 02/06/17
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#set backend to run on server
import matplotlib
matplotlib.use('Agg')

#import modules
from matplotlib import pyplot as plt
import numpy as np
import os
import matplotlib
import pandas
from LSDPlottingTools import LSDMap_PointTools as PointTools

def read_MChi_file(DataDirectory, csv_name, kp_threshold):
    """
    This function reads in the MChi file using pandas
    file structure:
    latitude	longitude chi	elevation	flow distance	drainage area	m_chi	b_chi	source_key	basin_key	segmented_elevation	knickpoints	knickpoint_sign	segment_length	file_from_combine
    FJC 29/03/17
    """
    df = pandas.read_csv(DataDirectory+csv_name, sep=",")
    df = df[df.knickpoints >= kp_threshold]
    return df

def get_data_columns_from_csv(DataDirectory, csv_name, kp_threshold, columns):
    """
    This function returns lists of specified column names from the MChi csv file.
    Must be strings equal to the column headers.
    User can specify a knickpoint threshold (values with knickpoint magnitudes
    below this will be excluded)
    FJC 29/03/17
    """
    column_lists = []
    df = read_MChi_file(DataDirectory, csv_name, kp_threshold)
    for column_name in columns:
        print("I'm returning the "+column_name+" values as a list")
        column_values = list(df[column_name])
        column_lists.append(column_values)
    return column_lists

def make_cumulative_plot(DataDirectory, csv_name, kp_threshold):
    print("Now printing the cumulative plot")
    sorted_data = read_MChi_file(DataDirectory, csv_name, kp_threshold)
    temp_count = 0
    x_cumul, y_cumul = np.unique(sorted_data[:,11],return_counts= True)
    #y_cumul = np.unique(sorted_data[:,11],return_counts= True)
    for i in range(1,x_cumul.size):
        y_cumul[i] = y_cumul[i]+y_cumul[i-1]

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(x_cumul, y_cumul, 'k--', linewidth=1)

    print("Saving the Cumulative plot")
    # tidy up the figure
    ax.grid(True)
    ax.set_title('Cumulative step histograms')
    ax.set_xlabel('kinckpoint value')
    ax.set_ylabel('Cumulative %')
    #ax.set_ylim(0,100)
    ax.set_xlim(0,sorted_data[:,11].max())
    write_name = "kp_cumulative"
    file_ext = "png"
    plt.savefig(DataDirectory+write_name+"."+file_ext,dpi=300)
    plt.clf()

    #### Elevation against Knickpoints ####
def plot_knickpoint_elevations(DataDirectory, csv_name, kp_threshold):
    """
    This function creates a plot of knickpoint elevations against magnitude
    FJC 29/03/17, modified from code by BG.
    """
    # read in the data from the csv to lists
    elevation = get_data_column_from_csv(DataDirectory, csv_name, kp_threshold, "elevation")
    kp_magnitude = get_data_column_from_csv(DataDirectory, csv_name, kp_threshold, "knickpoints")

    # plot the figure
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(kp_magnitude, elevation, 'k+', linewidth=0.5)
    ax.grid(True)
    ax.set_title('Elevation against knickpoint value')
    ax.set_xlabel('Knickpoint magnitude')
    ax.set_ylabel('Elevation (m)')
    #ax.set_ylim(0,100)
    ax.set_xlim(0,1000)
    write_name = "knickpoint_elevation"
    file_ext = "png"
    plt.savefig(DataDirectory+write_name+"."+file_ext,dpi=300)
    plt.clf()

def knickpoint_plots_for_basins(DataDirectory, csv_name, kp_threshold):
    """
    This function creates subplots of knickpoint characteristics for each individual
    basin.
    FJC 29/03/17
    """
    # read in data from the csv to lists
    columns = ["elevation", "flow distance", "drainage area", "knickpoints", "knickpoint_sign", "file_from_combine"]
    column_lists = get_data_columns_from_csv(DataDirectory, csv_name, kp_threshold, columns)
    print(len(column_lists))
    elevation = column_lists[0]
    flow_distance = column_lists[1]
    drainage_area = column_lists[2]
    kp_magnitude = column_lists[3]
    kp_sign = column_lists[4]
    basin_id = column_lists[5]
    #list_of_lists = zip(elevation,flow_distance,basin_id)
    #print column_lists

    # loop through and get a plot for each basin id
    ids = set(basin_id)
    for id in ids:
        print("This basin id is: "+str(id))
        these_lists = [(a,b,c,d,e,f) for (a,b,c,d,e,f) in zip(elevation,flow_distance,drainage_area,kp_magnitude,kp_sign,basin_id) if f == id]
        this_elev, this_distance, this_area, this_magnitude, this_sign, this_id = zip(*these_lists)

        fig,ax = plt.subplots(figsize=(10,12))
        #ax = ax.ravel()
        #for i in range(len(ax)):
        ax.scatter(this_distance,this_elev,facecolors="None", edgecolors="k", s=this_magnitude)
        ax.set_xlabel('Flow distance (m)')
        ax.set_ylabel('Elevation (m)')

        write_name = "knickpoint_plots_basin_"
        file_ext = "png"
        plt.savefig(DataDirectory+write_name+str(id)+"."+file_ext,dpi=100)
        plt.close()

def knickpoint_plotter(DataDirectory, DEM_prefix, kp_threshold=0, FigFormat='pdf'):
    """
    Function to test LSDMap_KnickpointPlotting

    Args:
        DataDirectory (str): the data directory of the chi csv file
        DEM_prefix (str): DEM name without extension
        kp_threshold (int): threshold knickpoint magnitude, values below this will be removed.

    Returns:
        Plot of knickpoint magnitude for each basin

    Author: FJC
    """
    from LSDPlottingTools import LSDMap_PointTools as PointTools
    from LSDPlottingTools import LSDMap_KnickpointPlotting as KP

    # read in the raw csv file
    kp_csv_fname = DataDirectory+DEM_prefix+'_MChi.csv'
    print("I'm reading in the csv file "+kp_csv_fname)

    # get the point data objects
    PointData = PointTools.LSDMap_PointData(kp_csv_fname)

    # get the basin keys
    basin = PointData.QueryData('file_from_combine')
    basin = [int(x) for x in basin]
    Basin = np.asarray(basin)
    basin_keys = np.unique(Basin)
    print('There are %s basins') %(len(basin_keys))

    # loop through each basin and make the figure
    for basin_key in basin_keys:
        FileName = DEM_prefix+"_KP_elevations_%s.%s" %(str(basin_key),FigFormat)
        KP.plot_knickpoint_elevations(PointData, DataDirectory, DEM_prefix, basin_key, kp_threshold, FileName, FigFormat)

if __name__ == "__main__":

    DataDirectory = '/home/s0923330/DEMs_for_analysis/sierras_kn/'
    baseName = "combined"
    #csv_name = baseName + "_MChi.csv"
    kp_threshold = 100 # every knickpoint below this will be erased
    FigFormat = 'png'
    #knickpoint_plots_for_basins(DataDirectory,csv_name, kp_threshold)
    #get_data_column_from_csv(DataDirectory,csv_name,kp_threshold,column_name="latitude")
    knickpoint_plotter(DataDirectory,baseName,kp_threshold=kp_threshold,FigFormat=FigFormat)
