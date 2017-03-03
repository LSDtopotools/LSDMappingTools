"""
This script plots (or will hopefully plot) Max/Mean/Min Elevation vs distance
from the Swath(s) driver of LSDTopoTools.
It takes the csv files generated from the driver and and the original raster.
Written By B.Gailleton and F.J.Clubb
"""

#import modules
import numpy as np, matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.cm as cmx
#import seaborn as sns
import matplotlib.colors as colors
#from rand_cmap import rand_cmap
from matplotlib import gridspec
#from adjustText import adjust_text
from scipy.interpolate import spline, UnivariateSpline

# Go to the end of this script to enter your parameters
#
#
#
################### Development part ############################
######### Do not change these parameters unless you want to #####
#################################################################
def cm2inch(value):
    return value/2.54

def read_all_data(DataDirectory,DEM_prefix, info_file, data_file):
    """
    Read in the data with all the points
    """

    # read in the file
    FileName = DataDirectory + data_file+'.csv'
    data = np.loadtxt(FileName, skiprows = 1, delimiter = ",") # Distance,Mean,Min,Max
    no_lines = data.size
    print "Number of lines: ", no_lines
    # Set up the arrays
    distance = data[:,0]
    MeanValues = data[:,1]
    MinValues = data[:,2]
    MaxValues = data[:,3]

    FileName = DataDirectory + info_file+'.csv'
    data2 = np.loadtxt(FileName, skiprows = 1, delimiter = ",") # latitudeA,longitudeA,latitudeB,longitudeB,general_mean_elevation,general_min_elevation,general_max_elevation,Xmin,Ymin

    return distance, MeanValues, MinValues, MaxValues, no_lines, data2

def make_elevation_plots(DataDirectory, DEM_prefix, DEM_swathed,info_file, data_file,  field_site, smoothing, OutputFigureName):
    """
    Make nice plots of the elevations along the swath
    """

    # Set up fonts
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = 14

    distance, MeanValues, MinValues, MaxValues, Nlines, info = read_all_data(DataDirectory, DEM_prefix,info_file, data_file)

    # add the XY Plot
    fig = plt.figure(1, figsize=(cm2inch(25),cm2inch(15)))
    # set ratios for subplots
    ax = fig.add_subplot(111)


    distance = np.array(distance)
    if(smoothing):

        MeanValues = np.array(MeanValues)
        MinValues = np.array(MinValues)
        MaxValues = np.array(MaxValues)
        smooth_k = 0.7


        distance_interpolate = np.linspace(distance.min(),distance.max(),num = Nlines * smooth_k)
        MeanValues_smooth = spline(distance,MeanValues,distance_interpolate,order = 1)
        MinValues_smooth = spline(distance,MinValues,distance_interpolate,order = 1)
        MaxValues_smooth = spline(distance,MaxValues,distance_interpolate,order = 1)

        ax.plot(distance_interpolate, MeanValues_smooth, 'k-', lw=0.75, label = "Mean")
        ax.plot(distance_interpolate, MinValues_smooth, '-',c = 'gray', lw=0, label = "Max/min")
        ax.plot(distance_interpolate, MaxValues_smooth, '-',c = 'gray', lw=0,)
        ax.fill_between(distance_interpolate, MaxValues_smooth, MinValues_smooth, color = '#D5D5D5')
        plt.xlabel('Distance along swath (m)')
        plt.ylabel('Elevation (m)')
    else:
        ax.plot(distance, MeanValues, 'k-', lw=0.75, label = "Mean")
        ax.plot(distance, MinValues, '-',c = 'gray', lw=0, label = "Max/min")
        ax.plot(distance, MaxValues, '-',c = 'gray', lw=0,)
        ax.fill_between(distance, MaxValues, MinValues, color = '#D5D5D5')
        plt.xlabel('Distance along swath (m)')
        plt.ylabel('Elevation (m)')


    temp = np.array(MinValues)
    temp2 = np.array(MaxValues)

    print('I am plotting the data, the minimum altitude is %s and the maximum is %s.' %(temp.min(),temp2.max()))

    plt.ylim(temp.min() - temp.min()/20 ,temp2.max() + temp2.max()/20)
    plt.xlim(0,distance.max())
    plt.title(field_site)
    plt.legend(loc='upper right')

    # Save figure
    OutputFigureFormat = 'png'
    plt.savefig(DataDirectory + OutputFigureName + '.' + OutputFigureFormat, format=OutputFigureFormat, dpi=1000)
    #plt.tight_layout()
    plt.close()

#############################################################################################################
######################################## Parameters to change ###############################################
#############################################################################################################

if __name__ == "__main__":

    DEM_prefix = 'SRTM90_final_PP' # The original DEM
    DEM_swathed = DEM_prefix + "_swath_raster_0"  # the swath raster produce by the driver, change it if you manually change the outpu file
    info_file  = DEM_prefix + "_swath_genInfo_0" # the information file, change it if you manually change the outpu file
    data_file = DEM_prefix + "_swath_elevations_0" # the data file, change it if you manually change the outpu file
    field_site = "Swath profile" # Title of the figure
    smoothing = False # if true, slightly smooth the curve
    #DataDirectory = 'Z:\\DEMs_for_analysis\\eel_river\\'
    DataDirectory = '/home/s1675537/PhD/DataStoreBoris/GIS/Data/Carpathian/Best_DEM/'
    OutputFigureName = DEM_prefix+'_elevation_plots'
    make_elevation_plots(DataDirectory, DEM_prefix, DEM_swathed,info_file, data_file,  field_site, smoothing, OutputFigureName)
