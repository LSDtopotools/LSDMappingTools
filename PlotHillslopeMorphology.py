#=============================================================================
# Script to plot the m/n analysis using the MLE collinearity
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

#from __future__ import print_function
import sys
import os
from LSDPlottingTools import LSDMap_HillslopeMorphology as HS
from LSDMapFigure import PlottingHelpers as Helper

#=============================================================================
# This is just a welcome screen that is displayed if no arguments are provided.
#=============================================================================
def print_welcome():

    print("\n\n=======================================================================")
    print("Hello! I'm going to plot the hillslope-channel results for you.")
    print("You will need to tell me which directory to look in.")
    print("Use the -dir flag to define the working directory.")
    print("If you don't do this I will assume the data is in the same directory as this script.")
    print("For help type:")
    print("   python PlotMOverNAnalysis.py -h\n")
    print("=======================================================================\n\n ")

#=============================================================================
# This is the main function that runs the whole thing
#=============================================================================
def main(argv):

    # print("On some windows systems you need to set an environment variable GDAL_DATA")
    # print("If the code crashes here it means the environment variable is not set")
    # print("Let me check gdal enviroment for you. Currently is is:")
    # print(os.environ['GDAL_DATA'])
    #os.environ['GDAL_DATA'] = os.popen('gdal-config --datadir').read().rstrip()
    #print("Now I am going to get the updated version:")
    #print(os.environ['GDAL_DATA'])

    # If there are no arguments, send to the welcome screen
    if not len(sys.argv) > 1:
        full_paramfile = print_welcome()
        sys.exit()

    # Get the arguments
    import argparse
    parser = argparse.ArgumentParser()

    # The location of the data files
    parser.add_argument("-dir", "--base_directory", type=str, help="The base directory with the m/n analysis. If this isn't defined I'll assume it's the same as the current directory.")
    parser.add_argument("-fname", "--fname_prefix", type=str, help="The prefix of your DEM WITHOUT EXTENSION!!! This must be supplied or you will get an error (unless you're running the parallel plotting).")

    # The basins you want to plot
    parser.add_argument("-basin_keys", "--basin_keys",type=str,default = "", help = "This is a comma delimited string that gets the list of basins you want for the plotting. Default = no basins")

    # The analysis you want to do
    parser.add_argument("-Mchi", "--plot_mchi", type=bool, default=False, help="If this is true, I'll make some plots of the hillslope-channel data against Mchi")
    parser.add_argument("-CHT", "--plot_CHT", type=bool, default=False, help="If this is true, I'll make some plots of hilltop curvature against data from the channel segments.")
    parser.add_argument("-segments", "--plot_segments", type=bool, default=False, help="If this is true, I'll make some long profile plots of the channel segments.")
    parser.add_argument("-in_basin", "--plot_data_within_basin", type=bool, default=False, help="If this is true, I'll make plots of the hillslope data vs distance upstream for each basin")
    parser.add_argument("-means", "--plot_mean_basin_data", type=bool, default=False, help="If this is true I'll make plots of the mean hillslope data vs basin ID")

    args = parser.parse_args()

    if not args.fname_prefix:
        print("WARNING! You haven't supplied your DEM name. Please specify this with the flag '-fname'")
        sys.exit()

    # get the base directory
    if args.base_directory:
        this_dir = args.base_directory
        # check if you remembered a / at the end of your path_name
        if not this_dir.endswith("/"):
            print("You forgot the '/' at the end of the directory, appending...")
            this_dir = this_dir+"/"
    else:
        this_dir = os.getcwd()

    # check the basins
    print("You told me that the basin keys are: ")
    print(args.basin_keys)

    if len(args.basin_keys) == 0:
        print("No basins found, I will plot all the basins")
        df = HS.ReadChannelData(this_dir, args.fname_prefix)
        these_basin_keys = list(df['basin_key'].unique())
        print (these_basin_keys)
    else:
        these_basin_keys = [int(item) for item in args.basin_keys.split(',')]
        print("The basins I will plot are:")
        print(these_basin_keys)


    # separate the data into basin csvs
    HS.SaveHillslopeDataByBasin(this_dir, args.fname_prefix)
    HS.SaveChannelDataByBasin(this_dir, args.fname_prefix)

    # make the directory for saving the plots
    PlotDirectory = this_dir+'hillslope_plots/'
    if not os.path.isdir(PlotDirectory):
        os.makedirs(PlotDirectory)

    if args.plot_mchi:
        for basin_key in these_basin_keys:
            HS.PlotChiElevationMChi(this_dir, args.fname_prefix, PlotDirectory, basin_key)
            HS.PlotLongProfileMChi(this_dir, args.fname_prefix, PlotDirectory, basin_key)
    if args.plot_CHT:
        for basin_key in these_basin_keys:
            HS.PlotCHTAgainstChannelData(this_dir, args.fname_prefix, PlotDirectory, basin_key)
    if args.plot_segments:
        for basin_key in these_basin_keys:
            HS.PlotLongProfileSegments(this_dir, args.fname_prefix, PlotDirectory, basin_key)
            HS.PlotChiElevationSegments(this_dir, args.fname_prefix, PlotDirectory, basin_key)
    if args.plot_data_within_basin:
        for basin_key in these_basin_keys:
            print "This basin key is: ", basin_key
            # HS.PlotEStarRStar(this_dir, args.fname_prefix, PlotDirectory, basin_key)
            # HS.PlotRStarDistance(this_dir, args.fname_prefix, PlotDirectory,basin_key)
            # HS.PlotLHDistance(this_dir, args.fname_prefix, PlotDirectory, basin_key)
            HS.PlotHillslopeDataVsDistance(this_dir, args.fname_prefix, PlotDirectory, basin_key)
    if args.plot_mean_basin_data:
        HS.PlotHillslopeDataWithBasins(this_dir, args.fname_prefix, PlotDirectory)


#=============================================================================
if __name__ == "__main__":
    main(sys.argv[1:])
