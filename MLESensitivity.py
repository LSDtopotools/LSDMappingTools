#=============================================================================
# Runs visulisation on each sensitivity test for the MLE collinearity analysis.
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
from LSDPlottingTools import LSDMap_MOverNPlotting as MN

#=============================================================================
# This is just a welcome screen that is displayed if no arguments are provided.
#=============================================================================
def print_welcome():

    print("\n\n=======================================================================")
    print("Hello! I'm going to plot the MLE sensitivity results for you.")
    print("You will need to tell me which directory to look in.")
    print("Use the -wd flag to define the working directory.")
    print("If you don't do this I will assume the data is in the same directory as this script.")
    print("For help type:")
    print("   python MLESensitivity.py -h\n")
    print("=======================================================================\n\n ")

#=============================================================================
# This is the main function that runs the whole thing
#=============================================================================
def main(argv):

    # If there are no arguments, send to the welcome screen
    if not len(sys.argv) > 1:
        full_paramfile = print_welcome()
        sys.exit()

    # Get the arguments
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-bd", "--base_directory", type=str, help="The base directory with the MLE analyses. If this isn't defined I'll assume it's the same as the current directory.")
    parser.add_argument("-fname", "--fname_prefix", type=str, help="The prefix of your DEM WITHOUT EXTENSION!!! This must be supplied or you will get an error.")
    parser.add_argument("-PR", "--plot_rasters", type=bool, default=False, help="If this is true, I'll make raster plots of the m/n value and basin keys")
    parser.add_argument("-PC", "--plot_chi_profiles", type=bool, default=False, help="If this is true, I'll make chi-elevation plots for each basin coloured by the MLE")
    parser.add_argument("-PO", "--plot_outliers", type=bool, default=False, help="If this is true, I'll make chi-elevation plots with the outliers removed")
    parser.add_argument("-MLE", "--plot_MLE_movern", type=bool, default=False, help="If this is true, I'll make a plot of the MLE values for each m/n showing how the MLE values change as you remove the tributaries")
    parser.add_argument("-start_movern", "--start_movern", type=float, default=0.2, help="Define the starting m/n value for testing, default = 0.2")
    parser.add_argument("-d_movern", "--d_movern", type=float, default=0.1, help="Define the change in m/n value for testing, default = 0.1")
    parser.add_argument("-n_movern", "--n_movern", type=int, default=7, help="Define the number of m/n values for testing, default = 7")
    parser.add_argument("-fmt", "--FigFormat", type=str, default='png', help="Set the figure format for the plots. Default is png")
    parser.add_argument("-size", "--size_format", type=str, default='ESURF', help="Set the size format for the figure. Can be 'big' (16 inches wide), 'geomorphology' (6.25 inches wide), or 'ESURF' (4.92 inches wide) (defualt esurf).")
    args = parser.parse_args()

    if not args.fname_prefix:
        print("WARNING! You haven't supplied your DEM name. Please specify this with the flag '-fname'")
        sys.exit()

    # get the base directory
    if args.base_directory:
        Directory = args.base_directory
    else:
        Directory = os.getcwd()

    # loop through each sub-directory with the sensitivity results
    MLE_str = "Chi_analysis_sigma_"
    for subdir, dirs, files in os.walk(Directory):
        for dir in dirs:
            if MLE_str in dir:
                this_dir = Directory+"/"+dir+'/'
                # make the plots depending on your choices
                if args.plot_rasters:
                    MN.MakeRasterPlotsBasins(this_dir, args.fname_prefix, args.size_format, args.FigFormat)
                    MN.MakeRasterPlotsMOverN(this_dir, args.fname_prefix, args.n_movern, args.d_movern, args.size_format, args.FigFormat)
                if args.plot_chi_profiles:
                    MN.MakeChiPlotsMLE(this_dir, args.fname_prefix, basin_list=[], start_movern=args.start_movern, d_movern=args.d_movern, n_movern=args.n_movern, size_format=args.size_format, FigFormat = args.FigFormat)
                if args.plot_outliers:
                    MN.PlotProfilesRemovingOutliers(this_dir, args.fname_prefix, basin_list=[], start_movern=args.start_movern, d_movern=args.d_movern, n_movern=args.n_movern)
                if args.plot_MLE_movern:
                    MN.PlotMLEWithMOverN(this_dir, args.fname_prefix,basin_list=[], start_movern=args.start_movern, d_movern=args.d_movern, n_movern=args.n_movern, size_format=args.size_format, FigFormat = args.FigFormat)

    # collate all the results to get the final figure
    MN.PlotSensitivityResults(Directory, args.fname_prefix,FigFormat=args.FigFormat,size_format=args.size_format)

#=============================================================================
if __name__ == "__main__":
    main(sys.argv[1:])
