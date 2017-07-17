#=============================================================================
# Script to plot the terrace analysis
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
from LSDPlottingTools import LSDMap_TerracePlotting as Terraces

#=============================================================================
# This is just a welcome screen that is displayed if no arguments are provided.
#=============================================================================
def print_welcome():

    print("\n\n=======================================================================")
    print("Hello! I'm going to plot the terrace analysis for you.")
    print("You will need to tell me which directory to look in.")
    print("Use the -wd flag to define the working directory.")
    print("If you don't do this I will assume the data is in the same directory as this script.")
    print("For help type:")
    print("   python TerracePlotting.py -h\n")
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
    # The location of the data files
    parser.add_argument("-bd", "--base_directory", type=str, help="The base directory with the terrace analysis. If this isn't defined I'll assume it's the same as the current directory.")
    parser.add_argument("-fname", "--fname_prefix", type=str, help="The prefix of your DEM WITHOUT EXTENSION!!! This must be supplied or you will get an error.")

    # What sort of analyses you want
    parser.add_argument("-rasters", "--plot_rasters", type=bool, default=False, help="If this is true, I'll make raster plots of the terrace locations")
    parser.add_argument("-profiles", "--plot_profiles", type=bool, default=False, help="If this is true, I'll make plots of the terrace elevations with distance upstream along the channel.")

    # These control the format of your figures
    parser.add_argument("-fmt", "--FigFormat", type=str, default='png', help="Set the figure format for the plots. Default is png")
    parser.add_argument("-size", "--size_format", type=str, default='ESURF', help="Set the size format for the figure. Can be 'big' (16 inches wide), 'geomorphology' (6.25 inches wide), or 'ESURF' (4.92 inches wide) (defualt esurf).")
    args = parser.parse_args()

    if not args.fname_prefix:
        print("WARNING! You haven't supplied your DEM name. Please specify this with the flag '-fname'")
        sys.exit()

    # get the base directory
    if args.base_directory:
        this_dir = args.base_directory
    else:
        this_dir = os.getcwd()

    # make the plots depending on your choices
    if args.plot_rasters:
        #Terraces.MakeTerraceIDRasterPlot(this_dir, args.fname_prefix, args.size_format, args.FigFormat)
        Terraces.MakeTerraceElevationRasterPlot(this_dir, args.fname_prefix, args.size_format, args.FigFormat)
    if args.plot_profiles:
        Terraces.make_terrace_plots_all_files(this_dir, args.fname_prefix, args.size_format, args.FigFormat)

#=============================================================================
if __name__ == "__main__":
    main(sys.argv[1:])
