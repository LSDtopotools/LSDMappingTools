#=============================================================================
# Script to plot the knickpoint data produced with the LSDTT knickpoint picking algorithm 
#
# Authors:
#     Boris Gailleton, Fiona J. Clubb
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
from LSDPlottingTools import LSDMap_BasicMaps as BP
from LSDMapFigure import PlottingHelpers as Helper
from LSDPlottingTools import LSDMap_KnickpointPlotting as KP

#=============================================================================
# This is just a welcome screen that is displayed if no arguments are provided.
#=============================================================================
def print_welcome():

    print("\n\n=======================================================================")
    print("Hello! I'm going to plot some knickpoints stuffs for you.")
    print("You will need to tell me which directory to look in.")
    print("Use the -wd flag to define the working directory.")
    print("If you don't do this I will assume the data is in the same directory as this script.")
    print("I also need to know the common prefix of all your files generated whith LSDTopoTool")
    print("For help type:")
    print("   python PlotknickpointAnalysis.py -h\n")
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
    parser.add_argument("-dir", "--base_directory", type=str, help="The base directory with the m/n analysis. If this isn't defined I'll assume it's the same as the current directory.")
    parser.add_argument("-fname", "--fname_prefix", type=str, help="The prefix of your DEM WITHOUT EXTENSION!!! This must be supplied or you will get an error.")

    # What sort of analyses you want
    parser.add_argument("-mb", "--map_basic", type=bool, default = False, help="Turn to True to plot a basic knickpoint map on the top of the hillshade of the field")
    parser.add_argument("-bh", "--basic_hist", type=bool, default = False, help="Turn to True to plot a basic histogram of the knickpoint spreading")

    # Data sorting option
    parser.add_argument("-mancut", "--manual_cutoff", type=float, default = 0, help="Set a manual cutoff value for plotting the basic maps without automatic stat")

    # Basin
    # Basin selection stuffs
    parser.add_argument("-basin_keys", "--basin_keys",type=str,default = "", help = "This is a comma delimited string that gets the list of basins you want for the plotting. Default = no basins")
    
    

    # These control the format of your figures
    parser.add_argument("-fmt", "--FigFormat", type=str, default='png', help="Set the figure format for the plots. Default is png")
    parser.add_argument("-size", "--size_format", type=str, default='ESURF', help="Set the size format for the figure. Can be 'big' (16 inches wide), 'geomorphology' (6.25 inches wide), or 'ESURF' (4.92 inches wide) (defualt esurf).")
    args = parser.parse_args()

    print("You told me that the basin keys are: ")
    print(args.basin_keys)

    if len(args.basin_keys) == 0:
        print("No basins found, I will plot all of them")
        these_basin_keys = []
    else:
        these_basin_keys = [int(item) for item in args.basin_keys.split(',')]
        print("The basins I will plot are:")
        print(these_basin_keys)

    if not args.fname_prefix:
        print("WARNING! You haven't supplied your DEM name. Please specify this with the flag '-fname'")
        sys.exit()


##################### Plotting facilities

    if args.map_basic:
        KP.map_knickpoint_standard(args.base_directory, args.fname_prefix, basin_list = these_basin_keys, size_format=args.size_format, FigFormat=args.FigFormat, mancut = args.manual_cutoff)

    if args.basic_hist:
        KP.basic_hist(args.base_directory, args.fname_prefix,basin_list = these_basin_keys, size_format=args.size_format, FigFormat=args.FigFormat)

#=============================================================================
if __name__ == "__main__":
    main(sys.argv[1:])

