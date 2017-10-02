#=============================================================================
# Script to plot the litho information acquired from 
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

#=============================================================================
# This is just a welcome screen that is displayed if no arguments are provided.
#=============================================================================
def print_welcome():

    print("\n\n=======================================================================")
    print("Hello! I'm going to plot some litho stuffs for you.")
    print("You will need to tell me which directory to look in.")
    print("Use the -wd flag to define the working directory.")
    print("If you don't do this I will assume the data is in the same directory as this script.")
    print("I also need to know some basic info like the prefix of your base dem and the name of your lithologic DEM")
    print("For help type:")
    print("   python PlotLithoInfo.py -h\n")
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
    parser.add_argument("-fname", "--fname_prefix", type=str, help="The prefix of your DEM WITHOUT EXTENSION!!! This must be supplied or you will get an error.")

    # What sort of analyses you want
    parser.add_argument("-t", "--topo", type=bool, default=True, help="Turn to True to plot a basic topographic map with the Terrain colormap")

    # These control the format of your figures
    parser.add_argument("-fmt", "--FigFormat", type=str, default='png', help="Set the figure format for the plots. Default is png")
    parser.add_argument("-size", "--size_format", type=str, default='ESURF', help="Set the size format for the figure. Can be 'big' (16 inches wide), 'geomorphology' (6.25 inches wide), or 'ESURF' (4.92 inches wide) (defualt esurf).")
    args = parser.parse_args()

    if not args.fname_prefix:
        print("WARNING! You haven't supplied your DEM name. Please specify this with the flag '-fname'")
        sys.exit()

    if args.topo:
        BP.PlotTopoRaster(args.base_directory,args.fname_prefix, FigFormat = args.FigFormat, size_format = args.size_format)


#=============================================================================
if __name__ == "__main__":
    main(sys.argv[1:])

