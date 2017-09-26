"""
This function reads hillshades and elevation files in a directory


INSTRUCTIONS FOR ffmpeg
This will produce a series of images that you can convert to a movie using ffmpeg

Here is a typical command line
ffmpeg -framerate 5 -pattern_type glob -i '*.png' -vcodec libx264 -s 1230x566 -pix_fmt yuv420p movie7.mp4

Some gotchas:
The -framerate flag must come after ffmpeg (I think)
You can be more specific about the pattern but this seems to work
For the -s flag you must have pixel sizes that are divisible by 2.

Written by Simon Mudd
June 2017
git
GPL3
"""


import LSDChiMappingExamples as CME
import LSDMapWrappers as MW
import os
from sys import platform
import sys


def get_filenames(root):
    # Create and empty list for the filenames

    print("The root is: "+root)

    these_filenames = []
    n = 0
    while 1:
        filename = "%s%d.bil" % (root, n+1)
        print("The filename is: "+filename)
        if os.path.isfile(filename):
             n += 1
             these_filenames.append(filename)
        else:
            return these_filenames

def run_plots(DataDirectory,Base_file, cbar_min_max = [0,400]):

    root = DataDirectory+Base_file
    filenames = get_filenames(root)
    counter = 0
    for fname in filenames:

        counter = counter+1
        # remove the extension form the files
        this_base_file = fname[:-4]
        if "win" in platform:
            split_file = this_base_file.split("\\")
        elif "linux" in platform:
            split_file = this_base_file.split("/")
        print("Split file is: ")
        print split_file
        this_base_file = split_file[-1]
        print("This base file is: "+ this_base_file)

        print("I am getting figures for the animation.")
        MW.SimpleHillshadeForAnimation(DataDirectory,this_base_file,cmap = "terrain", dpi = 250, imgnumber=counter, full_basefile = root,custom_cbar_min_max = cbar_min_max)

def animate_plots(base_directory, fname_prefix):
    """
    This function creates a movie of the plots using ffmpeg

    Args:
        base_directory (str): the directory with the plots.
        fname_prefix (str): the filename for the model run

    Returns:
        none but creates mp4 from pngs in base directory

    Author: FJC
    """
    import subprocess

    # animate the pngs using ffmpeg
    system_call = "ffmpeg -framerate 5 -pattern_type glob -i '"+base_directory+"*.png' -vcodec libx264 -s 1230x566 -pix_fmt yuv420p "+base_directory+fname_prefix+"_movie.mp4"
    print system_call
    subprocess.call(system_call, shell=True)

#=============================================================================
# This is just a welcome screen that is displayed if no arguments are provided.
#=============================================================================
def print_welcome():

    print("\n\n=======================================================================")
    print("Hello! I'm going to plot a series of hillshades for you.")
    print("You will need to tell me the directory and the base file name.")
    print("Use the -dir flag to define the working directory.")
    print("If you don't do this I will assume the data is in the same directory as this script.")
    print("Use the -fname flag to define the base file name.")
    print("For help type:")
    print("   python PlotMOverNAnalysis.py -h\n")
    print("=======================================================================\n\n ")

#=============================================================================
# This is the main function that runs the whole thing
#=============================================================================
def main(argv):

    """
    This is just a few lines for keeping track of how long the program is taking.
    You can ignore it.
    """
    import time
    tic = time.clock()

    # If there are no arguments, send to the welcome screen
    if not len(sys.argv) > 1:
        full_paramfile = print_welcome()
        sys.exit()

    # Get the arguments
    import argparse
    parser = argparse.ArgumentParser()
    # The location of the data files
    parser.add_argument("-dir", "--base_directory", type=str, help="The base directory with the hillshades. If this isn't defined I'll assume it's the same as the current directory.")
    parser.add_argument("-fname", "--fname_prefix", type=str, help="The base file name of the hillshades.")
    parser.add_argument("-animate", "--animate", type=bool, default=False, help="If this is true I'll create a movie of the model run.")
    parser.add_argument("-zmax", "--maximum_elevation_for_plotting", type=float, default = 400, help="This is the maximum elevation in the colourbar of the landscape plot.")
    args = parser.parse_args()

    cbar_min_max = [0,args.maximum_elevation_for_plotting]
    
    run_plots(args.base_directory,args.fname_prefix,cbar_min_max)

    if (args.animate):
        animate_plots(args.base_directory, args.fname_prefix)

    toc = time.clock()
    print("This took: "+str(toc - tic)+" units of time")


#=============================================================================
if __name__ == "__main__":
    main(sys.argv[1:])
