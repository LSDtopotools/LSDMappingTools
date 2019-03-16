"""
Create 3d elevation plots of model runs and then animate them
"""


from LSDPlottingTools import LSDMap_GDALIO as IO
import matplotlib.pyplot as plt
from mayavi import mlab
from mayavi.modules.grid_plane import GridPlane
from matplotlib import cm
import numpy as np
import os
from sys import platform
import sys
from glob import glob


def get_filenames(root):
    # Create and empty list for the filenames
    these_filenames = []
    for filename in glob(root+'*.bil'):
        if not 'hs' in filename:
            if not 'Raster' in filename:
                these_filenames.append(filename)

    these_filenames.sort()
    print(these_filenames)
    return these_filenames

def run_plots(DataDirectory,Base_file):

    root = DataDirectory+Base_file
    filenames = get_filenames(root)
    counter = 0

    # create the plot for the initial raster
    initial_file = filenames[0]

    # read in the raster
    raster = IO.ReadRasterArrayBlocks(initial_file)

    f = mlab.figure(size=(1000,1000), bgcolor=(0.5,0.5,0.5))
    s = mlab.surf(raster, warp_scale=0.4, colormap='gist_earth', vmax=100)
    #mlab.outline(color=(0,0,0))

    #mlab.axes(s, color=(1,1,1), z_axis_visibility=True, y_axis_visibility=False, xlabel='', ylabel='', zlabel='', ranges=[0,500,0,1000,0,0])

    #@mlab.animate(delay=10)
    #def anim():
    # now loop through each file and update the z values
    for fname in filenames:
        this_rast = IO.ReadRasterArrayBlocks(fname)
        s.mlab_source.scalars = this_rast
        #f.scene.render()
        #
        mlab.savefig(fname[:-4]+'_3d.png')
        #mlab.clf()

    # for (x, y, z) in zip(xs, ys, zs):
    #     print('Updating scene...')
    #     plt.mlab_source.set(x=x, y=y, z=z)
    #     yield



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
    system_call = "ffmpeg -framerate 5 -pattern_type glob -i '"+base_directory+"*.png' -vcodec libx264 -s 1000x1000 -pix_fmt yuv420p "+base_directory+fname_prefix+"_movie.mp4"
    print(system_call)
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

    run_plots(args.base_directory,args.fname_prefix)

    if (args.animate):
        animate_plots(args.base_directory, args.fname_prefix)

    toc = time.clock()
    print("This took: "+str(toc - tic)+" units of time")


#=============================================================================
if __name__ == "__main__":
    main(sys.argv[1:])
