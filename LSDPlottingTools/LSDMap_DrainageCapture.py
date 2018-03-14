#=============================================================================
# These functions create figures for visualising the data from the drainage
# capture metrics
#
# It creates separate plots for each basin in the DEM.
#
# Authors:
#     Simon M. Mudd
#     Fiona J. Clubb
#=============================================================================
#=============================================================================
# IMPORT MODULES
#=============================================================================
import numpy as np
import LSDPlottingTools as LSDP
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.ticker as ticker
import pandas as pd
import os

from LSDMapFigure import PlottingHelpers as Helper

#=============================================================================
# PLOTTING FUNCTIONS
#=============================================================================
def PlotBasinPerimeter(DataDirectory, fname_prefix, size_format='ESURF', FigFormat='png'):
    """
    Make a plot of the basin perimeter ordered by the outlet

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): filename of the DEM without extension
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.

    Author: FJC
    """
    # check if a directory exists for the perimeter plots. If not then make it.
    this_dir = DataDirectory+'basin_perimeters/'
    if not os.path.isdir(this_dir):
        os.makedirs(this_dir)

    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        #l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        #l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
        #l_pad = -35

    PerimeterDF = Helper.ReadPerimeterCSV(DataDirectory, fname_prefix)

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.05,right=0.95,top=0.95)
    ax = fig.add_subplot(gs[5:100,10:95])

    # plot the data
    ax.plot(PerimeterDF['node_key'],PerimeterDF['elevation'])

    # set the axis labels
    ax.set_xlabel('Perimeter node ordered from outlet')
    ax.set_ylabel('Node elevation')

    newFilename = this_dir+fname_prefix+"_basin_perimeter."+FigFormat
    plt.savefig(newFilename,format=FigFormat,dpi=300)
    ax.cla()
    plt.close(fig)
