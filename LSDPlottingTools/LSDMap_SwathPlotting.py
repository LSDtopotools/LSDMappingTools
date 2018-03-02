## LSDMap_SwathPlotting.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## These functions are tools to deal with plotting swaths
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## SMM
## 20/02/2018
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from __future__ import absolute_import, division, print_function

import numpy as np
import os
from . import cubehelix
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
from matplotlib import colors
import LSDPlottingTools.LSDMap_GDALIO as LSDMap_IO
import LSDPlottingTools.LSDMap_BasicPlotting as LSDMap_BP
import LSDPlottingTools.LSDMap_PointTools as LSDMap_PD
from LSDMapFigure.PlottingRaster import MapFigure
from LSDMapFigure.PlottingRaster import BaseRaster
from LSDMapFigure import PlottingHelpers as Helper

def PlotSwath(swath_csv_name, FigFileName = 'Image.png', size_format = "geomorphology", fig_format = "png", dpi = 500, aspect_ratio = 2):
    """
    This plots a swath profile
    
    Args: 
        swath_csv_name (str): the name of the csv file (with path!)
        
    Author: SMM
    
    Date 20/02/2018
    """
    
    print("STARTING swath plot.")

    # Set up fonts for plots
    label_size = 12
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure,
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        fig_size_inches = 6.25
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        fig_size_inches = 16
        l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))
        fig_size_inches = 4.92126
        l_pad = -35
    
    # Note all the below parameters are overwritten by the figure sizer routine
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])

    print("Getting data from the file: "+swath_csv_name)
    thisPointData = LSDMap_PD.LSDMap_PointData(swath_csv_name)
    
    distance = thisPointData.QueryData('Distance').values
    mean_val = thisPointData.QueryData('Mean').values
    min_val = thisPointData.QueryData('Min').values
    max_val = thisPointData.QueryData('Max').values
    
    # Get the minimum and maximum distances
    X_axis_min = 0
    X_axis_max = distance[-1]
    n_target_tics = 5
    xlocs,new_x_labels = LSDMap_BP.TickConverter(X_axis_min,X_axis_max,n_target_tics)
    
    ax.fill_between(distance, min_val, max_val, facecolor='orange', alpha = 0.5, interpolate=True)
    ax.plot(distance, mean_val,"b", linewidth = 1)
    ax.plot(distance, min_val,"k",distance,max_val,"k",linewidth = 1)
    
    ax.spines['top'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)

    ax.set_ylabel("Elevation (m)")
    ax.set_xlabel("Distance along swath (km)")
    
    ax.set_xticks(xlocs)
    ax.set_xticklabels(new_x_labels,rotation=60)

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
    ax.tick_params(axis='both', width=1, pad = 2)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(2)

        
    # Lets try to size the figure
    cbar_L = "None"
    [fig_size_inches,map_axes,cbar_axes] = Helper.MapFigureSizer(fig_size_inches,aspect_ratio, cbar_loc = cbar_L, title = "None")
    
    fig.set_size_inches(fig_size_inches[0], fig_size_inches[1])
    ax.set_position(map_axes)  
    
    FigFormat = fig_format    
    print("The figure format is: " + FigFormat)
    if FigFormat == 'show':
        plt.show()
    elif FigFormat == 'return':
        return fig
    else:
        plt.savefig(FigFileName,format=FigFormat,dpi=dpi)
        fig.clf()
    
    