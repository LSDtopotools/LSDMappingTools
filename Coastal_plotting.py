# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 10:25:29 2016

@author: smudd
"""

import numpy as np
import LSDPlottingTools as LSDP
import LSDOSystemTools as LSDOst
from matplotlib import rcParams

def ElevationSwaths(path, filename, axis):
    
    # get the path to the raster file
    NewPath = LSDOst.AppendSepToDirectoryPath(path)    
    FileName = NewPath+filename
        
    # get the data vectors
    means,medians,std_deviations,twentyfifth_percentile,seventyfifth_percentile = LSDP.SimpleSwath(path, filename, axis)
    
    print "Means shape is: "
    print means.shape    
    
    x_vec,y_vec = LSDP.GetLocationVectors(FileName)
    
    
    print "X shape is: "
    print x_vec.shape
    
    print "Y shape is: "
    print y_vec.shape
    
    import matplotlib.pyplot as plt
    import matplotlib.lines as mpllines
    from mpl_toolkits.axes_grid1 import AxesGrid

    label_size = 20
    #title_size = 30
    axis_size = 28

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure, sized for a ppt slide
    fig = plt.figure(1, facecolor='white',figsize=(10,7.5)) 

    gs = plt.GridSpec(100,75,bottom=0.1,left=0.1,right=0.9,top=1.0)
    ax = fig.add_subplot(gs[10:100,10:75])
    
    if axis == 0:
        dir_vec = x_vec
    else:
        dir_vec = y_vec
        
    # get the distance from shore
    dist_from_shore = np.subtract(dir_vec[-1],dir_vec)        
        
    min_sd = np.subtract(means,std_deviations)
    plus_sd = np.add(means,std_deviations) 
        
    ax.plot(dist_from_shore,means, linewidth = 2.5, color = "black")
    #ax.fill_between(dist_from_shore, twentyfifth_percentile, seventyfifth_percentile, facecolor='green', alpha = 0.7, interpolate=True)
    ax.fill_between(dist_from_shore, min_sd, plus_sd, facecolor='blue', alpha = 0.25, interpolate=True)  
    
    ax.set_xlim(dist_from_shore[0],dist_from_shore[-1])

    ax.annotate('Standard deviation envelope', xy=(dist_from_shore[10],plus_sd[10]), xycoords='data',
                xytext=(0.1, 0.8), textcoords='axes fraction',
                size=label_size,
                # bbox=dict(boxstyle="round", fc="0.8"),
                arrowprops=dict(arrowstyle="simple",
                                fc="0.6", ec="none",
                                connectionstyle="arc3,rad=0.3"),
                )


    ax.spines['top'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2) 
    #ax.tick_params(axis='both', width=1) 

    plt.xlabel('Distance from shore (m)', fontsize = axis_size)
    plt.ylabel('Bed elevation relative to MSL (m)', fontsize = axis_size)

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap        
    ax.tick_params(axis='both', width=2, pad = 10)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(10)   
    
    plt.show()

def BedThickSwaths(path, filename, axis):
    
    # get the path to the raster file
    NewPath = LSDOst.AppendSepToDirectoryPath(path)    
    FileName = NewPath+filename
        
    # get the data vectors
    means,medians,std_deviations,twentyfifth_percentile,seventyfifth_percentile = LSDP.SimpleSwath(path, filename, axis)
    
    print "Means shape is: "
    print means.shape    
    
    x_vec,y_vec = LSDP.GetLocationVectors(FileName)
    
    
    print "X shape is: "
    print x_vec.shape
    
    print "Y shape is: "
    print y_vec.shape
    
    import matplotlib.pyplot as plt
    import matplotlib.lines as mpllines
    from mpl_toolkits.axes_grid1 import AxesGrid

    label_size = 20
    #title_size = 30
    axis_size = 28

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure, sized for a ppt slide
    fig = plt.figure(1, facecolor='white',figsize=(10,7.5)) 

    gs = plt.GridSpec(100,75,bottom=0.1,left=0.1,right=0.9,top=1.0)
    ax = fig.add_subplot(gs[10:100,10:75])
    
    if axis == 0:
        dir_vec = x_vec
    else:
        dir_vec = y_vec
        
    # get the distance from shore
    dist_from_shore = np.subtract(dir_vec[-1],dir_vec)        
        
    min_sd = np.subtract(means,std_deviations)
    plus_sd = np.add(means,std_deviations) 
        
    ax.plot(dist_from_shore,means, linewidth = 2.5, color = "black")
    #ax.fill_between(dist_from_shore, twentyfifth_percentile, seventyfifth_percentile, facecolor='green', alpha = 0.7, interpolate=True)
    ax.fill_between(dist_from_shore, min_sd, plus_sd, facecolor='red', alpha = 0.25, interpolate=True)  
    
    ax.set_xlim(dist_from_shore[0],dist_from_shore[-1])

    ax.annotate('Standard deviation envelope', xy=(dist_from_shore[10],plus_sd[10]), xycoords='data',
                xytext=(0.1, 0.8), textcoords='axes fraction',
                size=label_size,
                # bbox=dict(boxstyle="round", fc="0.8"),
                arrowprops=dict(arrowstyle="simple",
                                fc="0.6", ec="none",
                                connectionstyle="arc3,rad=0.3"),
                )


    ax.spines['top'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2) 
    #ax.tick_params(axis='both', width=1) 

    plt.xlabel('Distance from shore (m)', fontsize = axis_size)
    plt.ylabel('Bed thickness (m)', fontsize = axis_size)

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap        
    ax.tick_params(axis='both', width=2, pad = 10)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(10)   
    
    plt.show()    


if __name__ == "__main__":

    DataDirectory = "T:\\analysis_for_papers\\Beaches\\"
    #Filename1 = "BedThickness_050.asc"
    Filename2 = "BedThickness_100.asc"
    Filename1 = "20m_bl.asc"
    axis = 1

    ElevationSwaths(DataDirectory, Filename1, axis)
    BedThickSwaths(DataDirectory, Filename2, axis)
    