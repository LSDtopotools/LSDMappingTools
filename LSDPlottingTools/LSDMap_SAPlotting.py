## LSDMap_SAPlotting.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## These functions are tools to deal slope-area data
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## SMM
## 14/12/2016
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import matplotlib.pyplot as plt
#from cycler import cycler
from matplotlib import rcParams
#import LSDPlottingTools.LSDMap_GDALIO as LSDMap_IO
#import LSDMap_BasicManipulation as LSDMap_BM
#import LSDMap_OSystemTools as LSDOst
#import LSDPlottingTools.LSDMap_BasicPlotting as LSDMap_BP
#import LSDPlottingTools.LSDMap_PointTools as LSDMap_PD
#import LSDPlottingTools.LSDMap_BasicManipulation as LSDMap_BM
import LSDPlottingTools.statsutilities as LSDStats

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Slope-area functions
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def SlopeAreaPlot(PointData, DataDirectory, FigFileName = 'Image.pdf',
                       FigFormat = 'show',
                       size_format = "ESURF",
                       basin_key = '0'):
    """
    This function makes a slope-area plot from the chi mapping tool.
    This is a work in progress! Need to think about what we want to visualise:
    I think maybe a separate plot for each basin, and then colour-code by source ID.

    Args:
        PointData : LSDPointData object produced from the csv file with chi, chi slope, etc information. This file is produced by the chi_mapping_tool. It should have the extension "_SAvertical.csv"
        FigFileName (str): The name of the figure file
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        basin_key (int): the ID of the basin to make the plot for.

    Returns:
         Does not return anything but makes a plot.

    Author: FJC
    """
    import matplotlib.colors as colors

    label_size = 10

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])

    # Get the slope, drainage area, basin ID and source ID
    slope = PointData.QueryData('slope')
    slope = [float(x) for x in slope]
    area = PointData.QueryData('drainage area')
    area = [float(x) for x in area]
    basin = PointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    source = PointData.QueryData('source_key')
    source = [int(x) for x in source]

    # need to convert everything into arrays so we can mask different basins
    Slope = np.asarray(slope)
    Area = np.asarray(area)
    Basin = np.asarray(basin)
    Source = np.asarray(source)

    # mask to just get the data for the basin of interest
    m = np.ma.masked_where(Basin!=basin_key, Basin)
    maskSlope = np.ma.masked_where(np.ma.getmask(m), Slope)
    maskArea = np.ma.masked_where(np.ma.getmask(m), Area)
    maskSource = np.ma.masked_where(np.ma.getmask(m), Source)

    #colour by source - this is the same as the script to colour channels over a raster,
    # (BasicChannelPlotGridPlotCategories) so that the colour scheme should match
    # make a color map of fixed colors
    NUM_COLORS = 15
    this_cmap = plt.cm.Set1
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    plt.cm.ScalarMappable(norm=cNorm, cmap=this_cmap)
    channel_data = [x % NUM_COLORS for x in maskSource]

    # now make the slope area plot. Need to add a lot more here but just to test for now.
    ax.scatter(maskArea,maskSlope,c=channel_data,cmap=this_cmap,s=10,marker="+",lw=1)
    ax.set_xlabel('Drainage area (m$^2$)')
    ax.set_ylabel('Slope (m/m)')

    # log
    ax.set_xscale('log')
    ax.set_yscale('log')

    # set axis limits
    x_pad = 1000
    y_pad = 0.00001
    ax.set_ylim(np.min(maskSlope)-y_pad,0)
    ax.set_xlim(np.min(maskArea)-x_pad,np.max(maskArea)+y_pad)

    # return or show the figure
    print("The figure format is: " + FigFormat)
    if FigFormat == 'show':
        plt.show()
    elif FigFormat == 'return':
        return fig
    else:
        save_fmt = FigFormat
        plt.savefig(DataDirectory+FigFileName,format=save_fmt,dpi=500)
        fig.clf()

def BinnedSlopeAreaPlot(PointData, DataDirectory, FigFileName = 'Image.pdf',
                       FigFormat = 'show',
                       size_format = "ESURF",
                       basin_key = '0', x_param = 'midpoints', y_param = 'mean'):
    """
    This function makes a slope-area plot from the chi mapping tool using the binned data.

    Args:
        PointData : LSDPointData object produced from the csv file with chi, chi slope, etc information. This file is produced by the chi_mapping_tool. It should have the extension "_SAbinned.csv"
        FigFileName (str): The name of the figure file
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        basin_key (int): the ID of the basin to make the plot for.
        x_param (str): Key for which parameter to plot on the x axis, either 'midpoints' for the midpoints of the area data (default), or 'mean' for the mean of the area data.
        y_param (str): Key for which parameter to plot on the y axis, either 'mean' for the mean of the slope data (default), or 'median', for the median of the slope data.

    Returns:
         Does not return anything but makes a plot.

    Author: FJC
    """
    import matplotlib.colors as colors
    import matplotlib.ticker

    label_size = 10

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])

    # Get the slope, drainage area, basin ID and source ID
    mean_log_S = PointData.QueryData('mean_log_S')
    mean_log_S = [float(10**x) for x in mean_log_S]
    median_log_S = PointData.QueryData('median_log_S')
    median_log_S = [float(10**x) for x in median_log_S]
    mean_log_A = PointData.QueryData('mean_log_A')
    mean_log_A = [float(10**x) for x in mean_log_A]
    midpoints_A = PointData.QueryData('midpoints_log_A')
    midpoints_A = [float(10**x) for x in midpoints_A]
    basin = PointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    source = PointData.QueryData('source_key')
    source = [int(x) for x in source]

    # get the errors
    log_S_sterr = PointData.QueryData('logS_stdErr')
    log_S_sterr = [float(10**x) for x in log_S_sterr]
    log_A_sterr = PointData.QueryData('logA_stdErr')
    log_A_sterr = [float(10**x) for x in log_A_sterr]

    # need to convert everything into arrays so we can mask different basins
    MeanLogSlope = np.asarray(mean_log_S)
    MedianLogSlope = np.asarray(median_log_S)
    MeanLogArea = np.asarray(mean_log_A)
    MidpointsArea = np.asarray(midpoints_A)
    SlopeError = np.asarray(log_S_sterr)
    AreaError = np.asarray(log_A_sterr)
    Basin = np.asarray(basin)
    Source = np.asarray(source)

    # mask to just get the data for the basin of interest
    m = np.ma.masked_where(Basin!=basin_key, Basin)
    MeanLogSlope = np.ma.masked_where(np.ma.getmask(m), MeanLogSlope)
    MedianLogSlope = np.ma.masked_where(np.ma.getmask(m), MedianLogSlope)
    MeanLogArea = np.ma.masked_where(np.ma.getmask(m), MeanLogArea)
    MidpointsArea = np.ma.masked_where(np.ma.getmask(m), MidpointsArea)
    SlopeError = np.ma.masked_where(np.ma.getmask(m), SlopeError)
    AreaError = np.ma.masked_where(np.ma.getmask(m), AreaError)
    maskSource = np.ma.masked_where(np.ma.getmask(m), Source)

    #colour by source - this is the same as the script to colour channels over a raster,
    # (BasicChannelPlotGridPlotCategories) so that the colour scheme should match
    # make a color map of fixed colors
    NUM_COLORS = 15
    this_cmap = plt.cm.Set1
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    plt.cm.ScalarMappable(norm=cNorm, cmap=this_cmap)
    channel_data = [x % NUM_COLORS for x in maskSource]

    # now make the slope area plot. Need to add a lot more here but just to test for now.
    if x_param == 'mean' and y_param == 'mean':
        #plt.errorbar(MeanLogArea,MeanLogSlope,xerr=AreaError,yerr=SlopeError,fmt='o',ms=1,ecolor='k')
        ax.scatter(MeanLogArea,MeanLogSlope,c=channel_data,cmap=this_cmap,s=10,marker="o",lw=0.5,edgecolors='k',zorder=100)
    elif x_param == 'mean' and y_param == 'median':
        #plt.errorbar(MeanLogArea,MedianLogSlope,xerr=AreaError,yerr=SlopeError,fmt='o',ms=1,ecolor='k')
        ax.scatter(MeanLogArea,MedianLogSlope,c=channel_data,cmap=this_cmap,s=10,marker="o",lw=0.5,edgecolors='k',zorder=100)
    elif x_param == 'midpoints' and y_param == 'median':
        #plt.errorbar(MidpointsArea,MedianLogSlope,xerr=AreaError,yerr=SlopeError,fmt='o',ms=1,ecolor='k')
        ax.scatter(MidpointsArea,MedianLogSlope,c=channel_data,cmap=this_cmap,s=10,marker="o",lw=0.5,edgecolors='k',zorder=100)
    else:
        #plt.errorbar(MidpointsArea,MeanLogSlope,xerr=AreaError,yerr=SlopeError,fmt='o',ms=1,ecolor='k')
        ax.scatter(MidpointsArea,MeanLogSlope,c=channel_data,cmap=this_cmap,s=10,marker="o",lw=0.5,edgecolors='k',zorder=100)

    ax.set_xlabel('Drainage area (m$^2$)')
    ax.set_ylabel('Slope (m/m)')

    # log
    ax.set_xscale('log')
    ax.set_yscale('log')

    # set axis limits
    x_pad = 1000
    y_pad = 0.00001
    ax.set_ylim(np.min(MeanLogSlope)-y_pad,0)
    ax.set_xlim(np.min(MeanLogArea)-x_pad,np.max(MeanLogArea)+y_pad)

    # return or show the figure
    print("The figure format is: " + FigFormat)
    if FigFormat == 'show':
        plt.show()
    elif FigFormat == 'return':
        return fig # return the axes object so can make nice subplots with the other plotting tools?
    else:
        save_fmt = FigFormat
        plt.savefig(DataDirectory+FigFileName,format=save_fmt,dpi=500)
        fig.clf()

def SegmentedSlopeAreaPlot(PointData, DataDirectory, FigFileName = 'Image.pdf',
                       FigFormat = 'show',
                       size_format = "ESURF",
                       basin_key = '0'):
    """
    This function makes a slope-area plot from the chi mapping tool using the binned data.

    Args:
        PointData : LSDPointData object produced from the csv file with chi, chi slope, etc information. This file is produced by the chi_mapping_tool. It should have the extension "_SAbinned.csv"
        FigFileName (str): The name of the figure file
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        basin_key (int): the ID of the basin to make the plot for.

    Returns:
         Does not return anything but makes a plot.

    Author: SMM
    """
    import matplotlib.colors as colors
    import matplotlib.ticker

    label_size = 10

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])

    # Get the slope, drainage area, basin ID and source ID
    median_log_S = PointData.QueryData('median_log_S')
    median_log_S = [float(10**x) for x in median_log_S]
    median_log_A = PointData.QueryData('median_log_A')
    median_log_A = [float(10**x) for x in median_log_A]
    fitted_log_S = PointData.QueryData('segmented_log_S')
    fitted_log_S = [float(10**x) for x in fitted_log_S]    
    basin = PointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    segment_number = PointData.QueryData('segment_number')
    segment_number = [int(x) for x in segment_number]

    # get the errors
    firstquartile= PointData.QueryData('logS_FirstQuartile')
    firstquartile = [float(10**x) for x in firstquartile]
    thirdquartile= PointData.QueryData('logS_ThirdQuartile')
    thirdquartile = [float(10**x) for x in thirdquartile]
    
    #print("Size of quartiles: "+ str( len(firstquartile))+ " "+str( len(thirdquartile)))
    
    # need to convert everything into arrays so we can mask different basins
    MedianLogSlope = np.asarray(median_log_S)
    MedianLogArea = np.asarray(median_log_A)
    FirstQuartile = np.asarray(firstquartile)
    ThirdQuartile = np.asarray(thirdquartile)
    FittedLogS = np.asarray(fitted_log_S)
    Basin = np.asarray(basin)
    segment_number = np.asarray(segment_number)
    
    # Get the errors
    yerr_down= np.subtract(ThirdQuartile,median_log_S)
    yerr_up= np.subtract(median_log_S,FirstQuartile)    

    # mask to just get the data for the basin of interest
    m = np.ma.masked_where(Basin!=basin_key, Basin)
    MedianLogSlope = np.ma.masked_where(np.ma.getmask(m), MedianLogSlope)
    MedianLogArea = np.ma.masked_where(np.ma.getmask(m), MedianLogArea)
    FittedLogS = np.ma.masked_where(np.ma.getmask(m), FittedLogS)
    mask_segment_number = np.ma.masked_where(np.ma.getmask(m), segment_number)
    
    # now make the slope area plot. Need to add a lot more here but just to test for now.
    plt.errorbar(MedianLogArea,MedianLogSlope,yerr=[yerr_up,yerr_down],fmt='o',ms=1,ecolor='k')
    ax.scatter(MedianLogArea,MedianLogSlope,c="b",s=10,marker="o",lw=0.5,edgecolors='k',zorder=100)

    # now get the segments
    segments = np.unique(segment_number)  
    #n_segments = len(segments)
    #print("The unique segment numbers are: ")
    #print(segments)
    #print("There are: "+str(n_segments)+" of them")
    
    # Mask the data of the segments sequentially
    for segment in segments:
    # mask to just get the data for the basin of interest
        m = np.ma.masked_where(mask_segment_number!=segment, mask_segment_number)
        MedianLogSlope = np.ma.masked_where(np.ma.getmask(m), MedianLogSlope)
        SegmentMedianLogArea = np.ma.masked_where(np.ma.getmask(m), MedianLogArea) 
        SegmentFittedLogS = np.ma.masked_where(np.ma.getmask(m), FittedLogS) 
        ax.plot(SegmentMedianLogArea,SegmentFittedLogS)
    

    ax.set_xlabel('Drainage area (m$^2$)')
    ax.set_ylabel('Slope (m/m)')

    # log
    ax.set_xscale('log')
    ax.set_yscale('log')

    # set axis limits
    #x_pad = 1000
    #y_pad = 0.0000001
    #ax.set_ylim(np.min(MedianLogSlope)-y_pad,0)
    #ax.set_xlim(np.min(MeanLogArea)-x_pad,np.max(MeanLogArea)+y_pad)

    # return or show the figure
    print("The figure format is: " + FigFormat)
    if FigFormat == 'show':
        plt.show()
    elif FigFormat == 'return':
        return fig # return the axes object so can make nice subplots with the other plotting tools?
    else:
        save_fmt = FigFormat
        plt.savefig(DataDirectory+FigFileName,format=save_fmt,dpi=500)
        fig.clf()

def SegmentedWithRawSlopeAreaPlot(PointData, RawPointData, DataDirectory, FigFileName = 'Image.pdf',
                       FigFormat = 'show',size_format = "ESURF",basin_key = '0',
                       cmap = plt.cm.Set1, n_colours = 5):
    """
    This function makes a slope-area plot from the chi mapping tool using the binned data.
    It plots the main stem only, but has the raw data as semitrasparent in the background.
    It also plots the segments on the main stem, determined by the segmentation algorithm. 

    Args:
        PointData : LSDPointData object produced from the csv file with binned and segmented slope area data. Produced by chi mapping tool. It should have the extension "_SAsegmented.csv"
        RawPointData: LSDPointData object produced from the csv file with binned and segmented slope area data. Produced by chi mapping tool. It should have the extension "_SAvertical.csv"
        FigFileName (str): The name of the figure file
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        basin_key (int): the ID of the basin to make the plot for.
        colourmap (string or colormap object): The colourmap
        n_colour (int): The number of colours
        
    Returns:
         Does not return anything but makes a plot.

    Author: SMM
    """
    import matplotlib.colors as colors
    import matplotlib.ticker

    label_size = 10

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])
    
    # Get the raw point data    
    Raw_log_S =  RawPointData.QueryData('slope')
    Raw_log_S = [float(x) for x in Raw_log_S]
    Raw_log_A =  RawPointData.QueryData('drainage area')
    Raw_log_A = [float(x) for x in Raw_log_A]
    
    Raw_basin = RawPointData.QueryData('basin_key')
    Raw_basin = [int(x) for x in Raw_basin]
    Raw_source = RawPointData.QueryData('source_key')
    Raw_source = [int(x) for x in Raw_source]   
    
    # Convert to numpy array
    RS = np.asarray(Raw_log_S)
    RA = np.asarray(Raw_log_A)
    RB = np.asarray(Raw_basin)
    RSource = np.asarray(Raw_source)
    
    # Find the source for the basin using brute force. 
    found_basin = False
    counter = 0
    this_source = 0
    while not found_basin:
        if RB[counter] == float(basin_key):
            this_source = RSource[counter]
            found_basin = True
        else:
            counter+=1
            
    # Now mask the data to the correct source (this is the main stem)
    m = np.ma.masked_where(RSource!=this_source, RSource)
    RawSlope = np.ma.masked_where(np.ma.getmask(m), RS)
    RawArea = np.ma.masked_where(np.ma.getmask(m), RA)
    
    # Plot the raw data    
    ax.scatter(RawArea,RawSlope,c="k",s=4,marker="+",lw=0.5,edgecolors='k',zorder=-10,alpha = 0.3)
         
    

    # Get the slope, drainage area, basin ID and source ID
    median_log_S = PointData.QueryData('median_log_S')
    median_log_S = [float(10**x) for x in median_log_S]
    median_log_A = PointData.QueryData('median_log_A')
    median_log_A = [float(10**x) for x in median_log_A]
    fitted_log_S = PointData.QueryData('segmented_log_S')
    fitted_log_S = [float(10**x) for x in fitted_log_S]    
    basin = PointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    segment_number = PointData.QueryData('segment_number')
    segment_number = [int(x) for x in segment_number]

    # get the errors
    firstquartile= PointData.QueryData('logS_FirstQuartile')
    firstquartile = [float(10**x) for x in firstquartile]
    thirdquartile= PointData.QueryData('logS_ThirdQuartile')
    thirdquartile = [float(10**x) for x in thirdquartile]
    
    #print("Size of quartiles: "+ str( len(firstquartile))+ " "+str( len(thirdquartile)))
    
    # need to convert everything into arrays so we can mask different basins
    MedianLogSlope = np.asarray(median_log_S)
    MedianLogArea = np.asarray(median_log_A)
    FirstQuartile = np.asarray(firstquartile)
    ThirdQuartile = np.asarray(thirdquartile)
    FittedLogS = np.asarray(fitted_log_S)
    Basin = np.asarray(basin)
    segment_number = np.asarray(segment_number)
    
    # Get the errors
    yerr_down= np.subtract(ThirdQuartile,median_log_S)
    yerr_up= np.subtract(median_log_S,FirstQuartile)    

    # mask to just get the data for the basin of interest
    m = np.ma.masked_where(Basin!=basin_key, Basin)
    MedianLogSlope = np.ma.masked_where(np.ma.getmask(m), MedianLogSlope)
    MedianLogArea = np.ma.masked_where(np.ma.getmask(m), MedianLogArea)
    FittedLogS = np.ma.masked_where(np.ma.getmask(m), FittedLogS)
    mask_segment_number = np.ma.masked_where(np.ma.getmask(m), segment_number)

    # make a color map of fixed colors
    NUM_COLORS = n_colours
    # First we set the colourmap
    this_cmap = cmap
    # then we use a normalization to map the colours between 0 and NUM_COLORS-1
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    # Now we make a scalar map. This is used to convert values in your dataset
    # to values between 0 and 1 that can be called to convert to rgba
    scalarMap = plt.cm.ScalarMappable(norm=cNorm, cmap=this_cmap)
    # If you want RGBA from this you use:  rgba_color = scalarMap.to_rgba(this_data)    
  
    segments = np.unique(segment_number)   
  
    # Mask the data of the segments sequentially
    for segment in segments:
    # mask to just get the data for the basin of interest
        m = np.ma.masked_where(mask_segment_number!=segment, mask_segment_number)
        SegmentMedianLogSlope = np.ma.masked_where(np.ma.getmask(m), MedianLogSlope)
        SegmentMedianLogArea = np.ma.masked_where(np.ma.getmask(m), MedianLogArea) 
        SegmentFittedLogS = np.ma.masked_where(np.ma.getmask(m), FittedLogS) 
        
        # Now add the colours for the segments
        tps_color = scalarMap.to_rgba(segment) 
        ax.plot(SegmentMedianLogArea,SegmentFittedLogS,c=tps_color,zorder=-10)
        ax.scatter(SegmentMedianLogArea,SegmentMedianLogSlope,c=tps_color,s=10,marker="o",lw=0.5,edgecolors='k',zorder=100)
        plt.errorbar(SegmentMedianLogArea,SegmentMedianLogSlope,yerr=[yerr_up,yerr_down],fmt='o',ms=1,ecolor=tps_color,zorder=0)
    

    ax.set_xlabel('Drainage area (m$^2$)')
    ax.set_ylabel('Slope (m/m)')

    # log
    ax.set_xscale('log')
    ax.set_yscale('log')

    # set axis limits
    #x_pad = 1000
    #y_pad = 0.0000001
    #ax.set_ylim(np.min(MedianLogSlope)-y_pad,0)
    #ax.set_xlim(np.min(MeanLogArea)-x_pad,np.max(MeanLogArea)+y_pad)

    # return or show the figure
    print("The figure format is: " + FigFormat)
    if FigFormat == 'show':
        plt.show()
    elif FigFormat == 'return':
        return fig # return the axes object so can make nice subplots with the other plotting tools?
    else:
        save_fmt = FigFormat
        plt.savefig(DataDirectory+FigFileName,format=save_fmt,dpi=500)
        fig.clf()

def BinnedWithRawSlopeAreaPlot(BinnedPointData, RawPointData, DataDirectory, FigFileName = 'Image.pdf',
                       FigFormat = 'show',size_format = "ESURF",basin_key = '0',
                       cmap = plt.cm.Set1, n_colours = 5):
    """
    This function makes a slope-area plot from the chi mapping tool using the binned data.
    It plots all the sources and all the raw data as semitransparent background.

    Args:
        PointData : LSDPointData object produced from the csv file with binned and segmented slope area data. Produced by chi mapping tool. It should have the extension "_SAsegmented.csv"
        RawPointData: LSDPointData object produced from the csv file with binned and segmented slope area data. Produced by chi mapping tool. It should have the extension "_SAvertical.csv"
        FigFileName (str): The name of the figure file
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        basin_key (int): the ID of the basin to make the plot for.
        colourmap (string or colormap object): The colourmap
        n_colour (int): The number of colours

    Returns:
         Does not return anything but makes a plot.

    Author: SMM
    """
    import matplotlib.colors as colors
    import matplotlib.ticker

    label_size = 10

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])
    
    # Get the raw point data    
    Raw_log_S =  RawPointData.QueryData('slope')
    Raw_log_S = [float(x) for x in Raw_log_S]
    Raw_log_A =  RawPointData.QueryData('drainage area')
    Raw_log_A = [float(x) for x in Raw_log_A]
    
    Raw_basin = RawPointData.QueryData('basin_key')
    Raw_basin = [int(x) for x in Raw_basin]
    Raw_source = RawPointData.QueryData('source_key')
    Raw_source = [int(x) for x in Raw_source]   
    
    # Convert to numpy array
    RS = np.asarray(Raw_log_S)
    RA = np.asarray(Raw_log_A)
    RB = np.asarray(Raw_basin)
    RSource = np.asarray(Raw_source)    
        
    # Now mask the data to the correct source (this is the main stem)
    this_basin = float(basin_key)    
    m = np.ma.masked_where(RB!=this_basin,RB)
    RawSlope = np.ma.masked_where(np.ma.getmask(m), RS)
    RawArea = np.ma.masked_where(np.ma.getmask(m), RA)
    RawSource = np.ma.masked_where(np.ma.getmask(m), RSource)
             
    # Get the slope, drainage area, basin ID and source ID
    median_log_S = BinnedPointData.QueryData('median_log_S')
    median_log_S = [float(10**x) for x in median_log_S]
    median_log_A = BinnedPointData.QueryData('median_log_A')
    median_log_A = [float(10**x) for x in median_log_A] 
    basin = BinnedPointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    source_number = BinnedPointData.QueryData('source_key')
    source_number = [int(x) for x in source_number]
    
    # get the errors
    firstquartile= BinnedPointData.QueryData('logS_FirstQuartile')
    firstquartile = [float(10**x) for x in firstquartile]
    thirdquartile= BinnedPointData.QueryData('logS_ThirdQuartile')
    thirdquartile = [float(10**x) for x in thirdquartile]

    print("The lengths of the data vectors:")
    print(len(median_log_S))
    print(len(median_log_A))
    print(len(basin))
    print(len(source_number))    
    print("Size of quartiles: "+ str( len(firstquartile))+ " "+str( len(thirdquartile)))
    
    # need to convert everything into arrays so we can mask different basins
    MedianLogSlope = np.asarray(median_log_S)
    MedianLogArea = np.asarray(median_log_A)
    FirstQuartile = np.asarray(firstquartile)
    ThirdQuartile = np.asarray(thirdquartile)
    Basin = np.asarray(basin)
    SourceNumber = np.asarray(source_number)
    
    # Get the errors
    yerr_down= np.subtract(ThirdQuartile,median_log_S)
    yerr_up= np.subtract(median_log_S,FirstQuartile)    

    # mask to just get the data for the basin of interest
    m = np.ma.masked_where(Basin!=basin_key, Basin)
    MedianLogSlope = np.ma.masked_where(np.ma.getmask(m), MedianLogSlope)
    MedianLogArea = np.ma.masked_where(np.ma.getmask(m), MedianLogArea)
    SourceNumber = np.ma.masked_where(np.ma.getmask(m), SourceNumber)

    # make a color map of fixed colors
    NUM_COLORS = n_colours
    # First we set the colourmap
    #this_cmap = plt.cm.Set1
    this_cmap = cmap
    # then we use a normalization to map the colours between 0 and NUM_COLORS-1
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    # Now we make a scalar map. This is used to convert values in your dataset
    # to values between 0 and 1 that can be called to convert to rgba
    scalarMap = plt.cm.ScalarMappable(norm=cNorm, cmap=this_cmap)
    # If you want RGBA from this you use:  rgba_color = scalarMap.to_rgba(this_data)    

    # now get the sources
    sources = np.unique(SourceNumber)  

    
    # Mask the data of the segments sequentially
    for idx,source in enumerate(sources):
    # mask to just get the data for the basin of interest
        m = np.ma.masked_where(SourceNumber!=source, SourceNumber)
        SourceMedianLogSlope = np.ma.masked_where(np.ma.getmask(m), MedianLogSlope)
        SourceMedianLogArea = np.ma.masked_where(np.ma.getmask(m), MedianLogArea) 
        
        # Now add the colours for the segments
        source_colour = idx%n_colours
        tps_color = scalarMap.to_rgba(source_colour) 
        ax.scatter(SourceMedianLogArea,SourceMedianLogSlope,c=tps_color,s=10,marker="o",lw=0.5,edgecolors='k',zorder=100)
        plt.errorbar(SourceMedianLogArea,SourceMedianLogSlope,yerr=[yerr_up,yerr_down],fmt='o',ms=1,ecolor=tps_color,zorder=0)
    
        # Plot the raw data 
        m2 = np.ma.masked_where(RawSource!=source, RawSource)
        SourceRawMedianS = np.ma.masked_where(np.ma.getmask(m2), RawSlope)
        SourceRawMedianA = np.ma.masked_where(np.ma.getmask(m2), RawArea)
        ax.scatter(SourceRawMedianA,SourceRawMedianS,c=tps_color,s=4,marker="+",lw=0.5,edgecolors='k',zorder=-10,alpha = 0.3)
    

    ax.set_xlabel('Drainage area (m$^2$)')
    ax.set_ylabel('Slope (m/m)')

    # log
    ax.set_xscale('log')
    ax.set_yscale('log')

    # set axis limits
    #x_pad = 1000
    #y_pad = 0.0000001
    #ax.set_ylim(np.min(MedianLogSlope)-y_pad,0)
    #ax.set_xlim(np.min(MeanLogArea)-x_pad,np.max(MeanLogArea)+y_pad)

    # return or show the figure
    print("The figure format is: " + FigFormat)
    if FigFormat == 'show':
        plt.show()
    elif FigFormat == 'return':
        return fig # return the axes object so can make nice subplots with the other plotting tools?
    else:
        save_fmt = FigFormat
        plt.savefig(DataDirectory+FigFileName,format=save_fmt,dpi=500)
        fig.clf()





def BinnedRegression(BinnedPointData, basin_key):
    """
    This function makes a slope-area plot from the chi mapping tool using the binned data.
    It plots all the sources and all the raw data as semitransparent background.

    Args:
        PointData : LSDPointData object produced from the csv file with binned and segmented slope area data. Produced by chi mapping tool. It should have the extension "_SAsegmented.csv"
        FigFileName (str): The name of the figure file
        basin_key (int): the ID of the basin to make the plot for.#
        
    Returns:
         Does not return anything but makes a plot.

    Author: SMM
    """
        
    # Get the slope, drainage area, basin ID and source ID
    median_log_S = BinnedPointData.QueryData('median_log_S')
    median_log_S = [float(x) for x in median_log_S]
    median_log_A = BinnedPointData.QueryData('median_log_A')
    median_log_A = [float(x) for x in median_log_A] 
    basin = BinnedPointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    source_number = BinnedPointData.QueryData('source_key')
    source_number = [int(x) for x in source_number]
    
    # get the errors
    firstquartile= BinnedPointData.QueryData('logS_FirstQuartile')
    firstquartile = [float(x) for x in firstquartile]
    thirdquartile= BinnedPointData.QueryData('logS_ThirdQuartile')
    thirdquartile = [float(x) for x in thirdquartile]

    #print("The lengths of the data vectors:")
    #print(len(median_log_S))
    #print(len(median_log_A))
    #print(len(basin))
    #print(len(source_number))    
    #print("Size of quartiles: "+ str( len(firstquartile))+ " "+str( len(thirdquartile)))
    
    # need to convert everything into arrays so we can mask different basins
    MedianLogSlope = np.asarray(median_log_S)
    MedianLogArea = np.asarray(median_log_A)
    FirstQuartile = np.asarray(firstquartile)
    ThirdQuartile = np.asarray(thirdquartile)
    Basin = np.asarray(basin)
    SourceNumber = np.asarray(source_number)
    
    # Get the errors
    yerr_down= np.subtract(ThirdQuartile,median_log_S)
    yerr_up= np.subtract(median_log_S,FirstQuartile)    

    # mask to just get the data for the basin of interest
    print("The basin key is: "+str(basin_key))
    m = np.ma.masked_where(Basin!=basin_key, Basin)
    MedianLogSlope = np.ma.masked_where(np.ma.getmask(m), MedianLogSlope)
    MedianLogArea = np.ma.masked_where(np.ma.getmask(m), MedianLogArea)
    SourceNumber = np.ma.masked_where(np.ma.getmask(m), SourceNumber)
    
    # get the compressed data
    SlopeCompressed = np.ma.compressed(MedianLogSlope)
    AreaCompressed = np.ma.compressed(MedianLogArea)  
    SourceCompressed = np.ma.compressed(SourceNumber)
    
    # get main stem source
    mainstem_source = SourceCompressed[0]
    
    m2 = np.ma.masked_where(SourceCompressed!=mainstem_source, SourceCompressed)
    MS_Slope =  np.ma.masked_where(np.ma.getmask(m2), SlopeCompressed)
    MS_Area = np.ma.masked_where(np.ma.getmask(m2), AreaCompressed)
    
    MSSlopeCompressed = np.ma.compressed(MS_Slope)
    MSAreaCompressed = np.ma.compressed(MS_Area)     

    # get the regression from the main stem
    [MSresiduals,m,b,r,pvalue,stderr]= LSDStats.linregress_residuals(MSAreaCompressed,MSSlopeCompressed)   
    print("slope of mainstem regression is: "+str(m))

    # see if there are any outlying residuals
    [new_x,new_y, is_outlier_vec, m,b]= LSDStats.remove_outlying_residuals(MSAreaCompressed,MSSlopeCompressed,MSresiduals)
    print("Removed outlier slope from mainstem data: " +str(m))

    # get the regression from all the data
    [residuals,m,b,r,pvalue,stderr]= LSDStats.linregress_residuals(AreaCompressed,SlopeCompressed)   
    print("slope of all data is: "+str(m))
    
    # see if there are any outlying residuals
    [new_x,new_y, is_outlier_vec, m,b]= LSDStats.remove_outlying_residuals(AreaCompressed,SlopeCompressed,residuals)
    print("Removed outlier slope from all data: " +str(m))
