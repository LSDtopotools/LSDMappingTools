"""
This is a script to make slope-area plots of channels using the data
from the chi mapping tool.

It creates a separate plot for each basin, and the points are coloured by the
source node.

At the moment it uses the raw data - will also work on scripts for the binned data.

FJC
02/06/16
"""

import numpy as np
import LSDPlottingTools as LSDP
import matplotlib.pyplot as plt

def MakeRawSlopeAreaPlot(DataDirectory, DEM_prefix, FigFormat = 'show',
                         size_format = 'ESURF'):
    """
    This function makes a slope-area plot based on the raw data, using the
    channel file generated from the chi mapping tool (ends in the extension
    "_SAvertical.csv".)  It has a separate plot for each basin and colour codes
    the points by source node so different tributaries can be identified.

    Args:
        DataDirectory (str): the path to the directory with the csv file
        DEM_prefix (str): name of your DEM without extension
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).

    Returns:
        Slope-area plot for each basin

    Author: FJC
    """
    from LSDPlottingTools import LSDMap_PointTools as PointTools

    # read in the csv file
    chi_csv_fname = DataDirectory+DEM_prefix+'_SAvertical.csv'
    print("I'm reading in the csv file "+chi_csv_fname)

    # get the point data object
    thisPointData = PointTools.LSDMap_PointData(chi_csv_fname)

    # get the basin keys
    basin = thisPointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    Basin = np.asarray(basin)
    basin_keys = np.unique(Basin)
    print('There are %s basins') %(len(basin_keys))

    for basin_key in basin_keys:
        FileName = DEM_prefix+'_SA_plot_basin%s.%s' %(str(basin_key),FigFormat)
        LSDP.LSDMap_SAPlotting.SlopeAreaPlot(thisPointData, DataDirectory, FigFileName=FileName, FigFormat=FigFormat, size_format=size_format, basin_key=basin_key)

def MakeBinnedSlopeAreaPlot(DataDirectory, DEM_prefix, FigFormat = 'show',
                         size_format = 'ESURF', x_param='midpoints', y_param='mean'):
    """
    This function makes a slope-area plot based on the raw data, using the
    channel file generated from the chi mapping tool (ends in the extension
    "_SAvertical.csv".)  It has a separate plot for each basin and colour codes
    the points by source node so different tributaries can be identified.

    Args:
        DataDirectory (str): the path to the directory with the csv file
        DEM_prefix (str): name of your DEM without extension
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        x_param (str): Key for which parameter to plot on the x axis, either 'midpoints' for the midpoints of the area data (default), or 'mean' for the mean of the area data.
        y_param (str): Key for which parameter to plot on the y axis, either 'mean' for the mean of the slope data (default), or 'median', for the median of the slope data.

    Returns:
        Slope-area plot for each basin

    Author: FJC
    """
    from LSDPlottingTools import LSDMap_PointTools as PointTools

    # read in the csv file
    chi_csv_fname = DataDirectory+DEM_prefix+'_SAbinned.csv'
    print("I'm reading in the csv file "+chi_csv_fname)

    # get the point data object
    thisPointData = PointTools.LSDMap_PointData(chi_csv_fname)

    # get the basin keys
    basin = thisPointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    Basin = np.asarray(basin)
    basin_keys = np.unique(Basin)
    print('There are %s basins') %(len(basin_keys))

    for basin_key in basin_keys:
        FileName = DEM_prefix+'_SA_plot_binned_basin%s.%s' %(str(basin_key),FigFormat)
        LSDP.LSDMap_SAPlotting.BinnedSlopeAreaPlot(thisPointData, DataDirectory, FigFileName=FileName, FigFormat=FigFormat, size_format=size_format, basin_key=basin_key)

def MakeSegmentedSlopeAreaPlot(DataDirectory, DEM_prefix, FigFormat = 'show',
                         size_format = 'ESURF',basin_keys = []):
    """
    This function makes a slope-area plot based on the raw data, using the
    channel file generated from the chi mapping tool (ends in the extension
    "_SAsegmented.csv".)  It has a separate plot for each basin and colour codes
    the points by source node so different tributaries can be identified.

    Args:
        DataDirectory (str): the path to the directory with the csv file
        DEM_prefix (str): name of your DEM without extension
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        x_param (str): Key for which parameter to plot on the x axis, either 'midpoints' for the midpoints of the area data (default), or 'mean' for the mean of the area data.
        y_param (str): Key for which parameter to plot on the y axis, either 'mean' for the mean of the slope data (default), or 'median', for the median of the slope data.
        basin_key (list): A list of the basin keys to plot. If empty, plot all the basins.

    Returns:
        Slope-area plot for each basin

    Author: SMM
    """
    from LSDPlottingTools import LSDMap_PointTools as PointTools

    # read in the csv file
    chi_csv_fname = DataDirectory+DEM_prefix+'_SAsegmented.csv'
    print("I'm reading in the csv file "+chi_csv_fname)

    # get the point data object
    thisPointData = PointTools.LSDMap_PointData(chi_csv_fname)

    # get the basin keys
    basin = thisPointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    Basin = np.asarray(basin)
    these_basin_keys = np.unique(Basin)

    final_basin_keys = [] 
    # A bit of logic for checking keys
    if (len(basin_keys) == 0):
        final_basin_keys = these_basin_keys
    else:               
        for basin in basin_keys:
            if basin not in these_basin_keys:
                print("You were looking for basin "+str(basin)+ " but it isn't in the basin keys.")
            else:
                final_basin_keys.append()
            

    
    print('There are %s basins') %(len(basin_keys))        
    #basin_keys.append(0)
    # Loop through the basin keys, making a plot for each one    
    for basin_key in final_basin_keys:
        FileName = DEM_prefix+'_SA_plot_segmented_basin%s.%s' %(str(basin_key),FigFormat)
        LSDP.LSDMap_SAPlotting.SegmentedSlopeAreaPlot(thisPointData, DataDirectory, FigFileName=FileName, FigFormat=FigFormat, size_format=size_format, basin_key=basin_key)


def MakeSegmentedWithRawSlopeAreaPlot(DataDirectory, DEM_prefix, FigFormat = 'show',
                         size_format = 'ESURF', basin_keys = []):
    """
    This function makes a slope-area plot based on the raw data, using the
    channel file generated from the chi mapping tool (ends in the extension
    "_SAsegmented.csv".)  It has a separate plot for each basin and colour codes
    the points by source node so different tributaries can be identified.

    Args:
        DataDirectory (str): the path to the directory with the csv file
        DEM_prefix (str): name of your DEM without extension
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        basin_key (list): A list of the basin keys to plot. If empty, plot all the basins.

    Returns:
        Slope-area plot for each basin

    Author: SMM
    """
    from LSDPlottingTools import LSDMap_PointTools as PointTools

    # read in the csv file
    segmented_csv_fname = DataDirectory+DEM_prefix+'_SAsegmented.csv'
    print("I'm reading in the csv file "+segmented_csv_fname)
    all_csv_fname = DataDirectory+DEM_prefix+'_SAvertical.csv'

    # get the point data object
    segmentedPointData = PointTools.LSDMap_PointData(segmented_csv_fname)
    allPointData = PointTools.LSDMap_PointData(all_csv_fname)


    # get the basin keys
    basin = segmentedPointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    Basin = np.asarray(basin)
    these_basin_keys = np.unique(Basin)

    final_basin_keys = [] 
    # A bit of logic for checking keys
    if (len(basin_keys) == 0):
        final_basin_keys = these_basin_keys
    else:               
        for basin in basin_keys:
            if basin not in these_basin_keys:
                print("You were looking for basin "+str(basin)+ " but it isn't in the basin keys.")
            else:
                final_basin_keys.append(basin)
            

    
    print('There are %s basins') %(len(basin_keys))        
    #basin_keys.append(0)
    this_cmap = plt.cm.Set1
    # Loop through the basin keys, making a plot for each one    
    for basin_key in final_basin_keys:
        FileName = DEM_prefix+'_SA_plot_raw_and_segmented_basin%s.%s' %(str(basin_key),FigFormat)
        LSDP.LSDMap_SAPlotting.SegmentedWithRawSlopeAreaPlot(segmentedPointData, allPointData,
                                                              DataDirectory, FigFileName=FileName, 
                                                              FigFormat=FigFormat, size_format=size_format, 
                                                              basin_key=basin_key,cmap = this_cmap, n_colours = 10)

def MakeBinnedWithRawSlopeAreaPlot(DataDirectory, DEM_prefix, FigFormat = 'show',
                         size_format = 'ESURF', basin_keys = []):
    """
    This function makes a slope-area plot based on the raw data, using the
    channel file generated from the chi mapping tool (ends in the extension
    "_SAvertical.csv".)  It has a separate plot for each basin and colour codes
    the points by source node so different tributaries can be identified.

    Args:
        DataDirectory (str): the path to the directory with the csv file
        DEM_prefix (str): name of your DEM without extension
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        basin_key (list): A list of the basin keys to plot. If empty, plot all the basins.

    Returns:
        Slope-area plot for each basin

    Author: SMM
    """
    from LSDPlottingTools import LSDMap_PointTools as PointTools

    # read in the csv file
    binned_csv_fname = DataDirectory+DEM_prefix+'_SAbinned.csv'
    print("I'm reading in the csv file "+binned_csv_fname)
    all_csv_fname = DataDirectory+DEM_prefix+'_SAvertical.csv'

    # get the point data object
    binnedPointData = PointTools.LSDMap_PointData(binned_csv_fname)
    allPointData = PointTools.LSDMap_PointData(all_csv_fname)


    # get the basin keys
    basin = binnedPointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    Basin = np.asarray(basin)
    these_basin_keys = np.unique(Basin)

    final_basin_keys = [] 
    # A bit of logic for checking keys
    if (len(basin_keys) == 0):
        final_basin_keys = these_basin_keys
    else:               
        for basin in basin_keys:
            if basin not in these_basin_keys:
                print("You were looking for basin "+str(basin)+ " but it isn't in the basin keys.")
            else:
                final_basin_keys.append(basin)
            

    this_cmap = plt.cm.Set1
    print('There are %s basins') %(len(basin_keys))        
    #basin_keys.append(0)
    # Loop through the basin keys, making a plot for each one    
    for basin_key in final_basin_keys:
        FileName = DEM_prefix+'_SA_plot_raw_and_binned_basin%s.%s' %(str(basin_key),FigFormat)
        LSDP.LSDMap_SAPlotting.BinnedWithRawSlopeAreaPlot(binnedPointData, allPointData,
                                                              DataDirectory, FigFileName=FileName, 
                                                              FigFormat=FigFormat, size_format=size_format, 
                                                              basin_key=basin_key, n_colours = 10,
                                                              cmap = this_cmap)

def BinnedRegressionDriver(DataDirectory, DEM_prefix, basin_keys = []):
    """
    This function analyses slope-area data based on the binned data, using the
    channel file generated from the chi mapping tool (ends in the extension
    "_SAbinned.csv".)  It has a separate plot for each basin and colour codes
    the points by source node so different tributaries can be identified.

    Args:
        DataDirectory (str): the path to the directory with the csv file
        DEM_prefix (str): name of your DEM without extension
        basin_key (list): A list of the basin keys to plot. If empty, plot all the basins.

    Returns:
        Slope-area plot for each basin

    Author: SMM
    """
    from LSDPlottingTools import LSDMap_PointTools as PointTools

    # read in the csv file
    binned_csv_fname = DataDirectory+DEM_prefix+'_SAbinned.csv'
    print("I'm reading in the csv file "+binned_csv_fname)

    # get the point data object
    binnedPointData = PointTools.LSDMap_PointData(binned_csv_fname)

    # get the basin keys
    basin = binnedPointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    Basin = np.asarray(basin)
    these_basin_keys = np.unique(Basin)

    final_basin_keys = [] 
    # A bit of logic for checking keys
    if (len(basin_keys) == 0):
        final_basin_keys = these_basin_keys
    else:               
        for basin in basin_keys:
            if basin not in these_basin_keys:
                print("You were looking for basin "+str(basin)+ " but it isn't in the basin keys.")
            else:
                final_basin_keys.append(basin)
            
    # Loop through the basin keys, making a plot for each one    
    for basin_key in final_basin_keys:
        LSDP.LSDMap_SAPlotting.BinnedRegression(binnedPointData, basin_key=basin_key)


def MakeChannelsMap(DataDirectory, DEM_prefix, FigFormat = 'show',
                         size_format = 'ESURF'):
    """
    Function to make a raster map with the channels
    colour-coded by source node. The colours should match up between the slope-area plot and the map.

    Args:
        DataDirectory (str): the path to the directory with the csv file
        DEM_prefix (str): name of your DEM without extension
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).

    Returns:
        Figure with map and slope-area plot for each basin

    Author: FJC
    """
    from LSDPlottingTools import LSDMap_PointTools as PointTools

    # read in the raw csv file
    raw_csv_fname = DataDirectory+DEM_prefix+'_SAvertical.csv'
    print("I'm reading in the csv file "+raw_csv_fname)

    # read in the binned csv file
    binned_csv_fname = DataDirectory+DEM_prefix+'_SAbinned.csv'
    print("I'm reading in the csv file "+binned_csv_fname)

    # get the point data objects
    RawPointData = PointTools.LSDMap_PointData(raw_csv_fname)
    BinnedPointData = PointTools.LSDMap_PointData(binned_csv_fname)

    # get the basin keys
    basin = BinnedPointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    Basin = np.asarray(basin)
    basin_keys = np.unique(Basin)
    print('There are %s basins') %(len(basin_keys))

    # loop through each basin and make the figure
    for basin_key in basin_keys:
        # set up the figure
        FileName = DEM_prefix+'_raster_SA_basin%s.%s' %(str(basin_key),FigFormat)

        # get the channels map
        raster_fname = DataDirectory+DEM_prefix+".bil"
        hs_fname = DataDirectory+DEM_prefix+"_hs.bil"
        LSDP.LSDMap_ChiPlotting.BasicChannelPlotByBasin(raster_fname, hs_fname, raw_csv_fname, size_format=size_format, basin_key=basin, FigFileName=DataDirectory+FileName,FigFormat=FigFormat)

if __name__ == "__main__":

    DataDirectory = "T:\\analysis_for_papers\\Xian\\"
    #DataDirectory = "T:\\analysis_for_papers\\movern_testing\\"
    #DataDirectory = "C:\\VagrantBoxes\\LSDTopoTools\\Topographic_projects\\LSDTT_chi_examples\\"
    DEM_prefix = "Xian3"
    #DataDirectory = '/home/s0923330/DEMs_for_analysis/mid_bailey_run_10m/'
    #DEM_prefix = 'bailey_dem_10m'
    FigFormat='png'
    x_param='mean'
    #MakeRawSlopeAreaPlot(DataDirectory, DEM_prefix, FigFormat='png', size_format = 'ESURF')
    #MakeBinnedSlopeAreaPlot(DataDirectory, DEM_prefix, FigFormat, x_param)
    #MakeChannelsMap(DataDirectory, DEM_prefix, FigFormat=FigFormat)
    #MakeSegmentedSlopeAreaPlot(DataDirectory, DEM_prefix, FigFormat)
    
    #these_basin_keys = [0,1,2]
    these_basin_keys = []
    #MakeSegmentedWithRawSlopeAreaPlot(DataDirectory, DEM_prefix, FigFormat, basin_keys = these_basin_keys)
    #MakeBinnedWithRawSlopeAreaPlot(DataDirectory, DEM_prefix, FigFormat, basin_keys = these_basin_keys)    
    #BinnedRegressionDriver(DataDirectory, DEM_prefix, basin_keys = these_basin_keys)
    LSDP.SAPlotDriver(DataDirectory, DEM_prefix, FigFormat = FigFormat,
                      show_raw = True, show_segments = False,
                      cmap = plt.cm.Set1, n_colours = 5,
                      basin_keys = these_basin_keys)
    
    #adict = LSDP.BinnedRegressionDriver(DataDirectory, DEM_prefix, basin_keys = these_basin_keys)
    
    from LSDPlottingTools import LSDMap_MOverNPlotting as MN
    MN.CompareChiAndSAMOverN(DataDirectory, DEM_prefix, basin_list=these_basin_keys, start_movern=0.1, d_movern=0.1, n_movern=8)
    
    