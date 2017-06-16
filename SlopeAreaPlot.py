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
        LSDP.LSDMap_ChiPlotting.SlopeAreaPlot(thisPointData, DataDirectory, FigFileName=FileName, FigFormat=FigFormat, size_format=size_format, basin_key=basin_key)

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
        LSDP.LSDMap_ChiPlotting.BinnedSlopeAreaPlot(thisPointData, DataDirectory, FigFileName=FileName, FigFormat=FigFormat, size_format=size_format, basin_key=basin_key)

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

    #DataDirectory = "T:\\analysis_for_papers\\movern_testing\\"
    DataDirectory = "C:\\VagrantBoxes\\LSDTopoTools\\Topographic_projects\\Irian_jaya\\"
    DEM_prefix = "Irian_Jaya_PP"
    #DataDirectory = '/home/s0923330/DEMs_for_analysis/mid_bailey_run_10m/'
    #DEM_prefix = 'bailey_dem_10m'
    FigFormat='pdf'
    x_param='mean'
    #MakeRawSlopeAreaPlot(DataDirectory, DEM_prefix, FigFormat='png', size_format = 'ESURF')
    MakeBinnedSlopeAreaPlot(DataDirectory, DEM_prefix, FigFormat, x_param)
    MakeChannelsMap(DataDirectory, DEM_prefix, FigFormat=FigFormat)
