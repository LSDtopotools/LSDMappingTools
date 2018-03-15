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
import pandas as pd
import os
from scipy import stats
#import LSDPlottingTools.LSDMap_GDALIO as LSDMap_IO
#import LSDMap_BasicManipulation as LSDMap_BM
#import LSDMap_OSystemTools as LSDOst
#import LSDPlottingTools.LSDMap_BasicPlotting as LSDMap_BP
#import LSDPlottingTools.LSDMap_PointTools as LSDMap_PD
#import LSDPlottingTools.LSDMap_BasicManipulation as LSDMap_BM
import LSDPlottingTools.statsutilities as LSDStats
from LSDMapFigure import PlottingHelpers as Helper


##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Regressions
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def LinearRegressionRawData(DataDirectory, DEM_prefix, basin_list=[],parallel=False):
    """
    This function performs a linear regression on all of the slope-area data.
    It returns a dataframe with the linear regression info for each basin key.

    Args:
        DataDirectory (str): the data directory
        DEM_prefix (str): the prefix of the DEM_prefix
        basin_list: a list of the basins to analyse, default = empty (all basins)

    Returns:
        pandas dataframe with the linear regression info

    Author: SMM and FJC
    """
    # read in binned data
    if not parallel:
      df = Helper.ReadRawSAData(DataDirectory, DEM_prefix)
    else:
      df = Helper.AppendRawSAData(DataDirectory, DEM_prefix)

    # get a list of the basins if needed
    if basin_list == []:
        print ("You didn't give me a basin list so I will analyse all the basins")
        basin_list = df['basin_key'].unique()

    # now do a linear regression for each basin
    columns = ['basin_key', 'regression_slope', 'std_err', 'R2', 'p_value']
    OutDF = pd.DataFrame(columns=columns)

    for basin_key in basin_list:
        df_slope = (df[df['basin_key']==basin_key]['slope']).values
        df_area = (df[df['basin_key']==basin_key]['drainage_area']).values

        logS = np.log10(df_slope[df_area!=0])
        logA = np.log10(df_area[df_area!=0])

        slope, intercept, r_value, p_value, std_err = stats.linregress(logA,logS)

        #print("Slope: " +str(slope)+ " std_err: "+str(std_err)+ " R2 is: " + str(r_value**2) + " p value is: " + str(p_value))
        this_key = int(basin_key)
        this_row = [this_key,abs(slope),std_err,r_value**2,p_value]
        OutDF.loc[basin_key] = this_row

    return OutDF

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Regressions
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def LinearRegressionRawDataByChannel(DataDirectory, DEM_prefix, basin_list=[]):
    """
    This function performs a linear regression on the raw slope-area data separated.
    by channel.
    It returns a dataframe with the linear regression info for each basin key.

    Args:
        DataDirectory (str): the data directory
        DEM_prefix (str): the prefix of the DEM_prefix
        basin_list: a list of the basins to analyse, default = empty (all basins)

    Returns:
        pandas dataframe with the linear regression info

    Author: SMM and FJC
    """
    # read in binned data


    df = Helper.ReadRawSAData(DataDirectory, DEM_prefix)

    # get a list of the basins if needed
    if basin_list == []:
        print ("You didn't give me a basin list so I will analyse all the basins")
        basin_list = df['basin_key'].unique()
        print (basin_list)


    # now do a linear regression for each basin
    columns = ['basin_key', 'source_key', 'regression_slope', 'std_err', 'R2', 'p_value']
    OutDF = pd.DataFrame(columns=columns)
    counter = 0

    for basin_key in basin_list:
        this_df = df[df['basin_key'] == basin_key]
        # get the sources for this basin
        these_sources = this_df['source_key'].unique()

        for source in these_sources:
            #mask the dataframe for this source
            df_slope = (df[df['source_key']==source]['slope']).values
            df_area = (df[df['source_key']==source]['drainage_area']).values

            logS = np.log10(df_slope[df_area!=0])
            logA = np.log10(df_area[df_area!=0])

            slope, intercept, r_value, p_value, std_err = stats.linregress(logA,logS)

            #print("Slope: " +str(slope)+ " std_err: "+str(std_err)+ " R2 is: " + str(r_value**2) + " p value is: " + str(p_value))
            this_basin_key = int(basin_key)
            this_source_key = int(source)
            this_row = [this_basin_key,this_source_key,slope,std_err,r_value**2,p_value]
            OutDF.loc[counter] = this_row
            counter+=1

    return OutDF

def LinearRegressionSegmentedData(DataDirectory, DEM_prefix, basin_list=[]):
    """
    This function performs a linear regression on each of the segments for the
    SA data.

    Args:
        DataDirectory (str): the data directory
        DEM_prefix (str): the prefix of the DEM_prefix
        basin_list: a list of the basins to analyse, default = empty (all basins)

    Returns:
        pandas dataframe with the linear regression info

    Author: FJC
    """

    df = Helper.ReadSegmentedSAData(DataDirectory, DEM_prefix)

    # get a list of the basins if needed
    if basin_list == []:
        print ("You didn't give me a basin list so I will analyse all the basins")
        basin_list = df['basin_key'].unique()

    columns = ['basin_key', 'segment_number', 'regression_slope', 'std_err', 'R2', 'p_value']
    OutDF = pd.DataFrame(columns=columns)
    counter=0

    # get the segments for each basin
    for basin_key in basin_list:
        print ("THIS BASIN IS: "+str(basin_key))
        SegmentDF = df[df['basin_key'] == basin_key]

        # get the data for each individual segment number
        segments = SegmentDF['segment_number'].unique()

        for segment_no in segments:
            print ("Segment number is: "+str(segment_no))
            ThisDF = SegmentDF[SegmentDF['segment_number']==segment_no]
            #now regress the data for this segment to get the best fit m/n
            median_log_S = ThisDF['median_log_S']
            median_log_A = ThisDF['median_log_A']
            slope, intercept, r_value, p_value, std_err = stats.linregress(median_log_A,median_log_S)
            print("Slope: " +str(slope)+ " std_err: "+str(std_err)+ " R2 is: " + str(r_value**2) + " p value is: " + str(p_value) + " intercept is: " +str(intercept))

            this_key = int(basin_key)
            this_row = [this_key,int(segment_no),slope,std_err,r_value**2,p_value]
            OutDF.loc[counter] = this_row
            counter+=1

    return OutDF

def GetRangeMOverNSegmentedData(DataDirectory,fname_prefix,basin_list=[]):
    """
    This function performs a linear regression on each segment of the segmented
    SA data and finds the median and interquartile range from the segments for each basin.
    The number of removed segments is also reported (segments are removed if they
    predict a positive relationship between S and A).

    Args:
        DataDirectory (str): the data directory
        DEM_prefix (str): the prefix of the DEM_prefix
        basin_list: a list of the basins to analyse, default = empty (all basins)

    Returns:
        pandas dataframe with the basin SA info

    Author: FJC
    """

    # get the linear regression info for the segments
    SegmentDF = LinearRegressionSegmentedData(DataDirectory,fname_prefix,basin_list)

    # now for each basin, get the median, min and max m/n.
    # get a list of the basins if needed
    if basin_list == []:
        basin_list = SegmentDF['basin_key'].unique()

    #set up the output dataframe
    columns = ['basin_key', 'median_movern', 'min_movern', 'max_movern', 'FirstQ_movern', 'ThirdQ_movern', 'removed_segments']
    OutDF = pd.DataFrame(columns=columns)

    for basin_key in basin_list:
        # get the segment info for this basin
        ThisSegmentDF = SegmentDF[SegmentDF['basin_key'] == basin_key]
        moverns = ThisSegmentDF['regression_slope'][ThisSegmentDF['regression_slope'] < 0].tolist()
        removed_segments = ThisSegmentDF['segment_number'][ThisSegmentDF['regression_slope'] >= 0].tolist()

        if moverns !=[]:
            #find the median m/n
            median_movern = np.median(moverns)
            min_movern = np.min(moverns)
            max_movern = np.max(moverns)
            q75, q25 = np.percentile(moverns, [75 ,25])
            print ("The median is: " +str(median_movern)+", the max is: "+str(max_movern)+", the min is: "+str(min_movern), ", the q25 is: "+str(q25)+", the q75 is: " +str(q75))
            this_key = int(basin_key)
            this_row = [this_key,abs(median_movern),abs(min_movern),abs(max_movern),abs(q75),abs(q25),removed_segments]
            OutDF.loc[this_key] = this_row

    return OutDF

def GetRangeMOverNRawDataByChannel(DataDirectory,fname_prefix,basin_list=[]):
    """
    This function performs a linear regression on each tributary of the raw data
     and finds the median and interquartile range of the predicted m/n values.

    Args:
        DataDirectory (str): the data directory
        DEM_prefix (str): the prefix of the DEM_prefix
        basin_list: a list of the basins to analyse, default = empty (all basins)

    Returns:
        pandas dataframe with the basin SA info

    Author: FJC
    """

    # get the linear regression info for the segments
    RawDF = LinearRegressionRawDataByChannel(DataDirectory,fname_prefix,basin_list)
    print (RawDF)

    # now for each basin, get the median, min and max m/n.
    # get a list of the basins if needed
    if basin_list == []:
        basin_list = RawDF['basin_key'].unique()

    #set up the output dataframe
    columns = ['basin_key', 'median_movern', 'min_movern', 'max_movern', 'FirstQ_movern', 'ThirdQ_movern', 'removed_sources']
    OutDF = pd.DataFrame(columns=columns)

    for basin_key in basin_list:
        # get the segment info for this basin
        print ("Basin key is: "+str(basin_key))
        ThisDF = RawDF[RawDF['basin_key'] == basin_key]
        moverns = ThisDF['regression_slope'][ThisDF['regression_slope'] < 0].tolist()
        removed_sources = ThisDF['source_key'][ThisDF['regression_slope'] >= 0].tolist()

        if moverns != []:
            #find the median m/n
            median_movern = np.median(moverns)
            min_movern = np.min(moverns)
            max_movern = np.max(moverns)
            q75, q25 = np.percentile(moverns, [75 ,25])
            print ("The median is: " +str(median_movern)+", the max is: "+str(max_movern)+", the min is: "+str(min_movern), ", the q25 is: "+str(q25)+", the q75 is: " +str(q75))
            this_key = int(basin_key)
            this_row = [this_key,abs(median_movern),abs(min_movern),abs(max_movern),abs(q75),abs(q25),removed_sources]
            OutDF.loc[this_key] = this_row

    return OutDF

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Slope-area functions
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def SAPlotDriver(DataDirectory, DEM_prefix, FigFormat = 'show', size_format = "ESURF",
                 show_raw = True, show_segments = True,
                 cmap = plt.cm.Set1, n_colours = 10,
                 basin_keys = [], parallel=False):
    """
    This is a driver function that manages plotting of Slope-Area data

    Args:
        DataDirectory (str): the path to the directory with the csv file
        DEM_prefix (str): name of your DEM without extension
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        show_raw (bool): If true show raw data in background
        show_segments (bool): If true, show the segmented main stem,
        cmap (string or colourmap): the colourmap use to colour tributaries
        n_colours (int): The number of coulours used in plotting tributaries
        basin_keys (list): A list of the basin keys to plot. If empty, plot all the basins.
        parallel (bool): If true the data is in multiple files and must be merged

    Returns:
        Slope-area plot for each basin

    Author: SMM
    """
    from LSDPlottingTools import LSDMap_PointTools as PointTools

    print("These basin keys are: ")
    print(basin_keys)

    # read in binned data
    binned_csv_fname = DataDirectory+DEM_prefix+'_SAbinned.csv'
    print("I'm reading in the csv file "+binned_csv_fname)

    if not parallel:
        binnedPointData = PointTools.LSDMap_PointData(binned_csv_fname)
    else:
        binnedPointData = Helper.AppendSABinnedCSVs(DataDirectory, DEM_prefix)
        binnedPointData = PointTools.LSDMap_PointData(binned_csv_fname)

    # Read in the raw data
    if(show_raw):
        print("I am going to show the raw data.")
        all_csv_fname = DataDirectory+DEM_prefix+'_SAvertical.csv'
        if not parallel:
            allPointData = PointTools.LSDMap_PointData(all_csv_fname)
        else:
            allPointData = Helper.AppendSAVerticalCSVs(DataDirectory, DEM_prefix)
            allPointData = PointTools.LSDMap_PointData(all_csv_fname)

    # Read in the segmented data
    if(show_segments):
        print("I am going to show segments on the main stem.")
        segmented_csv_fname = DataDirectory+DEM_prefix+'_SAsegmented.csv'
        if not parallel:
            segmentedPointData = PointTools.LSDMap_PointData(segmented_csv_fname)
        else:
            segmentedPointData = Helper.AppendSASegmentedCSVs(DataDirectory, DEM_prefix)
            segmentedPointData = PointTools.LSDMap_PointData(segmented_csv_fname)

    # get the basin keys and check if the basins in the basin list exist
    basin = binnedPointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    Basin = np.asarray(basin)
    these_basin_keys = np.unique(Basin)

    print("The unique basin keys are: ")
    print(these_basin_keys)

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

    print("The final basin keys are:")
    print(final_basin_keys)


    this_cmap = cmap
    print("There are "+str(len(final_basin_keys))+"basins that I will plot")
    #basin_keys.append(0)
    # Loop through the basin keys, making a plot for each one
    for basin_key in final_basin_keys:
        print("I am making a plot for basin: "+str(basin_key))
        if(show_segments):
            if(show_raw):
                FileName = DEM_prefix+'_SA_plot_raw_and_segmented_basin%s.%s' %(str(basin_key),FigFormat)
                SegmentedWithRawSlopeAreaPlot(segmentedPointData, allPointData,
                                                              DataDirectory, FigFileName=FileName,
                                                              FigFormat=FigFormat, size_format=size_format,
                                                              basin_key=basin_key,cmap = this_cmap, n_colours = n_colours)
            else:
                FileName = DEM_prefix+'_SA_plot_segmented_basin%s.%s' %(str(basin_key),FigFormat)
                SegmentedSlopeAreaPlot(segmentedPointData, DataDirectory,
                                       FigFileName=FileName, FigFormat=FigFormat,
                                       size_format=size_format, basin_key=basin_key)

        else:
            if(show_raw):
                FileName = DEM_prefix+'_SA_plot_raw_and_binned_basin%s.%s' %(str(basin_key),FigFormat)
                BinnedWithRawSlopeAreaPlot(DataDirectory, DEM_prefix, FigFileName=FileName,
                                          FigFormat=FigFormat, size_format=size_format,
                                          basin_key=basin_key, n_colours = n_colours,
                                          cmap = this_cmap)
            else:
                print("You selected an option that doesn't produce any plots. Turn either show raw or show segments to True.")


def BinnedRegressionDriver(DataDirectory, DEM_prefix, basin_keys = []):
    """
    This function goes through a basin list and reports back the best fit
    m/n values for mainstem data, all data, and both of these with outliers removed

    Args:
        DataDirectory (str): the path to the directory with the csv file
        DEM_prefix (str): name of your DEM without extension
        basin_keys (list): A list of the basin keys to plot. If empty, plot all the basins.

    Author: SMM
    """
    from LSDPlottingTools import LSDMap_PointTools as PointTools

    print("These basin keys are: ")
    print(basin_keys)

    # read in binned data
    binned_csv_fname = DataDirectory+DEM_prefix+'_SAbinned.csv'
    print("I'm reading in the csv file "+binned_csv_fname)
    binnedPointData = PointTools.LSDMap_PointData(binned_csv_fname)

    # get the basin keys and check if the basins in the basin list exist
    basin = binnedPointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    Basin = np.asarray(basin)
    these_basin_keys = np.unique(Basin)

    #print("The unique basin keys are: ")
    #print(these_basin_keys)

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

    #print("The final basin keys are:")
    #print(final_basin_keys)

    print("There are "+str(len(final_basin_keys))+"basins that I will plot")
    mn_by_basin_dict = {}
    #basin_keys.append(0)
    # Loop through the basin keys, making a plot for each one
    for basin_key in final_basin_keys:

        (m1,m2,m3,m4) = BinnedRegression(binnedPointData, basin_key)
        this_basin_SA_mn = []
        this_basin_SA_mn.append(m1)
        this_basin_SA_mn.append(m2)
        this_basin_SA_mn.append(m3)
        this_basin_SA_mn.append(m4)

        mn_by_basin_dict[basin_key] = this_basin_SA_mn

    return mn_by_basin_dict




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

    # check if a directory exists for the SA plots. If not then make it.
    SA_directory = DataDirectory+'SA_plots/'
    if not os.path.isdir(SA_directory):
        os.makedirs(SA_directory)

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
        plt.savefig(SA_directory+FigFileName,format=save_fmt,dpi=500)
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

    # check if a directory exists for the SA plots. If not then make it.
    SA_directory = DataDirectory+'SA_plots/'
    if not os.path.isdir(SA_directory):
        os.makedirs(SA_directory)

    # Set up fonts for plots
    label_size = 10
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
    Raw_log_A =  RawPointData.QueryData('drainage_area')
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
        plt.savefig(SA_directory+FigFileName,format=save_fmt,dpi=500)
        fig.clf()

def BinnedWithRawSlopeAreaPlot(DataDirectory, fname_prefix, FigFileName = 'Image.pdf',
                       FigFormat = 'show',size_format = "ESURF",basin_key = '0',
                       cmap = plt.cm.Set1, n_colours = 5):
    """
    This function makes a slope-area plot from the chi mapping tool using the binned data.
    It plots all the sources and all the raw data as semitransparent background.

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the prefix of your DEM
        FigFileName (str): The name of the figure file
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        basin_key (int): the ID of the basin to make the plot for.
        colourmap (string or colormap object): The colourmap
        n_colour (int): The number of colours

    Returns:
         Does not return anything but makes a plot.

    Author: SMM (panda-fied by FJC)
    """
    import matplotlib.colors as colors
    import matplotlib.ticker

    # check if a directory exists for the SA plots. If not then make it.
    SA_directory = DataDirectory+'SA_plots/'
    if not os.path.isdir(SA_directory):
        os.makedirs(SA_directory)

    # Set up fonts for plots
    label_size = 10
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

    # Get the raw data and mask to the basin
    RawDF = Helper.ReadRawSAData(DataDirectory,fname_prefix)
    RawDF = RawDF[RawDF['basin_key']==basin_key]

    # get the binend data
    BinnedDF = Helper.ReadBinnedSAData(DataDirectory,fname_prefix)
    BinnedDF = BinnedDF[BinnedDF['basin_key'] == basin_key]

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
    sources = np.unique(RawDF['source_key'].as_matrix())
    print (sources)

    # get the regression stats for this basin
    RegressionDF = LinearRegressionRawData(DataDirectory, fname_prefix)
    this_movern = round(RegressionDF.get_value(basin_key,'regression_slope'),2)

    # Mask the data of the segments sequentially
    for idx,source in enumerate(sources):
        print ("Source key: "+str(source))
        # mask to just get the data for the basin of interest
        RawDFMask = RawDF[RawDF['source_key'] == source]
        BinnedDFMask = BinnedDF[BinnedDF['source_key'] == source]
        SourceMedianLogSlope = 10**(BinnedDFMask['median_log_S'].as_matrix())
        SourceMedianLogArea = 10**(BinnedDFMask['median_log_A'].as_matrix())
        print (SourceMedianLogArea)

        # Now add the colours for the segments
        source_colour = idx%n_colours
        tps_color = scalarMap.to_rgba(source_colour)
        ax.scatter(SourceMedianLogArea,SourceMedianLogSlope,c=tps_color,s=10,marker="o",lw=0.5,edgecolors='k',zorder=100)

        # get the errors
        first_quartile = 10**(BinnedDFMask['logS_FirstQuartile'].as_matrix())
        third_quartile = 10**(BinnedDFMask['logS_ThirdQuartile'].as_matrix())
        yerr_up = third_quartile-SourceMedianLogSlope
        yerr_down = SourceMedianLogSlope-first_quartile

        plt.errorbar(SourceMedianLogArea,SourceMedianLogSlope,yerr=[yerr_down,yerr_up],fmt='o',ms=1,ecolor=tps_color,zorder=0)

        # Plot the raw data
        SourceRawMedianS = RawDFMask['slope']
        SourceRawMedianA = RawDFMask['drainage_area']
        ax.scatter(SourceRawMedianA,SourceRawMedianS,c=tps_color,s=4,marker="+",lw=0.5,edgecolors='k',zorder=-10,alpha = 0.3)

    ax.set_xlabel('Drainage area (m$^2$)')
    ax.set_ylabel('Slope (m/m)')

    # log
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_title('Basin '+str(basin_key)+', best fit $m/n$='+str(this_movern), fontsize=10)

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
        plt.savefig(SA_directory+FigFileName,format=save_fmt,dpi=500)
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
    [MSresiduals,m_ms,b,r,pvalue,stderr]= LSDStats.linregress_residuals(MSAreaCompressed,MSSlopeCompressed)
    #print("slope of mainstem regression is: "+str(m))

    # see if there are any outlying residuals
    [new_x,new_y, is_outlier_vec, m_ms_remove_outlier,b]= LSDStats.remove_outlying_residuals(MSAreaCompressed,MSSlopeCompressed,MSresiduals)
    #print("Removed outlier slope from mainstem data: " +str(m))

    # get the regression from all the data
    [residuals,m_all,b,r,pvalue,stderr]= LSDStats.linregress_residuals(AreaCompressed,SlopeCompressed)
    #print("slope of all data is: "+str(m))

    # see if there are any outlying residuals
    [new_x,new_y, is_outlier_vec, m_all_remove_outlier,b]= LSDStats.remove_outlying_residuals(AreaCompressed,SlopeCompressed,residuals)
    #print("Removed outlier slope from all data: " +str(m))

    return(m_ms,m_ms_remove_outlier,m_all,m_all_remove_outlier)
