#=============================================================================
# These functions create figures for visualising the m/n selection using the chi
# mapping tool.
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
# set backend to run on server
import matplotlib
matplotlib.use('Agg')

import numpy as np
import LSDPlottingTools as LSDP
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.ticker as ticker
import pandas as pd
from matplotlib import colors
import math
import os
import subprocess
#from shapely.geometry import Polygon
from LSDMapFigure import PlottingHelpers as Helper
from LSDMapFigure.PlottingRaster import MapFigure
from LSDMapFigure.PlottingRaster import BaseRaster
from LSDPlottingTools import LSDMap_SAPlotting as SA
from LSDPlottingTools import joyplot



#=============================================================================
#=============================================================================
# ANALYSIS FUNCTIONS
# These functions analyse the data (for example, doing the outlier checking)
#=============================================================================
def SimpleMaxMLECheck(BasinDF):
    """
    This function checks through the basin dataframe and returns a dict with the basin key and the best fit
    m/n

    Args:
        BasinDF: pandas dataframe from the basin csv file.

    Returns:
        m_over_n_dict: dictionary with the best fit m/n for each basin, the key is the basin key
        and the value is the best fit m/n

    Author: FJC
    """
    # read in the basin csv
    basin_keys = list(BasinDF['basin_key'])
    #print BasinDF

    # remove the first 2 columns from DF
    del BasinDF['basin_key']
    del BasinDF['outlet_jn']

    # now find the index and value of the max MLE in each row
    MOverNs = list(BasinDF.idxmax(axis=1))
    MOverNs = [float(x.split()[-1]) for x in MOverNs]
    print ("MAX MOVERNS")
    print (MOverNs)

    # zip into a dictionary
    MOverNDict = dict(zip(basin_keys, MOverNs))
    return MOverNDict

def GetMOverNRangeMCPoints(BasinDF, start_movern=0.2, d_movern=0.1, n_movern=7):
    """
    This function checks through the MC points basin dataframe and returns the best fit and
    range of m/n values. The range is given as a list in the column (all values of m/n where
    the 3rd quartile is above the first quartile value for the best fit m/n)

    Args:
        BasinDF: pandas dataframe from the basin MC points csv file.
        start_movern (float): starting m/n value
        d_movern (float): step in m/n value
        n_movern (float): number of m/n values

    Returns:
        dataframe with the basin key, best fit m/n, and range of m/n for each basin

    Author: FJC
    """
    from scipy import stats

    # get the medians from the dataframe
    MedianDF = BasinDF.filter(regex='median')

    # find the median with the highest MLE for each basin
    MaxMedians = list(MedianDF.idxmax(axis=1))
    MaxMedians = [str(x) for x in MaxMedians]
    Median_MOverNs = []
    for median in MaxMedians:
        if "=" in median:
            Median_MOverNs.append(float(median.split("=")[-1]))
        else:
            Median_MOverNs.append(np.nan)
    #Median_MOverNs = [float(x.split("=")[-1]) for x in MaxMedians]

    # now find the first quartile that corresponds to this median
    FirstQDF = BasinDF.filter(regex='FQ')
    FirstQF_MLEs = []
    for i, median in enumerate(Median_MOverNs):
        if not np.isnan(median):
            # find the right column
            this_col = "FQ_MLE_m_over_n="+str(median)
            FirstQDF_mask = FirstQDF[this_col]
            FirstQF_MLEs.append(float(FirstQDF_mask.iloc[i]))
        else:
            FirstQF_MLEs.append(np.nan)

    # now, for each basin, find the columns in the 3rd quartile which are higher than the first Q MLE
    ThirdQDF = BasinDF.filter(regex='TQ')

    # change the column names to just have the m/n values
    column_names = list(ThirdQDF)
    column_names = [x.split("=")[-1] for x in column_names]
    ThirdQDF.columns = column_names

    # add the threshold first Q MLEs to the dataframe
    ThirdQDF['threshold'] = pd.Series(FirstQF_MLEs, index=ThirdQDF.index)
    print (ThirdQDF)
    # change DF to a boolean where values are greater than the threshold
    TempDF = ThirdQDF.drop('threshold', 1).gt(ThirdQDF['threshold'], 0)
    # get the column names where the values are greater than the threshold for each basin
    TempDF['Range_MOverNs'] = TempDF.apply(lambda x: ','.join(x.index[x]),axis=1)

    end_movern = start_movern+d_movern*(n_movern-1)

    # get the m/n values greater than the threshold to a list, then get the min and max
    Min_MOverNs = []
    Max_MOverNs = []
    Range_MOverNs = list(TempDF['Range_MOverNs'])
    for i in range (len(Range_MOverNs)):
        movern_str = Range_MOverNs[i].split(",")
        if '' in movern_str:
            Min_MOverNs.append(np.nan)
            Max_MOverNs.append(np.nan)
        else:
            movern_floats = [float(x) for x in movern_str]
            # we need to linearly interpolate between the nearest m/ns for the min
            # and the max.
            Min_MOverN = min(movern_floats)
            Max_MOverN = max(movern_floats)

            # ok, for the minimum, get the previous m/n value
            if Min_MOverN - d_movern < start_movern:
                new_min_movern = Min_MOverN
            else:
                movern_list = [Min_MOverN-d_movern, Min_MOverN]
                mle_list = [ThirdQDF[str(Min_MOverN-d_movern)][i],ThirdQDF[str(Min_MOverN)][i]]
                slope, intercept, r_value, p_value, std_err = stats.linregress(movern_list, mle_list)
                new_min_movern = (ThirdQDF['threshold'][i] - intercept)/slope

            # for the maximum get the next m/n value
            if Max_MOverN + d_movern > end_movern:
                new_max_movern = Max_MOverN
            else:
                movern_list = [Max_MOverN, Max_MOverN+d_movern]
                mle_list = [ThirdQDF[str(Max_MOverN)][i], ThirdQDF[(str(Max_MOverN+d_movern))][i]]
                slope, intercept, r_value, p_value, std_err = stats.linregress(movern_list, mle_list)
                new_max_movern = (ThirdQDF['threshold'][i] - intercept)/slope

            Min_MOverNs.append(new_min_movern)
            Max_MOverNs.append(new_max_movern)


    # write the output dataframe
    OutputDF = pd.DataFrame()
    OutputDF['basin_key'] = BasinDF['basin_key']
    OutputDF['Median_MOverNs'] = pd.Series(Median_MOverNs)
    OutputDF['FirstQ_threshold'] = pd.Series(FirstQF_MLEs)
    OutputDF['Min_MOverNs'] = pd.Series(Min_MOverNs)
    OutputDF['Max_MOverNs'] = pd.Series(Max_MOverNs)

    return OutputDF

def GetRangeMOverNChiResiduals(DataDirectory, fname_prefix, basin_list=[0], parallel=False):
    """
    This function reads in the CSV files with the chi residuals data and
    calculates the median best fit m/n along with the min and max from the
    first and third quartiles. These are returned as a pandas dataframe because
    I love pandas <3 <3

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_list: a list of the basins to make the plots for. If an empty list is passed then
        all the basins will be analysed. Default = basin 0.

    Returns:
        pandas dataframe with m/n data from the chi residuals analysis

    Author: FJC
    """
    # first let's read in the csv files with the residuals data
    if not parallel:
        dfs = Helper.ReadChiResidualsCSVs(DataDirectory,fname_prefix)
    else:
        dfs = Helper.AppendChiResidualsCSVs(DataDirectory,fname_prefix)

    # get the best fit m/n from the dataframes
    movern_data = []
    for i in range(len(dfs)):
        dfs[i] = dfs[i][dfs[i]['basin_key'].isin(basin_list)]
        ThisDF = dfs[i][dfs[i] < 0]

        indices = []
        for i, row in ThisDF.iterrows():
            # find the first index where it is not a nan
            index = row.first_valid_index()
            if index == None:
                index = row.index[-1]
            indices.append(float(index.split("=")[-1]))
        movern_data.append(indices)

    # the movern_data is a list of lists containing the information.
    # movern_data[0] = median
    # movern_data[1] = first quartile
    # movern_data[2] = third quartile
    # the second index is the basin key.
    # so to get the median of basin 3 you would use movern_data[0][3]

    # write the output dataframe
    OutputDF = pd.DataFrame()
    OutputDF['basin_key'] = dfs[0]['basin_key']
    OutputDF['Median_MOverNs'] = pd.Series(movern_data[0])
    OutputDF['FirstQ_MOverNs'] = pd.Series(movern_data[1])
    OutputDF['ThirdQ_MOverNs'] = pd.Series(movern_data[2])

    return OutputDF

def GetBestFitMOverNFromDisorder(BasinDF):
    """
    This function looks in the basin DF from the disorder CSV files
    and returns a dict with the basin ID and the best-fit m/n for each one.
    The best fit m/n is the one with the minimum disorder statistic.

    Args:
        BasinDF: basin disorder dataframe

    Returns: dict with the best fit m/n for each basin, key is the basin ID and value is the best fit m/n

    Author: FJC
    """
    # read in the basin csv
    basin_keys = list(BasinDF['basin_key'])
    #print BasinDF

    # remove the first 2 columns from DF
    del BasinDF['basin_key']
    del BasinDF['outlet_jn']

    # now find the index and value of the max MLE in each row
    MOverNs = list(BasinDF.idxmin(axis=1))
    MOverNs = [float(x.split()[-1]) for x in MOverNs]
    print ("Best-fit m/ns disorder")
    print (MOverNs)

    # zip into a dictionary
    MOverNDict = dict(zip(basin_keys, MOverNs))
    return MOverNDict


def CompareChiAndSAMOverN(DataDirectory, fname_prefix, basin_list=[0], start_movern=0.2, d_movern=0.1, n_movern=7):
    """
    This function compiles

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_list: a list of the basins to make the plots for. If an empty list is passed then
        all the basins will be analysed. Default = basin 0.
        start_movern (float): the starting m/n value. Default is 0.2
        d_movern (float): the increment between the m/n values. Default is 0.1
        n_movern (float): the number of m/n values analysed. Default is 7.

    Author: SMM
    """


    from LSDPlottingTools import LSDMap_PointTools as PointTools
    from LSDPlottingTools import LSDMap_SAPlotting as SAPlot

    # read in binned SA data
    binned_csv_fname = DataDirectory+fname_prefix+'_SAbinned.csv'
    print("I'm reading in the csv file "+binned_csv_fname)
    binnedPointData = PointTools.LSDMap_PointData(binned_csv_fname)

    # get the basin keys and check if the basins in the basin list exist
    basin = binnedPointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    Basin = np.asarray(basin)
    these_basin_keys = np.unique(Basin)

    print("The unique basin keys are: ")
    print(these_basin_keys)

    final_basin_keys = []
    # A bit of logic for checking keys
    if (len(basin_list) == 0):
        final_basin_keys = these_basin_keys
    else:
        for basin in basin_list:
            if basin not in these_basin_keys:
                print("You were looking for basin "+str(basin)+ " but it isn't in the basin keys.")
            else:
                final_basin_keys.append(basin)

    print("The final basin keys are:")
    print(final_basin_keys)

    # Now get all the m/n values from the basin list
    SA_movern_dict = SAPlot.BinnedRegressionDriver(DataDirectory, fname_prefix,
                                                   basin_keys = final_basin_keys)

    # Now do the same with the chi-derived m/n data
    (Outlier_counter, removed_sources_dict, best_fit_movern_dict, MLEs_dict) = CheckMLEOutliers(DataDirectory, fname_prefix, basin_list=final_basin_keys,
                                 start_movern=start_movern, d_movern=d_movern,
                                 n_movern=n_movern)

    # Now print to files
    print("The best fit m.n for chi analysis are: ")
    print(best_fit_movern_dict)

    # Now plot the figure
    PlotMOverNDicts(DataDirectory, fname_prefix, SA_movern_dict,best_fit_movern_dict, FigFormat = "png", size_format = "ESURF")

def CompareMOverNEstimatesAllMethods(DataDirectory, fname_prefix, basin_list=[0], start_movern=0.2, d_movern=0.1, n_movern=7, parallel=False, Chi_disorder=False):
    """
    This function reads in all the files with the data for the various methods of estimating
    the best fit m/n and produces a summary csv file which has the best fit m/n and uncertainty
    for each basin. This csv can then be used to plot the data.

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_list: a list of the basins to make the plots for. If an empty list is passed then
        all the basins will be analysed. Default = basin 0.
        start_movern (float): the starting m/n value. Default is 0.2
        d_movern (float): the increment between the m/n values. Default is 0.1
        n_movern (float): the number of m/n values analysed. Default is 7.
        Chi_disorder (bool): If true, will include the chi disorder stats

    Returns:
        writes a csv with the best fit m/n info for each basin

    Author: FJC
    """
    # check if a directory exists for the summary plots. If not then make it.
    summary_directory = DataDirectory+'summary_plots/'
    if not os.path.isdir(summary_directory):
        os.makedirs(summary_directory)

    # read in the full chi dataframe
    if not parallel:
        FullChiBasinDF = Helper.ReadBasinStatsCSV(DataDirectory,fname_prefix)
    else:
        FullChiBasinDF = Helper.AppendBasinCSVs(DataDirectory,fname_prefix)

    # Let's get the list of basins
    if basin_list == []:
        print("You didn't supply a list of basins, so I'm going to analyse all of them.")
        basin_list = list(FullChiBasinDF['basin_key'])

    # We're going to write our summary csv using pandas, woooo
    # First of all we need to set up the dataframe
    OutDF = pd.DataFrame()
    OutDF['basin_key'] = pd.Series(basin_list)

    # get the best fit m/n for each basin in the list from the full chi method
    FullChiMOverNDict = SimpleMaxMLECheck(FullChiBasinDF)
    OutDF['Chi_MLE_full'] = OutDF['basin_key'].map(FullChiMOverNDict)

    # get the best fit m/n from the points method
    if not parallel:
        PointsChiBasinDF = Helper.ReadMCPointsCSV(DataDirectory,fname_prefix)
    else:
        PointsChiBasinDF = Helper.AppendBasinPointCSVs(DataDirectory,fname_prefix)
    PointsChiBasinDF = PointsChiBasinDF[PointsChiBasinDF['basin_key'].isin(basin_list)]

    UncertaintyDF = GetMOverNRangeMCPoints(PointsChiBasinDF, start_movern, d_movern, n_movern)
    OutDF['Chi_MLE_points'] = UncertaintyDF['Median_MOverNs']
    OutDF['Chi_MLE_points_min'] = UncertaintyDF['Min_MOverNs']
    OutDF['Chi_MLE_points_max'] = UncertaintyDF['Max_MOverNs']

    print ("Now getting the m/n from the chi residuals")

    # get the best fit m/n from the chi residuals method
    ResidualsDF = GetRangeMOverNChiResiduals(DataDirectory, fname_prefix, basin_list, parallel=parallel)
    OutDF['Chi_residuals'] = ResidualsDF['Median_MOverNs']
    OutDF['Chi_residuals_min'] = ResidualsDF['FirstQ_MOverNs']
    OutDF['Chi_residuals_max'] = ResidualsDF['ThirdQ_MOverNs']

    print ("Getting the m/n from the SA data")

    # get the best fit m/n from the raw SA data
    RawSADF = SA.LinearRegressionRawData(DataDirectory,fname_prefix,basin_list)
    OutDF['SA_raw'] = RawSADF['regression_slope']
    OutDF['SA_raw_sterr'] = RawSADF['std_err']
    OutDF['SA_raw_R2'] = RawSADF['R2']
    OutDF['SA_raw_p'] = RawSADF['p_value']

    # get the SA tributary information
    SATribsDF = SA.GetRangeMOverNRawDataByChannel(DataDirectory,fname_prefix,basin_list)
    OutDF['SA_tribs'] = SATribsDF['median_movern']
    OutDF['SA_tribs_min'] = SATribsDF['FirstQ_movern']
    OutDF['SA_tribs_max'] = SATribsDF['ThirdQ_movern']

    # get the best fit m/n from the segmented SA data
    SASegmentedDF = SA.GetRangeMOverNSegmentedData(DataDirectory,fname_prefix,basin_list)
    OutDF['SA_segments'] = SASegmentedDF['median_movern']
    OutDF['SA_segments_min'] = SASegmentedDF['FirstQ_movern']
    OutDF['SA_segments_max'] = SASegmentedDF['ThirdQ_movern']

    if Chi_disorder:
        DisorderDF = Helper.ReadDisorderUncertCSV(DataDirectory,fname_prefix)
        OutDF['Chi_disorder'] = DisorderDF['median']
        OutDF['Chi_disorder_min'] = DisorderDF['first_quartile']
        OutDF['Chi_disorder_max'] = DisorderDF['third_quartile']

    # print the SA segment data
    OutSAname = "_SA_segment_summary.csv"
    out_sa_name = summary_directory+fname_prefix+OutSAname
    SASegmentedDF.to_csv(out_sa_name,index=False)

    # now write the output dataframe to a csv file
    OutCSVname = "_movern_summary.csv"
    outname = summary_directory+fname_prefix+OutCSVname
    OutDF.to_csv(outname,index=False)


def CheckMLEOutliers(DataDirectory, fname_prefix, basin_list=[0], start_movern=0.2, d_movern=0.1, n_movern=7, parallel=False):
    """
    This function uses the fullstats files to search for outliers in the
    channels. It loops through m/n values and for each m/n value calculates which
    tributaries have outlying MLE and RMSE values. Each time a tributary has an outlying value,
    an outlier counter is incremented so that for each basin we end up with the number
    of times each of its tributaries has been selected as an outlier. This outlier counter is
    returned. From this outlier counter, the tributaries with outliers are ranked in
    descending order of outlier counts. These are iteratively removed. So for example,
    if one tributary has 4 outlier counts and another has 2, the first trib is removed,
    and MLE is recalculated, and then the second is also removed and MLE is recalculated again.
    The removed tributaries and the order in which they are removed is contained within
    the dictionary removed_sources_dict. Finally, a dictionary with the m/n value
    for each basin with the lowest MLE is returned. Each basin has a list where the first
    element is the MLE with no tributaries removed, the second with the first outlier
    count removed, etc. Note that if two tributaries ahve the same outlier counts
    then these are removed at the same time. The key into the dictionaries are the basin numbers.

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_list: a list of the basins to make the plots for. If an empty list is passed then
        all the basins will be analysed. Default = basin 0.
        start_movern (float): the starting m/n value. Default is 0.2
        d_movern (float): the increment between the m/n values. Default is 0.1
        n_movern (float): the number of m/n values analysed. Default is 7.

    Returns:
        Outlier_counter (dict): This is a dictionary where the key is the basin
        number and the values are lists that have the cumulative number of times
        each tributary in the basin is detected to be an outlier.
        removed_sources_dict (dict): This is a dictionary where the key is the basin
        number and the values are lists of lists where the top level list is the steps of removing channels
        and the second level are the channels themselves. The reason there are nested
        lists and not individual channels being removed is because sometimes two
        channels have the same outlier counts and so they are removed at the same time.
        best_fit_movern_dict (dict): This is a dictionary where the key is the basin
        number and the values lists of the m/n ratio that has the lowest MLE. These are lists beacuse
        each element in the list represents the number of outlying tributaries removed.
        The first element is where no tributaries are removed.
        MLEs_dict (dictionary of arrays): The key is the basin number and the value is an
        array containing all the MLE values for each m/n and
        each iterated removed tributary.
        To get the MLEs of, say, the first removal for the 3rd basin you would use MLEs_dict[2][:,1]

    Author: SMM
    """

    # Get a vector of the m over n values
    end_movern = start_movern+d_movern*(n_movern-1)
    m_over_n_values = np.linspace(start_movern,end_movern,n_movern)
    movern_strs = []
    for movern in m_over_n_values:
        movern_str = '%.2f'%movern
        if movern_str.endswith('0'):
            movern_str = movern_str[:-1]
            movern_strs.append(movern_str)

    # we open the first file just so that we can get a counter list
    print ("PARALLEL = ", parallel)
    if not parallel:
        FirstDF = Helper.ReadFullStatsCSV(DataDirectory,fname_prefix,movern_strs[0])
    else:
        FirstDF = Helper.AppendFullStatsCSVs(DataDirectory,movern_strs[0],fname_prefix)

    # get the number of basins
    basin_keys = list(FirstDF['basin_key'])
    basin_keys = [float(x) for x in basin_keys]

    # get the list of basins
    if basin_list == []:
        print("You didn't give me a list of basins, so I'll just run the analysis on all of them!")
        basin_list = basin_keys
        basin_set = set(basin_list)
        basin_list = list(basin_set)
        basin_list = [int(i) for i in basin_list]

    # make a data object that will hold the counters
    Outlier_counter = {}
    # loop through the basins
    for basin in basin_list:
        # mask the data so you only get the correct basin
        FirstDF_basin = FirstDF[FirstDF['basin_key'] == basin]

        trib_values = list(FirstDF_basin['test_source_key'])
        n_nodes = len(trib_values)
        # make the counter with zeros
        this_counter = np.zeros(n_nodes)
        Outlier_counter[basin] = this_counter

    # Now we loop through all the files, calculating the outliers
    for m_over_n in movern_strs:

        print("This_m_over_n is: "+m_over_n)

        #load the file
        if not parallel:
            FullStatsDF = Helper.ReadFullStatsCSV(DataDirectory,fname_prefix,m_over_n)
        else:
            FullStatsDF = Helper.AppendFullStatsCSVs(DataDirectory,m_over_n,fname_prefix)


        # loop through the basins
        for basin in basin_list:

            # mask the data so you only get the correct basin
            FullStatsDF_basin = FullStatsDF[FullStatsDF['basin_key'] == basin]

            # extract the relevant data
            MLE_values = list(FullStatsDF_basin['MLE'])
            RMSE_values = list(FullStatsDF_basin['RMSE'])
            trib_values = list(FullStatsDF_basin['test_source_key'])
            print ('N values: '+str(len(MLE_values))+' '+str(len(RMSE_values))+' '+str(len(trib_values)))

            # now get the outliers
            MLE_array = np.asarray(MLE_values)
            RMSE_array = np.asarray(RMSE_values)

            # Get the outliers using the MAD-based outlier function
            RMSE_outliers = LSDP.lsdstatsutilities.is_outlier(RMSE_array)
            MLE_outliers = LSDP.lsdstatsutilities.is_outlier(RMSE_array)

            # now check each of the outlier arrays to see if e need to flip the array
            RMSE_index_min = np.argmin(RMSE_array)
            #RMSE_min = np.min(RMSE_array)

            # if the max MLE is an outlier, flip the outlier vector
            if (RMSE_outliers[RMSE_index_min]):
                RMSE_outliers = [not i for i in RMSE_outliers]

            MLE_index_max = np.argmax(MLE_array)

            # if the max MLE is an outlier, flip the outlier vector
            if (MLE_outliers[MLE_index_max]):
                MLE_outliers = [not i for i in MLE_outliers]

            # turn the outliers vector into an integer
            int_Outlier = [int(i) for i in MLE_outliers]
            #print("Integer outliers are: ")
            #print(int_Outlier)

            # add this outlier counter to the outlier dict
            Outlier_counter[basin] = Outlier_counter[basin]+int_Outlier

    # Now try to calculate MLE by removing outliers

    # Set up dicts where the keys are the basin numbers and the
    # elements tell the m/n values of the best fit MLE, with each element
    # representing incrementally removed tributaries
    # The removed_sources_dict refers to the sources being removed each step
    best_fit_movern_dict = {}
    removed_sources_dict = {}
    MLEs_dict = {}
    for basin_number in basin_list:
        remove_list_index,movern_of_max_MLE,MLEs = Iteratively_recalculate_MLE_removing_outliers_for_basin(Outlier_counter,
                                                                                DataDirectory,
                                                                                fname_prefix,
                                                                                basin_number,
                                                                                start_movern,
                                                                                d_movern,
                                                                                n_movern, parallel)
        best_fit_movern_dict[basin_number] = movern_of_max_MLE
        removed_sources_dict[basin_number] = remove_list_index
        MLEs_dict[basin_number] = MLEs

    print("Here are the vitalstatisix, chief: ")
    print(best_fit_movern_dict)

    return Outlier_counter, removed_sources_dict, best_fit_movern_dict, MLEs_dict

def Iteratively_recalculate_MLE_removing_outliers_for_basin(Outlier_counter, DataDirectory, fname_prefix, basin_number, start_movern=0.2, d_movern=0.1, n_movern=7, parallel=False):
    """
    This function drives the calculations for removing outliers incrementally
    from the MLE calculations. This is specific to a basin.
    It calls functions for specific m/n values

    Args:
        Outlier_counter (dict): The dictionary containing the outlier lists for each basin
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_number (list): a list of the basins to make the plots for. If an empty list is passed then
        all the basins will be analysed. Default = basin 0.
        start_movern (float): the starting m/n value. Default is 0.2
        d_movern (float): the increment between the m/n values. Default is 0.1
        n_movern (float): the number of m/n values analysed. Default is 7.

    Returns:
        remove_list_index (list of list): This is the sequence of tributaries that will be removed
        movern_of_max_MLE (list): The m over n values of the maximum MLE sequentially for removed tributaries
        MLEs (array): This is an array containing all the MLE values for each m/n and
        each iterated removed tributary. To get the MLEs of, say, the first removal you would use MLEs[:,1]

    Author: SMM
    """

    # get the outlier counter for this basin
    thisBasinOutlierCounter = Outlier_counter[basin_number]

    print("I am going to recalculate MLE for basin: "+str(basin_number))
    #print("The counter list is:")
    #print(thisBasinOutlierCounter)

    # Get the sroted version and the indices into the sorted version
    sort_index = np.argsort(thisBasinOutlierCounter)
    sorted_outliers = np.sort(thisBasinOutlierCounter)

    # make sure the sourted outliers are ints
    int_sorted_outliers = [int(i) for i in sorted_outliers]

    # we need to reverse these lists so that the biggest outlier counts come first
    sort_index = sort_index[::-1]
    sorted_outliers = sorted_outliers[::-1]

    # get all the duplicates, with the total number of duplicates for each counter
    # This uses the unbelievably handy collections.Counter tool
    from collections import Counter
    all_counter_dict=Counter(int_sorted_outliers)

    # pop out the zero duplicates: we don't exclude non-outlier data
    all_counter_dict.pop(0, None)

    # now we need to iteratively remove the offending counters.
    remove_list = []
    remove_list_index = []
    last_count = -1

    # Only enter the loop if there are positive counts
    if all_counter_dict != 0:

        # Now loop through the sorted outlier count
        for idx,sorted_outlier_count in enumerate(sorted_outliers):

            # first check if this has more than one count
            if sorted_outlier_count == 0:
                all_removed = True
            else:
                # if this is the first element, add it to the remove list
                #if len(remove_list) == 0:
                #    remove_list.append([])
                #    remove_list_index.append([])

                # get the number of counts
                #n_counts = all_counter_dict[sorted_outlier_count]

                # either append the count and index to the current count
                # or make a new list for the next count
                if sorted_outlier_count != last_count:
                    remove_list.append([])
                    remove_list_index.append([])

                remove_list[-1].append(int(sorted_outlier_count))
                remove_list_index[-1].append(sort_index[idx])

                last_count = sorted_outlier_count

    # now print out the lists
    #print("remove list is: ")
    #print(remove_list)

    #print("and index of this list is: ")
    #print(remove_list_index)

    # now we loop through m over n, incrementally removing the tributaries
    # this is done by copying the MLE vector and then replacing
    # the offending channels with an MLE of 1
    # Get a vector of the m over n values
    end_movern = start_movern+d_movern*(n_movern-1)
    m_over_n_values = np.linspace(start_movern,end_movern,n_movern)

    # Now get the movern values
    movern_of_max_MLE, MLEs = Calculate_movern_after_iteratively_removing_outliers(m_over_n_values,
                                                                             DataDirectory,
                                                                             fname_prefix,
                                                                             basin_number,
                                                                             remove_list_index, parallel)

    # Returns the remove_list_index, which is a list where each element
    # is a list of tributaries removed in an iteration,
    # And the maximum MLE values in each iteration.
    return remove_list_index,movern_of_max_MLE, MLEs



def Calculate_movern_after_iteratively_removing_outliers(movern_list, DataDirectory,
                                                         fname_prefix, basin_number,
                                                         remove_list_index, parallel=False):
    """
    This function takes the remove list index, which contains information about
    the sequence of tributaries to be removed, and then recalculates MLE by incrementally
    removing tributaries.

    Args:
        movern_list (float): m/n value.
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_number (int): The basin you want
        remove_list_index (list of lists): This contains information about what tributaries to remove

    Returns:
        movern_of_max_MLE (list): A list containing the m/n values of the basin after outlying
        tributaries have been removed. Each element in the list represents an increment of
        tributary removed. The first element is with no tributaries removed.
        MLEs (array): This is an array containing all the MLE values for each m/n and
        each iterated removed tributary. To get the MLEs of, say, the first removal you would use MLEs[:,1]

    Author: SMM
    """
    # Loop through m over n values and recalculate MLE values after removing
    # the outlying data
    All_MLE = []
    for m_over_n in movern_list:
        movern_str = '%.2f'%m_over_n     # need to convert to a string here so files are read correctly. FJC 16/02/18
        if movern_str.endswith('0'):
            movern_str = movern_str[:-1]
    #    print movern_str

        MLE_vals = RecalculateTotalMLEWithRemoveList(DataDirectory, fname_prefix,
                                                     movern_str,basin_number, remove_list_index, parallel)
        All_MLE.append(MLE_vals)


    MLEs = np.asarray(All_MLE)
    index_of_maximums = np.argmax(MLEs,0)

    movern_of_max_MLE = []
    for index in index_of_maximums:
        movern_of_max_MLE.append(movern_list[index])

    # This is required because linspace gives floating point errors
    movern_of_max_MLE = np.around(movern_of_max_MLE,4)

    print("The MLEs for no removal are: ")
    print(MLEs[:,0])

    # Return the m/n ratio with the biggest MLE, but also the array
    # with all the MLE values
    return movern_of_max_MLE, MLEs

def RecalculateTotalMLEWithRemoveList(DataDirectory, fname_prefix,
                                      movern,basin_number, remove_list_index, parallel=False):
    """
    This function takes the remove list index and then recalculates MLE by incrementally
    removing tributaries

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        movern (float): m/n value.
        basin_number (int): The basin you want
        remove_list_index (list of lists): This contains information about what tributaries to remove

    Returns:
        MLE_vals (list): The MLE data with incrementally removed tributaries

    Author: SMM
    """
    #load the file
    if not parallel:
        FullStatsDF = Helper.ReadFullStatsCSV(DataDirectory,fname_prefix,movern)
    else:
        FullStatsDF = Helper.AppendFullStatsCSVs(DataDirectory,movern,fname_prefix)

    # mask the data so you only get the correct basin
    FullStatsDF_basin = FullStatsDF[FullStatsDF['basin_key'] == basin_number]

    # extract the relevant data
    MLE_values = list(FullStatsDF_basin['MLE'])
    #trib_values = list(FullStatsDF_basin['test_source_key'])

    # get the MLE_values as an array
    MLE_values = np.asarray(MLE_values)

    #print("The total MLE is: ")
    #print(np.prod(MLE_values))

    MLE_vals = []
    MLE_vals.append(np.prod(MLE_values))

    # now loop through the remove list
    remove_list = []
    for stuff_to_remove in remove_list_index:
        this_MLE = MLE_values

        # need to extend the remove list with these tribs
        remove_list.extend(stuff_to_remove)

        # now set the MLE of the removed tribs to 1
        for idx in remove_list:
            this_MLE[idx] = 1

        #print("\n REMOVING; The total MLE is: ")
        #print(np.prod(this_MLE))
        MLE_vals.append(np.prod(this_MLE))

    return MLE_vals

#=============================================================================
# PLOTTING FUNCTIONS
# Make plots of the m/n analysis
#=============================================================================

def MakePlotsWithMLEStats(DataDirectory, fname_prefix, basin_list = [0],
                  start_movern = 0.2, d_movern = 0.1, n_movern = 7, parallel=False):
    """
    This function makes a chi-elevation plot for each basin and each value of
    m/n and prints the MLE value between the tributaries and the main stem.
    The plot with the maximum value of MLE is highlighted in red (suggesting
    that this should be the appropriate m/n value for this basin). Channels
    are coloured by elevation.

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_list: a list of the basins to make the plots for. If an empty list is passed then
        all the basins will be analysed. Default = basin 0.
        start_movern (float): the starting m/n value. Default is 0.2
        d_movern (float): the increment between the m/n values. Default is 0.1
        n_movern (float): the number of m/n values analysed. Default is 7.

    Returns:
        Plot of each m/n value for each basin.

    Author: SMM, modified by FJC
    """

    # check if a directory exists for the chi plots. If not then make it.
    MLE_directory = DataDirectory+'basic_chi_plots/'
    if not os.path.isdir(MLE_directory):
        os.makedirs(MLE_directory)

    profile_suffix = "_movern.csv"
    basin_stats_suffix = "_movernstats_basinstats.csv"

    movern_profile_file = fname_prefix+profile_suffix
    movern_basin_stats_file = fname_prefix+basin_stats_suffix

    end_movern = start_movern+d_movern*(n_movern-1)
    m_over_n_values = np.linspace(start_movern,end_movern,n_movern)

    # get the maximum MLE of each basin
    if not parallel:
        pd_DF = Helper.ReadBasinStatsCSV(DataDirectory, fname_prefix)
    else:
        pd_DF = Helper.AppendBasinStatsCSVs(DataDirectory,fname_prefix)

    shp = pd_DF.shape
    max_MLEs = []
    max_MLEs_index = []
    for i in range(0,shp[0]):
        #print("I is: "+str(i))
        a = pd_DF.loc[[i]]
        b = np.asarray(a)
        c = b[0,2:]
        max_MLEs.append(max(c))
        max_MLEs_index.append(np.argmax(c))
    m_over_n_of_max = []

    for idx in max_MLEs_index:
        m_over_n_of_max.append(m_over_n_values[idx])

    print("The m over n of these max are: ")
    print(m_over_n_of_max)

    n_basins = len(max_MLEs)

    label_size = 12

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size
    size_format = "default"
    # make a figure,
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        #l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        #l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
        #l_pad = -35

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[5:100,10:95])

    # load the m_over_n data file
    thisPointData = LSDP.LSDMap_PointData(DataDirectory,movern_profile_file,parallel=parallel)
    allBasinStatsData = LSDP.LSDMap_PointData(DataDirectory+movern_basin_stats_file)

    print("m over n values are: ")
    print(m_over_n_values)

    mn_legends = []
    for mn in m_over_n_values:
        mn_legends.append("m_over_n = "+str(mn))

    print("The mn labels are: ")
    print(mn_legends)

    # get the data form the profiles
    elevation = thisPointData.QueryData('elevation')
    elevation = [float(x) for x in elevation]
    basin = thisPointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    source = thisPointData.QueryData('source_key')
    source = [int(x) for x in source]

    # get the basin keys
    allstats_basinkeys = allBasinStatsData.QueryData("basin_key")
    allstats_basinkeys = [int(x) for x in allstats_basinkeys]

    # need to convert everything into arrays so we can mask different basins
    Elevation = np.asarray(elevation)
    Basin = np.asarray(basin)
    #Source = np.asarray(source)

    # Loop through m/n values aggregating data
    for idx,mn in enumerate(m_over_n_values):

        counter = str(idx).zfill(3)
        print("Counter is: "+counter)
        # first get the chi values for this m_over_n
        #print ("mn is: " + str(mn))
        #print("index is: "+str(idx))
        mn_legend = "m_over_n = "+str(mn)
        #print("I am looking for the data element: "+mn_legend)
        this_chi = thisPointData.QueryData(mn_legend)

        # get the MLE value for this m/n
        this_MLE = allBasinStatsData.QueryData(mn_legend)

        # convert to a numpy array for masking
        Chi = np.asarray(this_chi)

        # some info about the chi and elevation values
        #max_chi = np.amax(Chi)
        #max_Elevation = np.amax(Elevation)
        #min_Elevation = np.amin(Elevation)

        #z_axis_min = int(min_Elevation/10)*10
        #z_axis_max = int(max_Elevation/10)*10+10
        #chi_axis_max = int(max_chi/5)*5+5

        # Now mask the data. Initially we will do only basin 0
        if basin_list == []:
            print("You didn't give me any basins so I assume you want all of them.")
            basin_list = range(0,n_basins-1)

        for basin_key in basin_list:
            # now we need to find out if this basin is in the allstats file,
            # and if so what index it is
            #this_basin_index = -99
            #for ii,bk in enumerate(allstats_basinkeys):
            #    #print("index: "+str(ii)+" and basin_key is: "+str(bk))
            #    if (bk == basin_key):
            #        this_basin_index = idx

            #if(this_basin_index != -99):
            #    MLE = this_MLE[this_basin_index]
            #else:
            #    MLE = "NaN"
            MLE = this_MLE[basin_key]
            #MLE_str = str(MLE)
            #short_MLE = str("%03.02e" % round(MLE,2))
            short_MLE = str(round(MLE,3))
            print("The short MLE is: "+short_MLE)

            #print("The MLE of this basin for this m over n is: "+short_MLE)

            # this gets the mask (for the chosen basin)
            m = np.ma.masked_where(Basin!=basin_key, Basin)

            # this is the masked chi value
            maskX = np.ma.masked_where(np.ma.getmask(m), Chi)
            maskElevation = np.ma.masked_where(np.ma.getmask(m), Elevation)

            # now plot the data with a colourmap
            ax.scatter(maskX,maskElevation,s=2.5, c=maskElevation,cmap="terrain",edgecolors='none')

            # some formatting of the figure
            ax.spines['top'].set_linewidth(1)
            ax.spines['left'].set_linewidth(1)
            ax.spines['right'].set_linewidth(1)
            ax.spines['bottom'].set_linewidth(1)

            # make the lables
            ax.set_xlabel("$\chi$ (m)")
            ax.set_ylabel("Elevation (m)")

            # This affects all axes because we set share_all = True.
            #ax.set_ylim(z_axis_min,z_axis_max)
            #ax.set_ylim(z_axis_min,z_axis_max)
            #ax.set_xlim(0,chi_axis_max)
            #ax.text("Basin = " +mn_legend+", $m/n$ = "+str(mn))

            #newline = "\n"
            title_string = "Basin "+str(basin_key)+", $m/n$ = "+str(mn)
            title_string2 = "MLE = "+short_MLE
            ax.text(0.05, 0.95, title_string,
                    verticalalignment='top', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='black', fontsize=10)
            # print("The basin index is: "+str(basin_key)+" and the max index is: "+str(max_MLEs_index[basin_key]))
            # if( idx == max_MLEs_index[basin_key]):
            #     print("This m/n is: "+str(mn)+" and it is the maximum MLE")
            #     ax.text(0.05, 0.88, title_string2+", maximum MLE in basin.",
            #         verticalalignment='top', horizontalalignment='left',
            #         transform=ax.transAxes,
            #         color='red', fontsize=10)
            # else:
            #     ax.text(0.05, 0.88, title_string2,
            #         verticalalignment='top', horizontalalignment='left',
            #         transform=ax.transAxes,
            #         color='black', fontsize=10)

            #save the plot
            newFilename = MLE_directory+"Chi_profiles_basin_"+str(basin_key)+"_"+counter+".png"

            # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
            ax.tick_params(axis='both', width=1, pad = 2)
            for tick in ax.xaxis.get_major_ticks():
                tick.set_pad(2)

            FigFormat = "png"
            plt.savefig(newFilename,format=FigFormat,dpi=300)
            ax.cla()

def MakeChiPlotsMLE(DataDirectory, fname_prefix, basin_list=[0], start_movern=0.2, d_movern=0.1, n_movern=7,
                    size_format='ESURF', FigFormat='png', animate=False, keep_pngs=False, parallel=False):
    """
    This function makes chi-elevation plots for each basin and each value of m/n
    where the channels are coloured by the MLE value compared to the main stem.
    The main stem is plotted in black.

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_list: a list of the basins to make the plots for. If an empty list is passed then
        all the basins will be analysed. Default = basin 0.
        start_movern (float): the starting m/n value. Default is 0.2
        d_movern (float): the increment between the m/n values. Default is 0.1
        n_movern (float): the number of m/n values analysed. Default is 7.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        animate (bool): If this is true then it creates a movie of the chi-elevation plots coloured by MLE.
        keep_pngs (bool): If this is false and the animation flag is true, then the pngs are deleted and just the video is kept.

    Returns:
        Plot of each m/n value for each basin.

    Author: FJC
    """
    from matplotlib.ticker import FormatStrFormatter

    # check if a directory exists for the chi plots. If not then make it.
    MLE_directory = DataDirectory+'chi_plots/'
    if not os.path.isdir(MLE_directory):
        os.makedirs(MLE_directory)

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

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[10:95,5:80])
    #colorbar axis
    ax2 = fig.add_subplot(gs[10:95,82:85])

    # read in the csv files
    if not parallel:
        ProfileDF = Helper.ReadChiProfileCSV(DataDirectory, fname_prefix)
        BasinStatsDF = Helper.ReadBasinStatsCSV(DataDirectory, fname_prefix)
    else:
        ProfileDF = Helper.AppendMovernCSV(DataDirectory, fname_prefix)
        BasinStatsDF = Helper.AppendBasinCSVs(DataDirectory, fname_prefix)

    # get the number of basins
    basin_keys = list(BasinStatsDF['basin_key'])
    basin_keys = [int(x) for x in basin_keys]

    # get the list of basins
    if basin_list == []:
        print("You didn't give me a list of basins, so I'll just run the analysis on all of them!")
        basin_list = basin_keys

    # loop through each m over n value
    end_movern = start_movern+d_movern*(n_movern-1)
    m_over_n_values = np.linspace(start_movern,end_movern,n_movern)

    # best fit moverns
    best_fit_moverns = SimpleMaxMLECheck(BasinStatsDF)

    for m_over_n in m_over_n_values:
        # read in the full stats file

        #Stupid floating point representation issues
        movern_str = "%.2f" % round(m_over_n,2)
        if movern_str.endswith('0'):
            movern_str = movern_str[:-1]

        print("This m/n is: "+movern_str)
        if not parallel:
            FullStatsDF = Helper.ReadFullStatsCSV(DataDirectory,fname_prefix,movern_str)
        else:
            FullStatsDF = Helper.AppendFullStatsCSVs(DataDirectory,movern_str,fname_prefix)

        # loop through all the basins in the basin list
        for basin_key in basin_list:
            print("This basin key is %s") %str(basin_key)

            # mask the data frames for this basin
            ProfileDF_basin = ProfileDF[ProfileDF['basin_key'] == basin_key]
            FullStatsDF_basin = FullStatsDF[FullStatsDF['basin_key'] == basin_key]

            #print FullStatsDF_basin
            print ("Getting the reference_source_key")

            print(FullStatsDF_basin.iloc[0]['reference_source_key'])

            # get the data frame for the main stem
            ProfileDF_MS = ProfileDF_basin[ProfileDF_basin['source_key'] == FullStatsDF_basin.iloc[0]['reference_source_key']]

            # get the data frame for the tributaries
            ProfileDF_basin = ProfileDF_basin[ProfileDF_basin['source_key'] != FullStatsDF_basin.iloc[0]['reference_source_key']]
            # merge with the full data to get the MLE for the tributaries
            ProfileDF_tribs = ProfileDF_basin.merge(FullStatsDF_basin, left_on = "source_key", right_on = "test_source_key")

            # get the chi and elevation data for the main stem
            movern_key = 'm_over_n = %s' % movern_str
            MainStemX = list(ProfileDF_MS[movern_key])
            MainStemElevation = list(ProfileDF_MS['elevation'])

            # get the chi, elevation, and MLE for the tributaries
            TributariesX = list(ProfileDF_tribs[movern_key])
            TributariesElevation = list(ProfileDF_tribs['elevation'])
            TributariesMLE = list(ProfileDF_tribs['MLE'])

            # get the colourmap to colour channels by the MLE value
            #NUM_COLORS = len(MLE)
            MLE_array = np.asarray(TributariesMLE)
            this_cmap = plt.cm.coolwarm
            cNorm  = colors.Normalize(vmin=np.min(MLE_array), vmax=np.max(MLE_array))
            plt.cm.ScalarMappable(norm=cNorm, cmap=this_cmap)

            # now plot the data with a colourmap
            sc = ax.scatter(TributariesX,TributariesElevation,c=TributariesMLE,cmap=this_cmap, norm=cNorm, s=2.5, edgecolors='none')
            ax.plot(MainStemX,MainStemElevation,lw=2, c='k')

            # some formatting of the figure
            ax.spines['top'].set_linewidth(1)
            ax.spines['left'].set_linewidth(1)
            ax.spines['right'].set_linewidth(1)
            ax.spines['bottom'].set_linewidth(1)

            # make the lables
            ax.set_xlabel("$\chi$ (m)")
            ax.set_ylabel("Elevation (m)")

            # the best fit m/n
            best_fit_movern = best_fit_moverns[basin_key]

            # label with the basin and m/n
            title_string = "Basin "+str(basin_key)+", "+ r"$\theta$ = "+movern_str
            if best_fit_movern == m_over_n:
                ax.text(0.05, 0.95, title_string,
                        verticalalignment='top', horizontalalignment='left',
                        transform=ax.transAxes,
                        color='red', fontsize=10)
            else:
                ax.text(0.05, 0.95, title_string,
                        verticalalignment='top', horizontalalignment='left',
                        transform=ax.transAxes,
                        color='black', fontsize=10)

            # add the colorbar
            colorbarlabel = "$MLE$"
            cbar = plt.colorbar(sc,cmap=this_cmap,spacing='uniform', orientation='vertical',cax=ax2)
            cbar.set_label(colorbarlabel, fontsize=10)
            ax2.set_ylabel(colorbarlabel, fontname='Arial', fontsize=10)
            ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

            #save the plot
            newFilename = MLE_directory+"MLE_profiles"+str(basin_key)+"_"+movern_str+"."+str(FigFormat)

            # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
            ax.tick_params(axis='both', width=1, pad = 2)
            for tick in ax.xaxis.get_major_ticks():
                tick.set_pad(2)

            plt.savefig(newFilename,format=FigFormat,dpi=300)
            ax.cla()
            ax2.cla()

    if animate:
        # animate the pngs using ffmpeg
        system_call = "ffmpeg -framerate 3 -pattern_type glob -i '"+MLE_directory+"MLE_profiles*.png' -y -vcodec libx264 -s 1230x566 -pix_fmt yuv420p "+MLE_directory+"MLE_profiles.mp4"
        print (system_call)
        subprocess.call(system_call, shell=True)
        # delete the pngs if you want
        if not keep_pngs:
            system_call = "rm "+MLE_directory+"MLE_profiles*.png"
            subprocess.call(system_call, shell=True)
    plt.close(fig)

def MakeChiPlotsColouredByK(DataDirectory, fname_prefix, basin_list=[0], start_movern=0.2, d_movern=0.1, n_movern=7,size_format='ESURF', FigFormat='png', animate=False, keep_pngs=False, parallel=False):
    """
    This function makes chi-elevation plots for each basin and each value of m/n
    where the channels are coloured by the K value (for model runs with spatially varying K).

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_list: a list of the basins to make the plots for. If an empty list is passed then
        all the basins will be analysed. Default = basin 0.
        start_movern (float): the starting m/n value. Default is 0.2
        d_movern (float): the increment between the m/n values. Default is 0.1
        n_movern (float): the number of m/n values analysed. Default is 7.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        animate (bool): If this is true then it creates a movie of the chi-elevation plots coloured by MLE.
        keep_pngs (bool): If this is false and the animation flag is true, then the pngs are deleted and just the video is kept.

    Returns:
        Plot of each m/n value for each basin.

    Author: FJC
    """
    from LSDPlottingTools import colours

    # check if a directory exists for the chi plots. If not then make it.
    K_directory = DataDirectory+'chi_plots_K/'
    if not os.path.isdir(K_directory):
        os.makedirs(K_directory)

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

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.95,top=1.0)
    ax = fig.add_subplot(gs[10:95,5:80])
    #colorbar axis
    ax2 = fig.add_subplot(gs[10:95,82:85])

    # read in the csv files
    if not parallel:
        ProfileDF = Helper.ReadChiProfileCSV(DataDirectory, fname_prefix)
        BasinStatsDF = Helper.ReadBasinStatsCSV(DataDirectory, fname_prefix)
    else:
        ProfileDF = Helper.AppendMovernCSV(DataDirectory,fname_prefix)
        BasinStatsDF = Helper.AppendBasinCSVs(DataDirectory,fname_prefix)

    # get the number of basins
    basin_keys = list(BasinStatsDF['basin_key'])
    basin_keys = [int(x) for x in basin_keys]

    # get the list of basins
    if basin_list == []:
        print("You didn't give me a list of basins, so I'll just run the analysis on all of them!")
        basin_list = basin_keys

    # loop through each m over n value
    end_movern = start_movern+d_movern*(n_movern-1)
    m_over_n_values = np.linspace(start_movern,end_movern,n_movern)

    # best fit moverns
    best_fit_moverns = SimpleMaxMLECheck(BasinStatsDF)

    for m_over_n in m_over_n_values:
        # read in the full stats file
        print("This m/n is: "+str(m_over_n))
        if not parallel:
            FullStatsDF = Helper.ReadFullStatsCSV(DataDirectory,fname_prefix,m_over_n)
        else:
            FullStatsDF = Helper.AppendFullStatsCSVs(DataDirectory,m_over_n,fname_prefix)

        # loop through all the basins in the basin list
        for basin_key in basin_list:
            print("This basin key is %s") %str(basin_key)

            # mask the data frames for this basin
            ProfileDF_basin = ProfileDF[ProfileDF['basin_key'] == basin_key]
            FullStatsDF_basin = FullStatsDF[FullStatsDF['basin_key'] == basin_key]

            #print FullStatsDF_basin
            print ("Getting the reference_source_key")

            print(FullStatsDF_basin.iloc[0]['reference_source_key'])

            # get the data frame for the main stem
            ProfileDF_MS = ProfileDF_basin[ProfileDF_basin['source_key'] == FullStatsDF_basin.iloc[0]['reference_source_key']]

            # get the data frame for the tributaries
            ProfileDF_basin = ProfileDF_basin[ProfileDF_basin['source_key'] != FullStatsDF_basin.iloc[0]['reference_source_key']]
            # merge with the full data to get the MLE for the tributaries
            ProfileDF_tribs = ProfileDF_basin.merge(FullStatsDF_basin, left_on = "source_key", right_on = "test_source_key")

            # get the chi and elevation data for the main stem
            movern_key = 'm_over_n = %s' %(str(m_over_n))
            MainStemX = list(ProfileDF_MS[movern_key])
            MainStemElevation = list(ProfileDF_MS['elevation'])
            MainStemK = list(ProfileDF_MS['K_value'])

            # get the chi, elevation, and MLE for the tributaries
            TributariesX = list(ProfileDF_tribs[movern_key])
            TributariesElevation = list(ProfileDF_tribs['elevation'])
            TributariesK = list(ProfileDF_tribs['K_value'])

            # get the colourmap to colour channels by the MLE value
            #NUM_COLORS = len(MLE)
            K_array = np.asarray(TributariesK)
            min_K = np.min(K_array)
            max_K = np.max(K_array)
            this_cmap = plt.cm.Spectral
            n_colours = 10
            this_cmap = colours.cmap_discretize(n_colours, this_cmap)
            cNorm  = colors.Normalize(vmin=min_K, vmax=max_K)
            plt.cm.ScalarMappable(norm=cNorm, cmap=this_cmap)

            # now plot the data with a colourmap
            sc = ax.scatter(TributariesX,TributariesElevation,c=TributariesK,cmap=this_cmap, norm=cNorm, s=2.5, edgecolors='none')
            sc = ax.scatter(MainStemX, MainStemElevation,c=MainStemK,cmap=this_cmap, norm=cNorm, s=2.5, edgecolors='none')

            # some formatting of the figure
            ax.spines['top'].set_linewidth(1)
            ax.spines['left'].set_linewidth(1)
            ax.spines['right'].set_linewidth(1)
            ax.spines['bottom'].set_linewidth(1)

            # make the labels
            ax.set_xlabel("$\chi$ (m)")
            ax.set_ylabel("Elevation (m)")

            # the best fit m/n
            best_fit_movern = best_fit_moverns[basin_key]
            print ("BEST FIT M/N IS: "+ str(best_fit_movern))
            print ("THIS M/N IS: "+str(m_over_n))

            # label with the basin and m/n
            title_string = "Basin "+str(basin_key)+", $m/n$ = "+str(m_over_n)
            if best_fit_movern == m_over_n:
                ax.text(0.05, 0.95, title_string,
                        verticalalignment='top', horizontalalignment='left',
                        transform=ax.transAxes,
                        color='red', fontsize=10)
            else:
                ax.text(0.05, 0.95, title_string,
                        verticalalignment='top', horizontalalignment='left',
                        transform=ax.transAxes,
                        color='black', fontsize=10)

            # add the colorbar
            colorbarlabel = "$K$"
            cbar = plt.colorbar(sc,cmap=this_cmap,spacing='uniform', orientation='vertical',cax=ax2)
            cbar.set_label(colorbarlabel, fontsize=10)
            ax2.set_ylabel(colorbarlabel, fontname='Arial', fontsize=10)

            #change labels to scientific notation
            colours.fix_colourbar_ticks(cbar,n_colours, cbar_type=float, min_value = min_K, max_value = max_K, cbar_label_rotation=0, cbar_orientation='vertical')
            # we need to get linear values between min and max K
            these_labels = np.linspace(min_K,max_K,n_colours)
            # now round these and convert to scientific notation
            these_labels = [str('{:.2e}'.format(float(x))) for x in these_labels]
            new_labels = []
            for label in these_labels:
                a,b = label.split("e")
                b = b.replace("0", "")
                new_labels.append(a+' x 10$^{%s}$' % b)

            ax2.set_yticklabels(new_labels, fontsize=8)

            #save the plot
            newFilename = K_directory+"Chi_profiles_by_K_"+str(basin_key)+"_"+str(m_over_n)+"."+str(FigFormat)

            # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
            ax.tick_params(axis='both', width=1, pad = 2)
            for tick in ax.xaxis.get_major_ticks():
                tick.set_pad(2)

            plt.savefig(newFilename,format=FigFormat,dpi=300)
            ax.cla()
            ax2.cla()

    if animate:
        # animate the pngs using ffmpeg
        system_call = "ffmpeg -framerate 3 -pattern_type glob -i '"+K_directory+"Chi_profiles_by_K*.png' -y -vcodec libx264 -s 1230x566 -pix_fmt yuv420p "+K_directory+"Chi_profiles_by_K.mp4"
        print (system_call)
        subprocess.call(system_call, shell=True)
        # delete the pngs if you want
        if not keep_pngs:
            system_call = "rm "+K_directory+"Chi_profiles_by_K*.png"
            subprocess.call(system_call, shell=True)
    plt.close(fig)

def MakeChiPlotsColouredByLith(DataDirectory, fname_prefix, basin_list=[0], start_movern=0.2, d_movern=0.1, n_movern=7,
                    size_format='ESURF', FigFormat='png', animate=False, keep_pngs=False, parallel=False):
    """
    This function makes chi-elevation plots for each basin and each value of m/n
    where the channels are coloured by the K value (for model runs with spatially varying K).

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_list: a list of the basins to make the plots for. If an empty list is passed then
        all the basins will be analysed. Default = basin 0.
        start_movern (float): the starting m/n value. Default is 0.2
        d_movern (float): the increment between the m/n values. Default is 0.1
        n_movern (float): the number of m/n values analysed. Default is 7.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        animate (bool): If this is true then it creates a movie of the chi-elevation plots coloured by MLE.
        keep_pngs (bool): If this is false and the animation flag is true, then the pngs are deleted and just the video is kept.

    Returns:
        Plot of each m/n value for each basin.

    Author: FJC
    """
    from LSDPlottingTools import colours
    print("WARNING DEPRECATED FUNCTION, USE THE MakeChiPlotsByLith FROM LSDMapLithoPlotting")

    # check if a directory exists for the chi plots. If not then make it.
    K_directory = DataDirectory+'chi_plots_Lith/'
    if not os.path.isdir(K_directory):
        os.makedirs(K_directory)

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

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.95,top=1.0)
    ax = fig.add_subplot(gs[10:95,5:80])
    #colorbar axis
    ax2 = fig.add_subplot(gs[10:95,82:85])

    # read in the csv files
    if not parallel:
        ProfileDF = Helper.ReadChiProfileCSV(DataDirectory, fname_prefix)
        BasinStatsDF = Helper.ReadBasinStatsCSV(DataDirectory, fname_prefix)
    else:
        ProfileDF = Helper.AppendMovernCSV(DataDirectory,fname_prefix)
        BasinStatsDF = Helper.AppendBasinCSVs(DataDirectory,fname_prefix)

    # get the number of basins
    basin_keys = list(BasinStatsDF['basin_key'])
    basin_keys = [int(x) for x in basin_keys]

    # get the list of basins
    if basin_list == []:
        print("You didn't give me a list of basins, so I'll just run the analysis on all of them!")
        basin_list = basin_keys

    # loop through each m over n value
    end_movern = start_movern+d_movern*(n_movern-1)
    m_over_n_values = np.linspace(start_movern,end_movern,n_movern)

    # best fit moverns
    best_fit_moverns = SimpleMaxMLECheck(BasinStatsDF)

    for m_over_n in m_over_n_values:
        # read in the full stats file
        print("This m/n is: "+str(m_over_n))
        if not parallel:
            FullStatsDF = Helper.ReadFullStatsCSV(DataDirectory,fname_prefix,m_over_n)
        else:
            FullStatsDF = Helper.AppendFullStatsCSVs(DataDirectory,m_over_n,fname_prefix)

        # loop through all the basins in the basin list
        for basin_key in basin_list:
            print("This basin key is %s") %str(basin_key)

            # mask the data frames for this basin
            ProfileDF_basin = ProfileDF[ProfileDF['basin_key'] == basin_key]
            FullStatsDF_basin = FullStatsDF[FullStatsDF['basin_key'] == basin_key]

            #print FullStatsDF_basin
            print ("Getting the reference_source_key")

            print(FullStatsDF_basin.iloc[0]['reference_source_key'])

            # get the data frame for the main stem
            ProfileDF_MS = ProfileDF_basin[ProfileDF_basin['source_key'] == FullStatsDF_basin.iloc[0]['reference_source_key']]

            # get the data frame for the tributaries
            ProfileDF_basin = ProfileDF_basin[ProfileDF_basin['source_key'] != FullStatsDF_basin.iloc[0]['reference_source_key']]
            # merge with the full data to get the MLE for the tributaries
            ProfileDF_tribs = ProfileDF_basin.merge(FullStatsDF_basin, left_on = "source_key", right_on = "test_source_key")

            # get the chi and elevation data for the main stem
            movern_key = 'm_over_n = %s' %(str(m_over_n))
            MainStemX = list(ProfileDF_MS[movern_key])
            MainStemElevation = list(ProfileDF_MS['elevation'])
            MainStemK = list(ProfileDF_MS[fname_prefix+"_geol"])


            # get the chi, elevation, and MLE for the tributaries
            TributariesX = list(ProfileDF_tribs[movern_key])
            TributariesElevation = list(ProfileDF_tribs['elevation'])
            TributariesK = list(ProfileDF_tribs[fname_prefix+"_geol"])


            # get the colourmap to colour channels by the MLE value
            #NUM_COLORS = len(MLE)
            K_array = np.asarray(TributariesK)
            min_K = np.min(K_array)
            max_K = np.max(K_array)
            this_cmap = plt.cm.Spectral
            n_colours = 10
            this_cmap = colours.cmap_discretize(n_colours, this_cmap)
            cNorm  = colors.Normalize(vmin=min_K, vmax=max_K)
            plt.cm.ScalarMappable(norm=cNorm, cmap=this_cmap)

            # now plot the data with a colourmap
            sc = ax.scatter(TributariesX,TributariesElevation,c=TributariesK,cmap=this_cmap, norm=cNorm, s=2.5, edgecolors='none')
            sc = ax.scatter(MainStemX, MainStemElevation,c=MainStemK,cmap=this_cmap, norm=cNorm, s=2.5, edgecolors='none')

            # some formatting of the figure
            ax.spines['top'].set_linewidth(1)
            ax.spines['left'].set_linewidth(1)
            ax.spines['right'].set_linewidth(1)
            ax.spines['bottom'].set_linewidth(1)

            # make the labels
            ax.set_xlabel("$\chi$ (m)")
            ax.set_ylabel("Elevation (m)")

            # the best fit m/n
            best_fit_movern = best_fit_moverns[basin_key]
            print ("BEST FIT M/N IS: "+ str(best_fit_movern))
            print ("THIS M/N IS: "+str(m_over_n))

            # label with the basin and m/n
            title_string = "Basin "+str(basin_key)+", $m/n$ = "+str(m_over_n)
            if best_fit_movern == m_over_n:
                ax.text(0.05, 0.95, title_string,
                        verticalalignment='top', horizontalalignment='left',
                        transform=ax.transAxes,
                        color='red', fontsize=10)
            else:
                ax.text(0.05, 0.95, title_string,
                        verticalalignment='top', horizontalalignment='left',
                        transform=ax.transAxes,
                        color='black', fontsize=10)

            # add the colorbar
            colorbarlabel = "$Lith$"
            cbar = plt.colorbar(sc,cmap=this_cmap,spacing='uniform', orientation='vertical',cax=ax2)
            cbar.set_label(colorbarlabel, fontsize=10)
            ax2.set_ylabel(colorbarlabel, fontname='Arial', fontsize=10)

            #change labels to scientific notation
            colours.fix_colourbar_ticks(cbar,n_colours, cbar_type=float, min_value = min_K, max_value = max_K, cbar_label_rotation=0, cbar_orientation='vertical')
            # we need to get linear values between min and max K
            these_labels = np.linspace(min_K,max_K,n_colours)
            # now round these and convert to scientific notation
            these_labels = [str('{:.2e}'.format(float(x))) for x in these_labels]
            new_labels = []
            for label in these_labels:
                a,b = label.split("e")
                b = b.replace("0", "")
                new_labels.append(a+' x 10$^{%s}$' % b)

            ax2.set_yticklabels(new_labels, fontsize=8)

            #save the plot
            newFilename = K_directory+"Chi_profiles_by_Lith_"+str(basin_key)+"_"+str(m_over_n)+"."+str(FigFormat)

            # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
            ax.tick_params(axis='both', width=1, pad = 2)
            for tick in ax.xaxis.get_major_ticks():
                tick.set_pad(2)

            plt.savefig(newFilename,format=FigFormat,dpi=300)
            ax.cla()
            ax2.cla()

    if animate:
        # animate the pngs using ffmpeg
        system_call = "ffmpeg -framerate 3 -pattern_type glob -i '"+K_directory+"Chi_profiles_by_Lith*.png' -y -vcodec libx264 -s 1230x566 -pix_fmt yuv420p "+K_directory+"Chi_profiles_by_Lith.mp4"
        print (system_call)
        subprocess.call(system_call, shell=True)
        # delete the pngs if you want
        if not keep_pngs:
            system_call = "rm "+K_directory+"Chi_profiles_by_Lith*.png"
            subprocess.call(system_call, shell=True)
    plt.close(fig)


def PlotProfilesRemovingOutliers(DataDirectory, fname_prefix, basin_list=[0], start_movern=0.2, d_movern=0.1, n_movern=7, size_format = "geomorphology", FigFormat="png", parallel=False):
    """
    This function is used to plot the chi profiles as they have outliers removed.
    It calls thefunction CheckMLEOutliers, which you should read to get details
    on how outliers are calculated and removed

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_list: a list of the basins to make the plots for. If an empty list is passed then
        all the basins will be analysed. Default = basin 0.
        start_movern (float): the starting m/n value. Default is 0.2
        d_movern (float): the increment between the m/n values. Default is 0.1
        n_movern (float): the number of m/n values analysed. Default is 7.

    Returns:
        Plots of chi profiles with basins removed

    Author: SMM
    """

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

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[10:95,5:80])
    #colorbar axis
    ax2 = fig.add_subplot(gs[10:95,82:85])

    # we open the first file just so that we can get a counter list
    if not parallel:
        FirstDF = Helper.ReadFullStatsCSV(DataDirectory,fname_prefix,start_movern)
    else:
        FirstDF = Helper.AppendFullStatsCSVs(DataDirectory,start_movern,fname_prefix)

    # get the number of basins
    basin_keys = list(FirstDF['basin_key'])
    basin_keys = [float(x) for x in basin_keys]

    # get the list of basins
    if basin_list == []:
        print("You didn't give me a list of basins, so I'll just run the analysis on all of them!")
        basin_list = basin_keys
        basin_set = set(basin_list)
        basin_list = list(basin_set)
        basin_list = [int(i) for i in basin_list]

    # First we get all the information about outliers, m/n values and MLE
    # values from the CheckMLEOutliers function
    Outlier_counter, removed_sources_dict, best_fit_movern_dict, MLEs_dict = CheckMLEOutliers(DataDirectory, fname_prefix, basin_list, start_movern, d_movern, n_movern, parallel=parallel)

    # Now get the chi profiles of all the basins and channels
    # Load from file and put into a pandas data frame
    if not parallel:
        ProfileDF = Helper.ReadChiProfileCSV(DataDirectory,fname_prefix)
    else:
        ProfileDF = Helper.AppendMovernCSV(DataDirectory,fname_prefix)

    # now we need to get a plot for each basin, showing the incremental removal of outlying tribs
    for basin_number in basin_list:

        # Get the removed sources indices and the MLEs for this particular basin
        these_removed_sources = removed_sources_dict[basin_number]
        these_MLEs = MLEs_dict[basin_number]

        # loop through the best fit moverns
        # each value represents the MLE for a given number of removed outlying tributaries
        removed_sources_list = []
        for idx,best_fit_movern in enumerate(best_fit_movern_dict[basin_number]):

            print("The best fit m/n is: "+ str(best_fit_movern)+" and the index is "+str(idx))

            # mask the data frames for this basin
            ProfileDF_basin = ProfileDF[ProfileDF['basin_key'] == basin_number]

            # We also need to get the fullstats MLE file. This will be used
            # to colour tribs as well as get the source numbers from
            # the outlier list
            if not parallel:
                FullStatsDF = Helper.ReadFullStatsCSV(DataDirectory,fname_prefix,best_fit_movern)
            else:
                FullStatsDF = Helper.AppendFullStatsCSVs(DataDirectory,best_fit_movern,fname_prefix)

            # mask the data so you only get the correct basin
            FullStatsDF_basin = FullStatsDF[FullStatsDF['basin_key'] == basin_number]

            # extract the relevant data
            #MLE_values = list(FullStatsDF_basin['MLE'])
            #RMSE_values = list(FullStatsDF_basin['RMSE'])
            trib_values = list(FullStatsDF_basin['test_source_key'])
            ref_values = list(FullStatsDF_basin['reference_source_key'])

            # Now get the excluded tributaries for this iteration
            #removed_MLE = these_MLEs[idx]

            # Note that index 0 is the basin with no removed tributaries.
            if idx != 0:
                removed_sources_list.extend(these_removed_sources[idx-1])

            # now you need to get the actual source numbers by indexing into the source list
            the_removed_sources = []
            for source_index in removed_sources_list:
                the_removed_sources.append( trib_values[source_index]  )

            print("The main stem is: ")
            print( ref_values[0])

            print("The removed tribs are: ")
            print(the_removed_sources)



            # get the data frame for the main stem
            # It does this because the main stem source is always the 0 element in the trib_values list
            ProfileDF_MS = ProfileDF_basin[ProfileDF_basin['source_key'] == ref_values[0]]

            # get the data frame for the tributaries
            ProfileDF_basin = ProfileDF_basin[ProfileDF_basin['source_key'] != ref_values[0]]

            # now split the tributaries into exluded and non excluded tribs
            #ProfileDF_outliers = ProfileDF_basin.filter(items=removed_sources_list)
            #ProfileDF_kept = ProfileDF_basin.filter(items=removed_sources_list)
            ProfileDF_outliers = ProfileDF_basin[ProfileDF_basin.source_key.isin(the_removed_sources)]
            ProfileDF_kept = ProfileDF_basin[~ProfileDF_basin.source_key.isin(the_removed_sources)]

            #ProfileDF_basin = ProfileDF_basin[ProfileDF_basin.source_key.isin(removed_sources_list)]
            #ProfileDF_basin = ProfileDF_basin[ProfileDF_basin.source_key.isin([2,3,4])]
            #print(ProfileDF_basin)

            # merge with the full data to get the MLE for the tributaries
            ProfileDF_trib_outliers = ProfileDF_outliers.merge(FullStatsDF_basin, left_on = "source_key", right_on = "test_source_key")
            ProfileDF_trib_kept = ProfileDF_kept.merge(FullStatsDF_basin, left_on = "source_key", right_on = "test_source_key")

            # get the chi and elevation data for the main stem
            movern_key = 'm_over_n = %s' %(str(best_fit_movern))
            MainStemX = list(ProfileDF_MS[movern_key])
            MainStemElevation = list(ProfileDF_MS['elevation'])

            # get the chi, elevation, and MLE for the tributaries
            TributariesX_outliers = list(ProfileDF_trib_outliers[movern_key])
            TributariesElevation_outliers = list(ProfileDF_trib_outliers['elevation'])
            #TributariesMLE_outliers = list(ProfileDF_trib_outliers['MLE'])

            TributariesX_kept = list(ProfileDF_trib_kept[movern_key])
            TributariesElevation_kept = list(ProfileDF_trib_kept['elevation'])
            TributariesMLE_kept = list(ProfileDF_trib_kept['MLE'])

            # now reset the

            # get the colourmap to colour channels by the MLE value
            #NUM_COLORS = len(MLE)
            MLE_array = np.asarray(TributariesMLE_kept)
            this_cmap = plt.cm.Reds
            cNorm  = colors.Normalize(vmin=np.min(MLE_array), vmax=np.max(MLE_array))
            plt.cm.ScalarMappable(norm=cNorm, cmap=this_cmap)

            # now plot the data with a colourmap
            sc = ax.scatter(TributariesX_kept,TributariesElevation_kept,c=TributariesMLE_kept,cmap=this_cmap, norm=cNorm, s=2.5, edgecolors='none')

            # Add the outliers if the basin has them
            if(len(removed_sources_list)>0):
                ax.scatter(TributariesX_outliers,TributariesElevation_outliers,c="b", norm=cNorm, s=2.5, edgecolors='none', alpha = 0.3)

            ax.plot(MainStemX,MainStemElevation,lw=2, c='k')

            # some formatting of the figure
            ax.spines['top'].set_linewidth(1)
            ax.spines['left'].set_linewidth(1)
            ax.spines['right'].set_linewidth(1)
            ax.spines['bottom'].set_linewidth(1)

            # make the lables
            ax.set_xlabel("$\chi$ (m)")
            ax.set_ylabel("Elevation (m)")

            # label with the basin and m/n
            title_string = "Basin "+str(basin_number)+", best fit "+r'$\theta$' + " = str(best_fit_movern)"
            ax.text(0.05, 0.95, title_string,
                    verticalalignment='top', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='black', fontsize=10)

            # add the colorbar
            colorbarlabel = "$MLE$"
            cbar = plt.colorbar(sc,cmap=this_cmap,spacing='uniform', orientation='vertical',cax=ax2)
            cbar.set_label(colorbarlabel, fontsize=10)
            ax2.set_ylabel(colorbarlabel, fontname='Arial', fontsize=10)

            #save the plot
            newFilename = DataDirectory+"MLE_profiles"+str(basin_number)+"_"+str(best_fit_movern)+"_removed_"+str(idx)+".png"

            # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
            ax.tick_params(axis='both', width=1, pad = 2)
            for tick in ax.xaxis.get_major_ticks():
                tick.set_pad(2)

            plt.savefig(newFilename,format=FigFormat,dpi=300)
            ax.cla()
            ax2.cla()

def PlotMLEWithMOverN(DataDirectory, fname_prefix, basin_list = [0], size_format='ESURF', FigFormat='png', start_movern=0.2, d_movern = 0.1, n_movern = 7, parallel=False):
    """
    This function makes a plot of the MLE values for each m/n showing how the MLE values change
    as you remove the tributaries.

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_list: a list of the basins to make the plots for. If an empty list is passed then
        all the basins will be analysed. Default = basin 0.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        start_movern (float): the starting m/n value. Default is 0.2
        d_movern (float): the increment between the m/n values. Default is 0.1
        n_movern (float): the number of m/n values analysed. Default is 7.

    Returns:
        Plots of MLE values for each m/n

    Author: FJC
    """
    from cycler import cycler
    from matplotlib import lines
    import matplotlib.patches as patches

    # check if a directory exists for the MLE plots. If not then make it.
    MLE_directory = DataDirectory+'MLE_plots/'
    if not os.path.isdir(MLE_directory):
        os.makedirs(MLE_directory)

    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # set colours for tributaries
    # plt.rc('axes', prop_cycle=(cycler('color', ['k', '0.5', '0.5', '0.5', '0.5']) +
    #                        cycler('linestyle', ['-', '--', ':', '-.', '--'])))

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

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.85,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])

    # we open the first file just so that we can get a counter list
    if not parallel:
        FirstDF = Helper.ReadFullStatsCSV(DataDirectory,fname_prefix,start_movern)
    else:
        FirstDF = Helper.AppendFullStatsCSVs(DataDirectory,start_movern,fname_prefix)

    # get the list of m over n values
    end_movern = start_movern+d_movern*(n_movern-1)
    m_over_n_values = np.linspace(start_movern,end_movern,n_movern)

    # get the number of basins
    basin_keys = list(FirstDF['basin_key'])
    basin_keys = [float(x) for x in basin_keys]

    # get the list of basins
    if basin_list == []:
        print("You didn't give me a list of basins, so I'll just run the analysis on all of them!")
        basin_list = basin_keys
        basin_set = set(basin_list)
        basin_list = list(basin_set)
        basin_list = [int(i) for i in basin_list]

    # First we get all the information about outliers, m/n values and MLE
    # values from the CheckMLEOutliers function
    Outlier_counter, removed_sources_dict, best_fit_movern_dict, MLEs_dict = CheckMLEOutliers(DataDirectory, fname_prefix, basin_list, start_movern, d_movern, n_movern, parallel=parallel)

    # Now get the chi profiles of all the basins and channels
    # Load from file and put into a pandas data frame
    if not parallel:
        ProfileDF = Helper.ReadChiProfileCSV(DataDirectory,fname_prefix)
    else:
        ProfileDF = Helper.AppendMovernCSV(DataDirectory,fname_prefix)

    # get list of line styles for plotting. This is hacky but not sure how else to do this.
    ls = lines.lineStyles.keys()
    ls = ls[3:]
    ls = ls + ls

    # loop through each basin and each number of removed tributaries
    for basin_number in basin_list:

        print ("This basin is: " +str(basin_number))

        # Get the removed sources indices and the MLEs for this particular basin
        basin_removed_sources = removed_sources_dict[basin_number]
        n_removed_sources = len(basin_removed_sources)
        basin_MLEs = MLEs_dict[basin_number]

        # get the best fit m over ns for this basin
        best_fit_moverns = best_fit_movern_dict[basin_number]
        print (best_fit_moverns)

        # colours for each iteration
        #colours = np.linspace(0.9,0.2,len(n_removed_sources)+1)


        # loop through the number of removed tributaires and get the MLE for each m/n for each iteration
        for i in range(n_removed_sources+1):

            # get the MLE of the best fit m over n
            best_fit_movern = best_fit_moverns[i]
            # get the index in the MLE list
            idx = int(round((best_fit_movern - start_movern)/d_movern,0))
            print (idx)
            best_fit_MLE = basin_MLEs[idx][i]
            print ("The best fit MLE is: "+str(best_fit_MLE)+", where "+ r'$\theta$'  ' = '+str(best_fit_movern))

            # get the MLEs for this iteration
            these_MLEs = basin_MLEs[:,i]

            # get the ratio of these MLEs to the best fit
            ratio_MLEs = [x/best_fit_MLE for x in these_MLEs]

            # no removed tributaries
            if i == 0:
                # plot the data
                ax.scatter(m_over_n_values,ratio_MLEs, label = str(i), c='k', s=5, zorder=100)
                ax.plot(m_over_n_values,ratio_MLEs, c='0.75', ls="--")

                # get the limits for the arrow
                max_MLE = max(ratio_MLEs)
                min_MLE = min(ratio_MLEs)
                dy = (max_MLE-min_MLE)/8
                spacing = 1.5
                # add arrow at best fit m/n
                ax.add_patch(
                    patches.Arrow(
                        best_fit_movern, #x
                        max_MLE-(dy*spacing), #y
                        0, #dx
                        dy, #dy
                        width = 0.05,
                        facecolor = 'r',
                        edgecolor = 'r')
                    )
                ax.text(best_fit_movern-0.09, max_MLE-1.5*(dy*spacing), "Best-fit " + r"$\theta$ = "+str(best_fit_movern),fontsize=8, color="r")
            #remove tribs
            #else:
                # plot the data
                #ax.scatter(m_over_n_values,ratio_MLEs, label = str(i), s=5, c=colours[i]) # different linestyle for each iteration?

        # set the axes labels
        ax.set_xlabel(r'$\theta$')
        ax.set_ylabel('MLE')

        # set the ylim
        ax.set_ylim(min_MLE,max_MLE+(dy*spacing))

        # add the legend
        #ax.legend(loc='right', bbox_to_anchor=(1.25,0.5), title = 'Iterations', frameon=False)

        # some formatting of the figure
        ax.spines['top'].set_linewidth(1)
        ax.spines['left'].set_linewidth(1)
        ax.spines['right'].set_linewidth(1)
        ax.spines['bottom'].set_linewidth(1)

        # label with the basin and m/n
        best_fit_movern = best_fit_movern_dict[basin_number][0]
        #title_string = "Basin "+str(basin_number)+"; Best fit $m/n$: "+str(best_fit_movern)
        #ax.set_title(title_string)
        # ax.text(0, 1.1, title_string,
        #     verticalalignment='top', horizontalalignment='left',
        #     transform=ax.transAxes,
        #     color='red', fontsize=10)

        #save the plot
        newFilename = MLE_directory+"MLE_fxn_movern_"+str(basin_number)+"."+FigFormat

        # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
        ax.tick_params(axis='both', width=1, pad = 2)
        for tick in ax.xaxis.get_major_ticks():
            tick.set_pad(2)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base=0.1))

        plt.savefig(newFilename,format=FigFormat,dpi=300)
        ax.cla()
    plt.close(fig)

def MakeMOverNSummaryPlot(DataDirectory, fname_prefix, basin_list=[], start_movern=0.2, d_movern=0.1, n_movern=7, size_format='ESURF', FigFormat='png', SA_channels=False, show_legend=True,parallel=False,Chi_disorder=False):
    """
    This function makes a summary plot of the best fit m/n from the different
    methods.

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_list: list of basin keys to analyse, default = [] (all basins)
        start_movern (float): the starting m/n value. Default is 0.2
        d_movern (float): the increment between the m/n values. Default is 0.1
        n_movern (float): the number of m/n values analysed. Default is 7.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        SA_channels (bool): If true, will include the SA data separated by channel
        show_legend (bool): If true, display the legend with the plot
        Chi_disorder (bool): if true, will include the data from the chi disorder method

    Returns:
        Makes a summary plot

    Author: FJC
    """
    # check if a directory exists for the summary plots. If not then make it.
    summary_directory = DataDirectory+'summary_plots/'
    if not os.path.isdir(summary_directory):
        os.makedirs(summary_directory)

    from matplotlib.ticker import FuncFormatter, MaxNLocator
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

    print("SHOW LEGEND", show_legend)

    if show_legend:
        gs = plt.GridSpec(100,100,bottom=0.15,left=0.05,right=0.75,top=0.9)
    else:
        gs = plt.GridSpec(100,100,bottom=0.15,left=0.05,right=0.95,top=0.95)

    ax = fig.add_subplot(gs[5:100,10:95])

    # read in the summary csv
    df = Helper.ReadMOverNSummaryCSV(summary_directory,fname_prefix)
    print (df)

    if basin_list != []:
        basin_keys = basin_list
    else:
        # get the basin keys
        basin_keys = df['basin_key'].tolist()
        print (basin_keys)

    df = df[df['basin_key'].isin(basin_keys)]

    # plot the full chi data
    full_chi_keys = df['basin_key'].as_matrix()-0.2
    ax.scatter(full_chi_keys, df['Chi_MLE_full'],marker='o', edgecolors='k', lw=0.5, facecolors='#e34a33', s=15, zorder=200, label='Chi all data')

    # plot the points data
    median_movern = df['Chi_MLE_points'].as_matrix()
    points_max_err = df['Chi_MLE_points_max'].as_matrix()
    points_max_err = points_max_err-median_movern
    points_min_err = df['Chi_MLE_points_min'].as_matrix()
    points_min_err = median_movern-points_min_err
    errors = np.array(zip(points_min_err, points_max_err)).T

    points_chi_keys = df['basin_key'].as_matrix()-0.1
    ax.errorbar(points_chi_keys, df['Chi_MLE_points'], s=15, marker='o', xerr=None, yerr=errors, ecolor='#fdbb84', fmt='none', elinewidth=1,label='_nolegend_')
    ax.scatter(points_chi_keys, df['Chi_MLE_points'], s=15, c='#fdbb84', marker='o', edgecolors='k', lw=0.5,facecolors='#fdbb84', label='Chi bootstrap',zorder=200)

    # plot the chi disorder data if you want it
    if Chi_disorder:
        median_movern = df['Chi_disorder'].as_matrix()
        points_max_err = df['Chi_disorder_max'].as_matrix()
        points_max_err = points_max_err-median_movern
        points_min_err = df['Chi_disorder_min'].as_matrix()
        points_min_err = median_movern-points_min_err
        errors = np.array(zip(points_min_err, points_max_err)).T

        disorder_chi_keys = df['basin_key'].as_matrix()-0.3
        ax.errorbar(disorder_chi_keys, df['Chi_disorder'], s=15, marker='o', xerr=None, yerr=errors, ecolor='#F06292', fmt='none', elinewidth=1,label='_nolegend_')
        ax.scatter(disorder_chi_keys, df['Chi_disorder'],marker='o', edgecolors='k', lw=0.5, facecolors='#F06292', s=15, zorder=100, label='Chi disorder')


    # plot the SA data
    SA_keys = df['basin_key'].as_matrix()
    SA_sterr = df['SA_raw_sterr'].as_matrix()
    ax.errorbar(SA_keys, df['SA_raw'], yerr=SA_sterr, c='#2b8cbe', elinewidth=1, fmt='none',label='_nolegend_')
    ax.scatter(SA_keys, df['SA_raw'], s=15, c='#2b8cbe', edgecolors='k',lw=0.5, label='S-A all data', zorder=100)

    if SA_channels:
        # plot the SA data by tribs
        median_movern = df['SA_tribs'].as_matrix()
        points_max_err = df['SA_tribs_max'].as_matrix()
        points_max_err = points_max_err-median_movern
        points_min_err = df['SA_tribs_min'].as_matrix()
        points_min_err = median_movern-points_min_err
        errors = np.array(zip(points_min_err, points_max_err)).T

        SA_tribs_keys = df['basin_key'].as_matrix()+0.1
        ax.errorbar(SA_tribs_keys, df['SA_tribs'], s=15, marker='D', facecolors='white', xerr=None, yerr=errors, edgecolors='r', fmt='none', elinewidth=1, linestyle = ":", ecolor='r',label='_nolegend_')
        ax.scatter(SA_tribs_keys, df['SA_tribs'], s=15, marker='D', facecolors='white', edgecolors='r', label='S-A by channel',zorder=100)

    # plot the segmented SA data
    median_movern = df['SA_segments'].as_matrix()
    points_max_err = df['SA_segments_max'].as_matrix()
    points_max_err = points_max_err-median_movern
    points_min_err = df['SA_segments_min'].as_matrix()
    points_min_err = median_movern-points_min_err
    errors = np.array(zip(points_min_err, points_max_err)).T

    SA_segment_keys = df['basin_key'].as_matrix()+0.2
    ax.errorbar(SA_segment_keys, df['SA_segments'], s=15, marker='o', facecolors='#a6bddb', xerr=None, yerr=errors, edgecolors='#a6bddb', fmt='none', elinewidth=1, linestyle = ":", ecolor='#a6bddb',label='_nolegend_')
    ax.scatter(SA_segment_keys, df['SA_segments'], s=15, marker='o', facecolors='#a6bddb', edgecolors='k', lw=0.5, label='Segmented S-A', zorder=100)

    # set the axis labels
    ax.set_xlabel('Basin key')
    ax.set_ylabel('Best fit '+r'$\theta$')

    if show_legend:
        print ("ADDING THE LEGEND")
        # sort both labels and handles by labels
        handles, labels = ax.get_legend_handles_labels()
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
        # add the legend
        ax.legend(handles, labels,fontsize=8, bbox_to_anchor=(1.0,0.7),bbox_transform=plt.gcf().transFigure)

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
    ax.tick_params(axis='both', width=1, pad = 2)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(2)

    # change tick spacing
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=2))
    ax.set_xlim(-1,)
    #ax.yaxis.set_major_locator(ticker.MultipleLocator(base=d_movern))

    # remove first tick from the x axis
    #ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))

    #set y axis lims
    #end_movern = end_movern = start_movern+d_movern*(n_movern-1)
    #ax.set_ylim(start_movern-d_movern,1.3)

    # change y axis to the moverns tested
    # end_movern = end_movern = start_movern+d_movern*(n_movern-1)
    # print end_movern
    # ax.yaxis.set_ticks(np.arange(start_movern, end_movern, d_movern))

    newFilename = summary_directory+fname_prefix+"_movern_summary."+FigFormat
    if SA_channels:
        newFilename = summary_directory+fname_prefix+"_movern_summary_with_SA_tribs."+FigFormat
    plt.savefig(newFilename,format=FigFormat,dpi=300)
    ax.cla()
    plt.close(fig)

def MakeMOverNPlotOneMethod(DataDirectory, fname_prefix, basin_list=[], start_movern=0.2, d_movern=0.1, n_movern=7, size_format='ESURF', FigFormat='png', movern_method='chi_points'):
        """
        This function makes a summary plot of the best fit m/n, you choose which method
        you want to plot.

        Args:
            DataDirectory (str): the data directory with the m/n csv files
            fname_prefix (str): The prefix for the m/n csv files
            basin_list: list of basin keys to analyse, default = [] (all basins)
            start_movern (float): the starting m/n value. Default is 0.2
            d_movern (float): the increment between the m/n values. Default is 0.1
            n_movern (float): the number of m/n values analysed. Default is 7.
            size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
            FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
            movern_method (str): the method you want to plot. Can be 'chi_all', 'chi_points', 'SA_raw', or 'SA_segments'. Default 'chi_points'
        Returns:
            Makes a summary plot

        Author: FJC
        """
        # check if a directory exists for the summary plots. If not then make it.
        summary_directory = DataDirectory+'summary_plots/'
        if not os.path.isdir(summary_directory):
            os.makedirs(summary_directory)

        from matplotlib.ticker import FuncFormatter, MaxNLocator
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

        gs = plt.GridSpec(100,100,bottom=0.15,left=0.05,right=0.95,top=0.95)

        ax = fig.add_subplot(gs[5:100,10:95])

        # read in the summary csv
        df = Helper.ReadMOverNSummaryCSV(summary_directory,fname_prefix)
        print (df)

        if basin_list != []:
            basin_keys = basin_list
        else:
            # get the basin keys
            basin_keys = df['basin_key'].tolist()
            print (basin_keys)

        df = df[df['basin_key'].isin(basin_keys)]

        # plot the full chi data
        if movern_method=='chi_all':
            full_chi_keys = df['basin_key'].as_matrix()-0.2
            ax.scatter(full_chi_keys, df['Chi_MLE_full'],marker='o', edgecolors='k', lw=0.5, facecolors='#e34a33', s=15, zorder=200, label='Chi all data')

        elif movern_method=='chi_points':
            # plot the points data
            median_movern = df['Chi_MLE_points'].as_matrix()
            points_max_err = df['Chi_MLE_points_max'].as_matrix()
            points_max_err = points_max_err-median_movern
            points_min_err = df['Chi_MLE_points_min'].as_matrix()
            points_min_err = median_movern-points_min_err
            errors = np.array(zip(points_min_err, points_max_err)).T

            points_chi_keys = df['basin_key'].as_matrix() # -0.1 I removed that to fix the movern plots - Boris
            ax.errorbar(points_chi_keys, df['Chi_MLE_points'], s=15, marker='o', xerr=None, yerr=errors, ecolor='#fdbb84', fmt='none', elinewidth=1,label='_nolegend_')
            ax.scatter(points_chi_keys, df['Chi_MLE_points'], s=15, c='#fdbb84', marker='o', edgecolors='k', lw=0.5,facecolors='#fdbb84', label='Chi bootstrap',zorder=200)

            # plot the SA data
        elif movern_method=='SA_raw':
                SA_keys = df['basin_key'].as_matrix()
                SA_sterr = df['SA_raw_sterr'].as_matrix()
                ax.errorbar(SA_keys, df['SA_raw'], yerr=SA_sterr, c='#2b8cbe', elinewidth=1, fmt='none',label='_nolegend_')
                ax.scatter(SA_keys, df['SA_raw'], s=15, c='#2b8cbe', edgecolors='k',lw=0.5, label='S-A all data', zorder=100)

        else:
            # plot the segmented SA data
            median_movern = df['SA_segments'].as_matrix()
            points_max_err = df['SA_segments_max'].as_matrix()
            points_max_err = points_max_err-median_movern
            points_min_err = df['SA_segments_min'].as_matrix()
            points_min_err = median_movern-points_min_err
            errors = np.array(zip(points_min_err, points_max_err)).T

            SA_segment_keys = df['basin_key'].as_matrix()+0.2
            ax.errorbar(SA_segment_keys, df['SA_segments'], s=15, marker='o', facecolors='#a6bddb', xerr=None, yerr=errors, edgecolors='#a6bddb', fmt='none', elinewidth=1, linestyle = ":", ecolor='#a6bddb',label='_nolegend_')
            ax.scatter(SA_segment_keys, df['SA_segments'], s=15, marker='o', facecolors='#a6bddb', edgecolors='k', lw=0.5, label='Segmented S-A', zorder=100)

        # set the axis labels
        ax.set_xlabel('Basin key')
        ax.set_ylabel('Best fit '+r'$\theta$')

        # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
        ax.tick_params(axis='both', width=1, pad = 2)
        for tick in ax.xaxis.get_major_ticks():
            tick.set_pad(2)

        # change tick spacing
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1))
        ax.set_xlim(-0.5,df['basin_key'].max()+0.5)

        newFilename = summary_directory+fname_prefix+"_movern_"+movern_method+"."+FigFormat
        plt.savefig(newFilename,format=FigFormat,dpi=300)
        ax.cla()
        plt.close(fig)

def MakeMOverNSummaryHistogram(DataDirectory, fname_prefix, basin_list=[], size_format='ESURF', FigFormat='png', start_movern=0.1, n_movern=7, d_movern=0.1, mn_method = "Chi", show_legend=True, parallel=False, Chi_disorder=False):
    """
    This function makes a joyplot showing the distribution of m/n values for each method

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        start_movern (float): the starting m/n value. Default is 0.2
        d_movern (float): the increment between the m/n values. Default is 0.1
        n_movern (float): the number of m/n values analysed. Default is 7.
        mn_method (str): the method of calculating m/n, either "Chi", or "SA". Default "Chi"
        show_legend (bool): if true, add a legend to the plot (default=True)
        Chi_disorder(bool): if true, add the chi disorder analysis

    Returns:
        Histogram of the m/n values

    Author: FJC
    """

    # check if a directory exists for the summary plots. If not then make it.
    summary_directory = DataDirectory+'summary_plots/'
    if not os.path.isdir(summary_directory):
        os.makedirs(summary_directory)

    from matplotlib.ticker import FuncFormatter, MaxNLocator
    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        #fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        figsize=(6.25,3.5)
        #l_pad = -40
    elif size_format == "big":
        #fig = plt.figure(1, facecolor='white',figsize=(16,9))
        figsize=(16,9)
        #l_pad = -50
    else:
        #fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
        figsize=(4.92126,3.2)
        #l_pad = -35

    # if show_legend:
    #     gs = plt.GridSpec(100,100,bottom=0.15,left=0.05,right=0.75,top=0.9)
    # else:
    #     gs = plt.GridSpec(100,100,bottom=0.1,left=0.05,right=0.95,top=0.95)
    #
    # ax = fig.add_subplot(gs[5:100,10:95])

    # read in the summary csv
    if not parallel:
      df = Helper.ReadMOverNSummaryCSV(summary_directory,fname_prefix)
    else:
      print ("MakeMOverNSummaryHistogram not parallelised yet")
      return

    # get the basin keys
    basin_keys = df['basin_key'].tolist()
    print (basin_keys)
    if Chi_disorder:
        columns = ['Chi_MLE_full', 'Chi_MLE_points', 'Chi_disorder', 'SA_raw', 'SA_segments']
        these_labels = ['Chi all data', 'Chi bootstrap', 'Chi disorder', 'S-A all data', 'Segmented S-A']
        colours = ['#e34a33', '#fdbb84', '#F06292', '#2b8cbe', '#a6bddb']

    else:
        columns = ['Chi_MLE_full', 'Chi_MLE_points', 'SA_raw', 'SA_segments']
        these_labels = ['Chi all data', 'Chi bootstrap', 'S-A all data', 'Segmented S-A']
        colours = ['#e34a33', '#fdbb84', '#2b8cbe', '#a6bddb']
    x_spacing = 0.1
    fig, ax = joyplot.joyplot(df, figsize=figsize, column=columns, label_strings=these_labels, x_range=[0,1],grid="x",color=colours,x_title='Best fit '+r'$\theta$' +' distribution',x_spacing=x_spacing)
    #plt.xlabel('Best fit $m/n$')

    newFilename = summary_directory+fname_prefix+"_movern_hist."+FigFormat

    plt.savefig(newFilename,format=FigFormat,dpi=300)
    plt.close(fig)

def PlotMOverNByBasin(DataDirectory, fname_prefix, basin_list = [], size_format='ESURF', FigFormat='png', colour_points=False, column_header="Lithology"):
    """
    This function makes a plot of the best-fit m/n value calculated using the chi Monte Carlo points analysis
    for each basin. The points are coloured by another value which is specified by the user (e.g. lithology).
    The column header must be specified by the user.

    This is a work in progress. I haven't tested it yet and wrote it on the plane.
    I need to think about what is the best way to read in the lithology - will we have a string
    for each lithology? Then we want to map each string to a different colour.

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_list: a list of the basins to analyse
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        colour_points (bool): boolean to decide whether to colour the points by another value, e.g. lithology. Default = False
        column_header (str): if colour_points=True, then use this to tell us what column header in the csv you want to use

    Returns:
        m/n plot by basin_list

    Author: FJC
    """

    # check if a directory exists for the summary plots. If not then make it.
    summary_directory = DataDirectory+'summary_plots/'
    if not os.path.isdir(summary_directory):
        os.makedirs(summary_directory)

    # read in the summary csv file
    df = Helper.ReadMOverNSummaryCSV(DataDirectory,fname_prefix)

    # get the basin keys
    if basin_list == []:
        basin_list = df['basin_keys'].tolist()

    # mask the dataframe for the required basin keys
    df = df[df['basin_key'].isin(basin_list)]

    columns = df.columns()

    # set up the figure
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

    gs = plt.GridSpec(100,100,bottom=0.1,left=0.05,right=0.95,top=0.95)
    ax = fig.add_subplot(gs[5:100,10:95])

    # check if the user wants to colour the points by a specific columns
    if colour_points:
        if not column_header.isin(columns):
            print ("The column header name does not exist! Please check your dataframe")
        else:
            values_to_colour = df[column_header].tolist()
            this_cmap = plt.cm.Set2
            ax.errorbar(basin_list, df['Chi_MLE_points'], s=20, marker='o', xerr=None, yerr=errors, ecolor=values_to_colour, cmap=this_cmap, fmt='none', elinewidth=1,label='_nolegend_')
            ax.scatter(basin_list, df['Chi_MLE_points'], s=20, c=values_to_colour, cmap=this_cmap, marker='o', edgecolors=None, lw=0.5, label='Chi bootstrap',zorder=200)
    else:
        ax.errorbar(basin_list, df['Chi_MLE_points'], s=20, marker='o', xerr=None, yerr=errors, ecolor='#fdbb84', fmt='none', elinewidth=1,label='_nolegend_')
        ax.scatter(basin_list, df['Chi_MLE_points'], s=20, c='#fdbb84', marker='o', edgecolors='k', lw=0.5,facecolors='#fdbb84', label='Chi bootstrap',zorder=200)

    # set the axis labels
    ax.set_xlabel('Basin key')
    ax.set_ylabel('Best fit '+r'$\theta$')

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
    ax.tick_params(axis='both', width=1, pad = 2)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(2)

    # change tick spacing
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1))

    newFilename = summary_directory+fname_prefix+"_movern_summary."+FigFormat
    if SA_channels:
        newFilename = summary_directory+fname_prefix+"_movern_summary_with_SA_tribs."+FigFormat
    plt.savefig(newFilename,format=FigFormat,dpi=300)
    ax.cla()
    plt.close(fig)



# def MakeBasinJoyplot(DataDirectory, fname_prefix, basin_list=[], size_format='ESURF', FigFormat='png'):
#     """
#     This function makes a joyplot showing the m/n value for a list of basins. m/n value
#     is calculated using the chi points method.
#
#     Args:
#         DataDirectory (str): the data directory with the m/n csv files
#         fname_prefix (str): The prefix for the m/n csv files
#         basin_list: a list of the basins to analyse
#         size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
#         FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
#
#     Returns:
#         Joyplot of the m/n values for each basin
#
#     Author: FJC
#     """
#
#     # check if a directory exists for the summary plots. If not then make it.
#     summary_directory = DataDirectory+'summary_plots/'
#     if not os.path.isdir(summary_directory):
#         os.makedirs(summary_directory)
#
#     from matplotlib.ticker import FuncFormatter, MaxNLocator
#     # Set up fonts for plots
#     label_size = 10
#     rcParams['font.family'] = 'sans-serif'
#     rcParams['font.sans-serif'] = ['arial']
#     rcParams['font.size'] = label_size
#
#     # make a figure
#     if size_format == "geomorphology":
#         #fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
#         figsize=(6.25,3.5)
#         #l_pad = -40
#     elif size_format == "big":
#         #fig = plt.figure(1, facecolor='white',figsize=(16,9))
#         figsize=(16,9)
#         #l_pad = -50
#     else:
#         #fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
#         figsize=(4.92126,3.2)
#         #l_pad = -35
#
#     # if show_legend:
#     #     gs = plt.GridSpec(100,100,bottom=0.15,left=0.05,right=0.75,top=0.9)
#     # else:
#     #     gs = plt.GridSpec(100,100,bottom=0.1,left=0.05,right=0.95,top=0.95)
#     #
#     # ax = fig.add_subplot(gs[5:100,10:95])
#
#     # read in the summary csv
#     df = Helper.ReadMOverNSummaryCSV(summary_directory,fname_prefix)
#     print df
#
#     # get the basin keys
#     if basin_list == []:
#         print "You didn't give me a list of basins so I'll analyse all of them"
#         basin_list = df['basin_key'].tolist()
#
#     df = df[df['basin_key'].isin(basin_list)]
#
#     columns = ['basin_key']
#     #these_labels = ['Chi all data', 'Chi Monte Carlo', 'S-A all data', 'Segmented S-A']
#     #colours = ['#e34a33', '#fdbb84', '#2b8cbe', '#a6bddb']
#     x_spacing = 0.2
#     fig, ax = joyplot.joyplot(df, figsize=figsize, column=columns, x_range=[0,1],grid="x",x_title='Best fit $m/n$ distribution',x_spacing=x_spacing)
#     #plt.xlabel('Best fit $m/n$')
#
#     newFilename = summary_directory+fname_prefix+"_basin_joyplots."+FigFormat
#
#     plt.savefig(newFilename,format=FigFormat,dpi=300)
#     plt.close(fig)

#=============================================================================
# RASTER PLOTTING FUNCTIONS
# Functions that interface with LSDMapFigure to plot the m/n analysis with
# spatial data.
#=============================================================================
def MakeRasterPlotsBasins(DataDirectory, fname_prefix, size_format='ESURF', FigFormat='png', parallel=False):
    """
    This function makes a shaded relief plot of the DEM with the basins coloured
    by the basin ID.

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.

    Returns:
        Shaded relief plot with the basins coloured by basin ID

    Author: FJC
    """
    # check if a directory exists for the chi plots. If not then make it.
    raster_directory = DataDirectory+'raster_plots/'
    if not os.path.isdir(raster_directory):
        os.makedirs(raster_directory)

    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # set figure sizes based on format
    if size_format == "geomorphology":
        fig_width_inches = 6.25
    elif size_format == "big":
        fig_width_inches = 16
    else:
        fig_width_inches = 4.92126

    # get the basin IDs to make a discrete colourmap for each ID
    if not parallel:
        BasinInfoDF = Helper.ReadBasinInfoCSV(DataDirectory, fname_prefix)
    else:
        BasinInfoDF = Helper.AppendBasinInfoCSVs(DataDirectory,fname_prefix)

    print (BasinInfoDF)

    basin_keys = list(BasinInfoDF['basin_key'])
    basin_keys = [int(x) for x in basin_keys]

    basin_junctions = list(BasinInfoDF['outlet_junction'])
    basin_junctions = [int(x) for x in basin_junctions]

    print ('Basin keys are: ')
    print (basin_keys)

    # get a discrete colormap
    cmap = plt.cm.jet

    # going to make the basin plots - need to have bil extensions.
    print("I'm going to make the basin plots. Your topographic data must be in ENVI bil format or I'll break!!")

    # get the rasters
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext
    BasinsName = fname_prefix+'_AllBasins'+raster_ext

    # create the map figure
    MF = MapFigure(HillshadeName, DataDirectory,coord_type="UTM_km", colourbar_location='None')

    # add the basins drape
    MF.add_basin_plot(BasinsName,fname_prefix,DataDirectory,
                      use_keys_not_junctions = True, show_colourbar = True,
                      discrete_cmap=True, n_colours=len(basin_keys), colorbarlabel = "Basin ID",
                      colourmap = cmap, adjust_text = False, parallel=parallel)

    # add the basin outlines ### need to parallelise
    if not parallel:
      Basins = LSDP.GetBasinOutlines(DataDirectory, BasinsName)
    else:
      Basins = LSDP.GetMultipleBasinOutlines(DataDirectory)

    MF.plot_polygon_outlines(Basins, linewidth=0.8)

    # add the channel network
    if not parallel:
        ChannelDF = Helper.ReadChiDataMapCSV(DataDirectory,fname_prefix)
    else:
        ChannelDF = Helper.AppendChiDataMapCSVs(DataDirectory,fname_prefix)
    ChannelPoints = LSDP.LSDMap_PointData(ChannelDF, data_type = "pandas", PANDEX = True)
    MF.add_point_data(ChannelPoints,show_colourbar="False", scale_points=True, column_for_scaling='drainage_area',alpha=0.5,zorder=100)

    # add the basin labelling
    label_dict = dict(zip(basin_junctions,basin_keys))

    if not parallel:
      Points = LSDP.GetPointWithinBasins(DataDirectory, BasinsName)
    else:
      Points = LSDP.GetPointsWithinMultipleBasins(DataDirectory, BasinsName)

    MF.add_text_annotation_from_shapely_points(Points, text_colour='k', label_dict=label_dict, zorder=200)

    # Save the figure
    ImageName = raster_directory+fname_prefix+'_basin_keys.'+FigFormat
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 300)

def MakeRasterPlotsMOverN(DataDirectory, fname_prefix, start_movern=0.2, n_movern=7, d_movern=0.1, movern_method='Chi_full', size_format='ESURF', FigFormat='png',lith = False, parallel=False):
    """
    This function makes a shaded relief plot of the DEM with the basins coloured
    by the best fit m/n

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        start_movern (int): The starting m/n, default = 0.2
        n_movern (int): The number of m/n values tested, default = 7.
        movern_method: the method of estimating m over n. Options are full chi "Chi_full", points "Chi_points", or slope-area "SA". Default is "Chi_full".
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        lith (bool): the lithologic information if available. You need to follow the lithologic documentation to generate the right files. (This last does not exists yet).

    Returns:
        Shaded relief plot with the basins coloured by best fit m/n

    Author: FJC
    """
    # check if a directory exists for the chi plots. If not then make it.
    raster_directory = DataDirectory+'raster_plots/'
    if not os.path.isdir(raster_directory):
        os.makedirs(raster_directory)

    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # set figure sizes based on format
    if size_format == "geomorphology":
        fig_width_inches = 6.25
    elif size_format == "big":
        fig_width_inches = 16
    else:
        fig_width_inches = 4.92126

    # get the basin IDs to make a discrete colourmap for each ID
    if not parallel:
        BasinDF = Helper.ReadBasinStatsCSV(DataDirectory,fname_prefix)
        # get the basin IDs to make a discrete colourmap for each ID
        BasinInfoDF = Helper.ReadBasinInfoCSV(DataDirectory, fname_prefix)
    else:
        BasinDF = Helper.AppendBasinCSVs(DataDirectory,fname_prefix)
        BasinInfoDF = Helper.AppendBasinInfoCSVs(DataDirectory,fname_prefix)

    basin_keys = list(BasinInfoDF['basin_key'])
    basin_keys = [int(x) for x in basin_keys]

    basin_junctions = list(BasinInfoDF['outlet_junction'])
    basin_junctions = [int(x) for x in basin_junctions]

    print ('Basin keys are: ')
    print (basin_keys)

    labeldict = {}
    # get the best fit m/n for each basin
    if movern_method == "Chi_full":
        print ('Making a plot for the full chi method')
        Outlier_counter, removed_sources_dict, best_fit_movern_dict, MLEs_dict = CheckMLEOutliers(DataDirectory, fname_prefix, basin_list=basin_keys, start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, parallel=parallel)
        #MOverNDict = SimpleMaxMLECheck(BasinDF)
        m_over_ns = [round(i[0],2) for i in best_fit_movern_dict.values()]
        print ('You tested these m_over_ns')
        print (m_over_ns)
        MOverNDict = dict(zip(basin_keys,m_over_ns))
        label_list = [str(i)+": "+str(j) for i,j in zip(basin_keys,m_over_ns)]
        labeldict = dict(zip(basin_junctions,label_list))
        ImageName = raster_directory+fname_prefix+'_basins_movern_chi_full.'+FigFormat
        title = "Chi analysis (all data)"
    elif movern_method == "Chi_points":
        print ('Making a plot for the chi points method')
        if not parallel:
            PointsChiBasinDF = Helper.ReadMCPointsCSV(DataDirectory,fname_prefix)
        else:
            PointsChiBasinDF = Helper.AppendBasinPointCSVs(DataDirectory,fname_prefix)

        PointsDF = GetMOverNRangeMCPoints(PointsChiBasinDF,start_movern,d_movern,n_movern)
        moverns = PointsDF['Median_MOverNs'].tolist()
        MOverNDict = dict(zip(basin_keys,moverns))
        # labelling the basins
        moverns = [round(i,2) for i in moverns]
        label_list = [str(i)+": "+str(j) for i,j in zip(basin_keys,moverns)]
        labeldict = dict(zip(basin_junctions,label_list))
        ImageName = raster_directory+fname_prefix+'_basins_movern_chi_points.'+FigFormat
        title = "Chi bootstrap analysis"
    elif movern_method == "SA":
        SlopeAreaDF = SA.LinearRegressionRawData(DataDirectory,fname_prefix,parallel=parallel)
        moverns = SlopeAreaDF['regression_slope'].tolist()
        MOverNDict = dict(zip(basin_keys,moverns))
        # labelling the basins
        moverns = [round(i,2) for i in moverns]
        label_list = [str(i)+": "+str(j) for i,j in zip(basin_keys,moverns)]
        labeldict = dict(zip(basin_junctions,label_list))
        ImageName = raster_directory+fname_prefix+'_basins_movern_SA.'+FigFormat
        title = "Slope-area analysis"
    elif movern_method == "Chi_disorder":
        #DisorderDF = Helper.ReadDisorderCSV(DataDirectory, fname_prefix)
        DisorderDF = Helper.ReadDisorderUncertCSV(DataDirectory, fname_prefix)
        m_over_ns = DisorderDF['median'].values()
        label_list = [str(i)+": "+str(j) for i,j in zip(basin_keys,m_over_ns)]
        labeldict = dict(zip(basin_junctions,label_list))
        ImageName = raster_directory+fname_prefix+'_basins_movern_chi_disorder.'+FigFormat
        title = "Chi disorder analysis"

    else:
        print ("You didn't select an appropriate movern method. Please choose either 'Chi_full', 'Chi_points, or 'SA'.")
        sys.exit()

    # get moverns for cbar plotting. We always want a spacing of 0.1.
    min_max_str = ['min', 'max']
    end_movern = start_movern+d_movern*(n_movern-1)
    all_moverns = np.linspace(start_movern,end_movern,n_movern)
    print ("END MOVERN:")
    print (end_movern)
    min_max_moverns = [start_movern, end_movern]
    cbar_dict = dict(zip(min_max_str,min_max_moverns))

    # work out how many colours we need for the colormap
    n_colours = len(all_moverns)
    print ("N colours is: "+str(n_colours))

    # get a discrete colormap
    mn_cmap = plt.cm.Reds

    # going to make the basin plots - need to have bil extensions.
    print("I'm going to make the basin m/n plots. Your topographic data must be in ENVI bil format or I'll break!!")

    # get the rasters
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext
    BasinsName = fname_prefix+'_AllBasins'+raster_ext

    print ("THE TITLE IS:"+title)
    # create the map figure
    MF = MapFigure(HillshadeName, DataDirectory,coord_type="UTM_km", plot_title=title)
    # add the basins drape
    if lith:

        LithoMap = fname_prefix+"_LITHRAST"+raster_ext
        df_litho = pd.read_csv(DataDirectory+fname_prefix+"_lithokey.csv")
        MF.add_drape_image(LithoMap,DataDirectory,colourmap = plt.cm.jet,
                                alpha=0.4,
                                show_colourbar = False,
                                colorbarlabel = "Colourbar", discrete_cmap=True, n_colours=df_litho.shape[0],
                                norm = "None",
                                colour_min_max = [],
                                modify_raster_values=False,
                                old_values=[], new_values=[], cbar_type=int,
                                NFF_opti = False, custom_min_max = [])

    else:
        MF.add_basin_plot(BasinsName,fname_prefix,DataDirectory, value_dict = MOverNDict,
                          use_keys_not_junctions = True, show_colourbar = False,
                          discrete_cmap=True, n_colours=n_colours, colorbarlabel = "$m/n$",
                          colourmap = mn_cmap, adjust_text = False, cbar_dict=cbar_dict, parallel=parallel)

    # add the basin outlines
    if not parallel:
      Basins = LSDP.GetBasinOutlines(DataDirectory, BasinsName)
    else:
      Basins = LSDP.GetMultipleBasinOutlines(DataDirectory)

    MF.plot_polygon_outlines(Basins, linewidth=0.8)

    # add the channel network
    if not parallel:
        ChannelDF = Helper.ReadChiDataMapCSV(DataDirectory,fname_prefix)
    else:
        ChannelDF = Helper.AppendChiDataMapCSVs(DataDirectory,fname_prefix)
    ChannelPoints = LSDP.LSDMap_PointData(ChannelDF, data_type = "pandas", PANDEX = True)
    MF.add_point_data(ChannelPoints,show_colourbar="False", scale_points=True, column_for_scaling='drainage_area',alpha=0.1,zorder=100)

    # add the basin labelling
    if not parallel:
      Points = LSDP.GetPointWithinBasins(DataDirectory, BasinsName)
    else:
      Points = LSDP.GetPointsWithinMultipleBasins(DataDirectory, BasinsName)

    print ("Adding labels, the label dict is:", labeldict)
    MF.add_text_annotation_from_shapely_points_v2(Points, text_colour='k', label_dict=labeldict, zorder=200)

    if(lith):
        ImageName = ImageName[0:-(len(FigFormat)+1)]+"_lith."+FigFormat # adding a particle to the file
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 300) # Save the figure


#=============================================================================
# UNCERTAINTY
# Functions that make plots of uncertainties on the m/n analysis and comparison
# between the different methods
#=============================================================================

def PlotMOverNDicts(DataDirectory,fname_prefix,SA_based_dict,Chi_based_dict, FigFormat = "png", size_format = "ESURF"):
    """
    This plots the estimates of m over n for basins contained in the SA and Chi based dicts

    Args:
        SA_based_dict (dict): A dictionary where the key is the basin and the values are lists of the best fit m/n for S-A plots
        Chi_based_dict (dict):  A dictionary where the key is the basin and the values are lists of the best fit m/n after removing outliers

    Author: SMM
    """

    # We need to deal with the colours. Each colour will indicate a different
    # Data point


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

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.85,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])
    ax.grid(zorder=-100)

    # Loop through the basins adding data
    for basin in SA_based_dict:
        this_SA_list = SA_based_dict[basin]
        this_chi_list = Chi_based_dict[basin]

        this_SA_list = np.multiply(this_SA_list,-1)
        full_list = np.concatenate((this_SA_list,this_chi_list))
        basin_list = [basin]*len(full_list)

        colours = []
        this_num = -2
        for element in this_SA_list:
            colours.append(this_num)
            this_num = this_num+0.25

        this_num = 1
        for element in this_chi_list:
            colours.append(this_num)
            this_num = this_num+0.25


        #print("Hey there pal, let me print out the two lists for you")
        #print(this_SA_list)
        #print(this_chi_list)

        #print("These are the colour and basin lists")
        #print(full_list)
        #print(basin_list)


        #colours = range(len(full_list))


        ax.scatter(basin_list,full_list,c=colours,edgecolors='k',s=15,cmap = plt.cm.PuOr,zorder = 100)



    # set the axes labels
    ax.set_xlabel('basin')
    ax.set_ylabel('$m/n$')
    # Save the figure
    ImageName = DataDirectory+fname_prefix+'_MOverNbasins.'+FigFormat
    plt.savefig(ImageName, format=FigFormat, dpi=300)
    fig.clf()

def plot_MCMC_analysis(DataDirectory,fname_prefix,basin_list=[],FigFormat='png',size_format='ESURF', parallel=False):
    """
    This function makes a plot of the MCMC uncertainty analysis for the MLE collinearity
    NOTE: haven't set this up for basins in parallel yet (FJC 19/10/17)

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the DEM name without extension
        basin_list: list of basins to analyse. If none are passed then all basins are plotted.
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        parallel (bool): If true the data exists in multiple files for each basin and must be recombined (Default is False)

    Author: FJC
    """
    # check if a directory exists for the chi plots. If not then make it.
    MCMC_directory = DataDirectory+'MCMC_plots/'
    if not os.path.isdir(MCMC_directory):
        os.makedirs(MCMC_directory)

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

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.85,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])

    # get the basin info
    if not parallel:
      basin_df = Helper.ReadBasinInfoCSV(DataDirectory,fname_prefix)
      basin_stats_df = Helper.ReadBasinStatsCSV(DataDirectory,fname_prefix)
    else:
      basin_df = Helper.AppendBasinInfoCSVs(DataDirectory,fname_prefix)
      basin_stats_df = Helper.AppendBasinStatsCSVs(DataDirectory,fname_prefix)

    basin_keys = basin_df['basin_key']

    MOverNDict = SimpleMaxMLECheck(basin_stats_df)

    if basin_list == []:
        print ("You didn't give me a basin list so I'm going to plot all of them!")
        basin_list = basin_keys

    for basin in basin_list:
        # read in the chain csv file
        chain_df = Helper.ReadChainCSV(DataDirectory,fname_prefix,int(basin))

        best_fit_movern = MOverNDict[int(basin)]
        print (("The best fit m/n is:", best_fit_movern))
        n_iterations = len(chain_df['i'])
        #print n_iterations
        best_fit_line = [best_fit_movern] * len(chain_df['i'])

        # find the accepted values
        chain_df['AcceptedMask'] = chain_df['NAccepted'].diff()
        masked_df = chain_df.where(chain_df['AcceptedMask'] == 1.0)
        #print chain_df

        # plot the parameter with number of iterations
        ax.plot(masked_df['i'], masked_df['movern_New'], c='0.5', lw=0.5, zorder=1)
        ax.plot(chain_df['i'], best_fit_line, 'k--', zorder=100)

        #set plot labels
        ax.set_xlabel('N iterations')
        ax.set_ylabel('$m/n$')
        ax.set_title('Basin '+str(basin))

        # save the figure
        ImageName = MCMC_directory+fname_prefix+'_MCMC_basin' +str(basin)+'.'+FigFormat
        plt.savefig(ImageName, format=FigFormat, dpi=300)
        ax.cla()

def PlotMCPointsUncertainty(DataDirectory,fname_prefix, basin_list=[0], FigFormat='png',size_format='ESURF', start_movern=0.2,d_movern=0.1,n_movern=7,parallel=False):
    """
    This function makes a plot showing the range in m/n calculated
    using the MC points analysis. Makes a separate plot for each basin.

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the DEM name without extension
        basin_list: list of basins to be analysed. If empty then will analyse all of them.
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).

    Returns:
        plot of uncertainty in m/n

    Author:
        FJC
    """
    import matplotlib.patches as patches

    # check if a directory exists for the chi plots. If not then make it.
    points_directory = DataDirectory+'MC_points_plots/'
    if not os.path.isdir(points_directory):
        os.makedirs(points_directory)

    # get the basin info
    if not parallel:
        basin_df = Helper.ReadBasinInfoCSV(DataDirectory,fname_prefix)
    else:
        basin_df = Helper.AppendBasinInfoCSVs(DataDirectory,fname_prefix)

    basin_keys = basin_df['basin_key']

    if basin_list == []:
        print ("You didn't give me a basin list so I'm going to plot all of them!")
        basin_list = basin_keys

    # get the uncertainty df
    if not parallel:
        PointsChiBasinDF = Helper.ReadMCPointsCSV(DataDirectory,fname_prefix)
    else:
        PointsChiBasinDF = Helper.AppendBasinPointCSVs(DataDirectory,fname_prefix)

    PointsChiBasinDF = PointsChiBasinDF[PointsChiBasinDF['basin_key'].isin(basin_list)]

    UncertaintyDF = GetMOverNRangeMCPoints(PointsChiBasinDF,start_movern,d_movern,n_movern)

    # set up the plot
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

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.85,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])

    if not parallel:
        BasinStatsDF = Helper.ReadBasinStatsCSV(DataDirectory,fname_prefix)
    else:
        BasinStatsDF = Helper.AppendBasinCSVs(DataDirectory,fname_prefix)
    # best fit moverns
    best_fit_moverns = SimpleMaxMLECheck(BasinStatsDF)

    for basin_key in basin_keys:
        ThisBasinDF = PointsChiBasinDF[PointsChiBasinDF['basin_key'] == basin_key]

        # get the m/ns tested
        end_movern = start_movern+d_movern*(n_movern-1)
        all_moverns = np.linspace(start_movern,end_movern,n_movern)

        # get the median m/ns for this basin
        Medians = (ThisBasinDF.filter(regex='median')).values.tolist()[-1]
        FirstQs = (ThisBasinDF.filter(regex='FQ')).values.tolist()[-1]
        ThirdQs = (ThisBasinDF.filter(regex='TQ')).values.tolist()[-1]

        # now get the threshold value
        ThisUncertaintyDF = UncertaintyDF[UncertaintyDF['basin_key'] == basin_key]
        print (ThisUncertaintyDF)

        #plot the median and quartiles
        ax.plot(all_moverns,Medians,c="k",lw=1,zorder=2)
        ax.plot(all_moverns,FirstQs,c="k", lw=0.5, ls="--",zorder=2)
        ax.plot(all_moverns,ThirdQs,c="k", lw=0.5, ls="--",zorder=2)
        ax.fill_between(all_moverns, FirstQs, ThirdQs, color="0.8",alpha=0.75,zorder=1)

        #add a line for the threshold
        threshold = ThisUncertaintyDF.iloc[0]['FirstQ_threshold']
        min_movern = ThisUncertaintyDF.iloc[0]['Min_MOverNs']
        max_movern = ThisUncertaintyDF.iloc[0]['Max_MOverNs']
        threshold_moverns = np.linspace(min_movern,max_movern, 4)
        threshold_MLEs = np.full((n_movern,),threshold)
        ax.plot(all_moverns,threshold_MLEs,c='r',zorder=3, ls="--",lw=1)

        # add shaded background over range of m/n values
        ax.axvspan(min_movern,max_movern, alpha=0.1, color='red',zorder=0.2)

        # add arrow at best fit m/n
        # get the limits for the arrow
        max_MLE = max(Medians)
        min_MLE = min(FirstQs)
        dy = (max_MLE-min_MLE)/8
        spacing = 1.5
        # add arrow at best fit m/n
        ax.add_patch(
            patches.Arrow(
                best_fit_moverns[basin_key], #x
                max_MLE-(dy*spacing), #y
                0, #dx
                dy, #dy
                width = 0.05,
                facecolor = 'k',
                edgecolor = 'k')
            )
        ax.text(best_fit_moverns[basin_key]-0.075, max_MLE-1.5*(dy*spacing), "Best-fit "+r"$\theta$",fontsize=8)

        # set the axes labels
        ax.set_xlabel(r'$\theta$')
        ax.set_ylabel('MLE')
        #ax.set_title('Basin '+str(basin_key))
        ax.set_xlim(start_movern,end_movern)
        ax.set_ylim(min(FirstQs),max(ThirdQs)+0.0005)

        # save the figure
        ImageName = points_directory+fname_prefix+'_MC_points' +str(basin_key)+'.'+FigFormat
        plt.savefig(ImageName, format=FigFormat, dpi=300)
        ax.cla()

    fig.clf()
    plt.close(fig)

#=============================================================================
# SENSITIVITY FUNCTIONS
# Functions that make plots of sensitivity tests on the m/n analysis
#=============================================================================
def PlotSensitivityResultsSigma(DataDirectory,fname_prefix, FigFormat = "png", size_format = "ESURF", movern_method = "all"):
    """
    This function makes a plot of the results of a sensitivity analysis on sigma in the MLE method
    of calculating m/n.
    You need to specify the base directory - will look in every folder in this base directory
    with the prefix "Chi_analysis_sigma_" and get the m/n analysis in each sub-directory.
    NOTE - not set up for parallel (FJC 19/10/17)

    Args:
        DataDirectory (str): the root data directory
        fname_prefix (str): the DEM name without extension
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        movern_method (str): the movern method you want to test. Can be "all" or "points"

    Author: FJC
    """
    import os

    # Try doing this as a dataframe and masking?!?!
    columns = ["sigma", "basin_key", "m_over_n"]
    combined_DF = pd.DataFrame(columns=columns)

    # loop through each sub-directory with the sensitivity results
    MLE_str = "sigma_"
    for subdir, dirs, files in os.walk(DataDirectory):
        for dir in dirs:
            if MLE_str in dir:
                this_dir = DataDirectory+"/"+dir+'/'
                print (this_dir)

                this_df = pd.DataFrame(columns=columns)

                # get this value of sigma
                this_sigma = (dir.split("_"))[-1]

                # get the best fit m/n dataframe
                if movern_method == 'all':
                    BasinDF = Helper.ReadBasinStatsCSV(this_dir,fname_prefix)
                elif movern_method == 'points':
                    BasinDF = Helper.ReadBasinStatsPointCSV(this_dir,fname_prefix)
                else:
                    print ("You didn't specify a correct m/n method, you need to specify either 'all' or 'points'")
                MOverNDict = SimpleMaxMLECheck(BasinDF)
                this_df['sigma'] = [int(this_sigma)] * len(MOverNDict)
                this_df['basin_key'] = MOverNDict.keys()
                this_df['m_over_n'] = MOverNDict.values()
                combined_DF = combined_DF.append(this_df, ignore_index=True)

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

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.85,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])
    #ax.grid(zorder=-100)

    # now get the basin keys
    basin_keys = combined_DF['basin_key'].unique()
    basin_keys = [int(x) for x in basin_keys]

    keys = []
    sigmas = []
    # loop through the basin keys and plot sigma and m/n for each basin
    for key in basin_keys:
        this_df = combined_DF.loc[combined_DF['basin_key'] == key]
        #sort the data for plotting
        this_df.sort_values('sigma',inplace=True)
        this_df = this_df.loc[combined_DF['m_over_n'] > 0.2]
        if not this_df.empty:
            keys.append(int(key))
            sigmas.append(this_df['sigma'].iloc[0])

    ax.scatter(keys, sigmas, c='0.75', edgecolor='k')
    #
    # set the axes labels
    ax.set_xlabel('Basin key')
    ax.set_ylabel('$\sigma$ where $m/n$ is invariant')
    ax.set_xticks(basin_keys)

    # Save the figure
    ImageName = DataDirectory+fname_prefix+'_sensitivity.'+FigFormat
    plt.savefig(ImageName, format=FigFormat, dpi=300)
    fig.clf()

# def PlotSensitivityResultsMLEWithSigma(DataDirectory,fname_prefix, FigFormat = "png", size_format = "ESURF"):
#     """
#     This function makes a plot of the results of a sensitivity analysis on sigma in the MLE method
#     of calculating m/n.
#     You need to specify the base directory - will look in every folder in this base directory
#     with the prefix "sigma_" and get the m/n analysis in each sub-directory.
#     Makes a plot of MLE with increasing sigma to show where it becomes invariant
#
#     Args:
#         DataDirectory (str): the root data directory
#         fname_prefix (str): the DEM name without extension
#         FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
#         size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
#
#     Author: FJC
#     """
#     import os
#
#     # Try doing this as a dataframe and masking?!?!
#     columns = ["sigma", "basin_key", "m_over_n"]
#     combined_DF = pd.DataFrame(columns=columns)
#
#     # loop through each sub-directory with the sensitivity results
#     MLE_str = "sigma_"
#     for subdir, dirs, files in os.walk(DataDirectory):
#         for dir in dirs:
#             if MLE_str in dir:
#                 this_dir = DataDirectory+"/"+dir+'/'
#                 print this_dir
#
#                 this_df = pd.DataFrame(columns=columns)
#
#                 # get this value of sigma
#                 this_sigma = (dir.split("_"))[-1]
#
#                 # get the best fit m/n dataframe
#                 BasinDF = Helper.ReadBasinStatsCSV(this_dir,fname_prefix)
#                 MOverNDict = SimpleMaxMLECheck(BasinDF)
#                 this_df['sigma'] = [int(this_sigma)] * len(MOverNDict)
#                 this_df['basin_key'] = MOverNDict.keys()
#                 this_df['m_over_n'] = MOverNDict.values()
#                 combined_DF = combined_DF.append(this_df, ignore_index=True)
#
#     # Set up fonts for plots
#     label_size = 10
#     rcParams['font.family'] = 'sans-serif'
#     rcParams['font.sans-serif'] = ['arial']
#     rcParams['font.size'] = label_size
#
#     # make a figure
#     if size_format == "geomorphology":
#         fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
#         #l_pad = -40
#     elif size_format == "big":
#         fig = plt.figure(1, facecolor='white',figsize=(16,9))
#         #l_pad = -50
#     else:
#         fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
#         #l_pad = -35
#
#     gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.85,top=0.9)
#     ax = fig.add_subplot(gs[5:100,10:95])
#     #ax.grid(zorder=-100)
#
#     # now get the basin keys
#     basin_keys = combined_DF['basin_key'].unique()
#     basin_keys = [int(x) for x in basin_keys]
#
#     keys = []
#     sigmas = []
#     MLEs = []
#     # loop through the basin keys and plot sigma and m/n for each basin
#     for key in basin_keys:
#         this_df = combined_DF.loc[combined_DF['basin_key'] == key]
#         #sort the data for plotting
#         this_df.sort_values('sigma',inplace=True)
#         this_df = this_df.loc[combined_DF['m_over_n'] > 0.2]
#         if not this_df.empty:
#             keys.append(int(key))
#             sigmas.append(this_df['sigma'].iloc[0])
#
#     ax.scatter(keys, sigmas, c='0.75', edgecolor='k')
#     #
#     # set the axes labels
#     ax.set_xlabel('Basin key')
#     ax.set_ylabel('$\sigma$ where $m/n$ is invariant')
#     ax.set_xticks(basin_keys)
#
#     # Save the figure
#     ImageName = DataDirectory+fname_prefix+'_sensitivity.'+FigFormat
#     plt.savefig(ImageName, format=FigFormat, dpi=300)
#     fig.clf()
