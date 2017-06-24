#=============================================================================
# This script creates figures for visualising the m/n selection using the chi
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
import numpy as np
import LSDPlottingTools as LSDP
import matplotlib.pyplot as plt
#import seaborn as sns
from matplotlib import rcParams
import pandas as pd
from matplotlib import colors
from shapely.geometry import Polygon
from LSDMapFigure import PlottingHelpers as Helper
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

    # zip into a dictionary
    MOverNDict = dict(zip(basin_keys, MOverNs))
    return MOverNDict


def CheckMLEOutliers(DataDirectory, fname_prefix, basin_list=[0], start_movern=0.2, d_movern=0.1, n_movern=7):
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

    # we open the first file just so that we can get a counter list
    FirstDF = Helper.ReadFullStatsCSV(DataDirectory,fname_prefix,m_over_n_values[0])

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
    for m_over_n in m_over_n_values:

        #load the file
        FullStatsDF = Helper.ReadFullStatsCSV(DataDirectory,fname_prefix,m_over_n)

        print("This_m_over_n is: "+str(m_over_n))

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
                                                                                n_movern)
        best_fit_movern_dict[basin_number] = movern_of_max_MLE
        removed_sources_dict[basin_number] = remove_list_index
        MLEs_dict[basin_number] = MLEs

    print("Here are the vitalstatisix, chief: ")
    print(best_fit_movern_dict)

    return Outlier_counter, removed_sources_dict, best_fit_movern_dict, MLEs_dict

def Iteratively_recalculate_MLE_removing_outliers_for_basin(Outlier_counter, DataDirectory, fname_prefix, basin_number, start_movern=0.2, d_movern=0.1, n_movern=7):
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

    print("I am going to recalcualte MLE for basin: "+str(basin_number))
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
                                                                             remove_list_index)

    # Returns the remove_list_index, which is a list where each element
    # is a list of tributaries removed in an iteration,
    # And the maximum MLE values in each iteration.
    return remove_list_index,movern_of_max_MLE, MLEs



def Calculate_movern_after_iteratively_removing_outliers(movern_list, DataDirectory,
                                                         fname_prefix, basin_number,
                                                         remove_list_index):
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
        MLE_vals = RecalculateTotalMLEWithRemoveList(DataDirectory, fname_prefix,
                                                     m_over_n,basin_number, remove_list_index)
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
                                      movern,basin_number, remove_list_index):
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
    FullStatsDF = Helper.ReadFullStatsCSV(DataDirectory,fname_prefix,movern)

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
                  start_movern = 0.2, d_movern = 0.1, n_movern = 7):
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

    profile_suffix = "_movern.csv"
    basin_stats_suffix = "_movernstats_basinstats.csv"

    movern_profile_file = fname_prefix+profile_suffix
    movern_basin_stats_file = fname_prefix+basin_stats_suffix

    end_movern = start_movern+d_movern*(n_movern-1)
    m_over_n_values = np.linspace(start_movern,end_movern,n_movern)

    # get the maximum MLE of each basin
    pd_DF = Helper.ReadBasinStatsCSV(DataDirectory, fname_prefix)
    shp = pd_DF.shape
    max_MLEs = []
    max_MLEs_index = []
    for i in range(0,shp[0]):
        #print("I is: "+str(i))
        a = pd_DF.loc[[i]]
        b = np.asarray(a)
        c = b[0,:]
        max_MLEs.append(max(c))
        max_MLEs_index.append(np.argmax(c))
    print("max_MLEs are: ")
    print(max_MLEs)
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
    thisPointData = LSDP.LSDMap_PointData(DataDirectory+movern_profile_file)
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
            #plt.title("Basin = " +mn_legend+", MLE = "+short_MLE)

            #newline = "\n"
            title_string = "Basin "+str(basin_key)+", $m/n$ = "+str(mn)
            title_string2 = "MLE = "+short_MLE
            ax.text(0.05, 0.95, title_string,
                    verticalalignment='top', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='black', fontsize=10)
            print("The basin index is: "+str(basin_key)+" and the max index is: "+str(max_MLEs_index[basin_key]))
            if( idx == max_MLEs_index[basin_key]):
                print("This m/n is: "+str(mn)+" and it is the maximum MLE")
                ax.text(0.05, 0.88, title_string2+", maximum MLE in basin.",
                    verticalalignment='top', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='red', fontsize=10)
            else:
                ax.text(0.05, 0.88, title_string2,
                    verticalalignment='top', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='black', fontsize=10)

            #save the plot
            newFilename = DataDirectory+"Chi_profiles_basin_"+str(basin_key)+"_"+counter+".png"

            # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
            ax.tick_params(axis='both', width=1, pad = 2)
            for tick in ax.xaxis.get_major_ticks():
                tick.set_pad(2)

            FigFormat = "png"
            plt.savefig(newFilename,format=FigFormat,dpi=300)
            ax.cla()

def MakeChiPlotsMLE(DataDirectory, fname_prefix, basin_list=[0], start_movern=0.2, d_movern=0.1, n_movern=7,
                    size_format='ESURF', FigFormat='png'):
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

    Returns:
        Plot of each m/n value for each basin.

    Author: FJC
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

    # read in the csv files
    ProfileDF = Helper.ReadChiProfileCSV(DataDirectory, fname_prefix)
    BasinStatsDF = Helper.ReadBasinStatsCSV(DataDirectory, fname_prefix)

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

    for m_over_n in m_over_n_values:
        # read in the full stats file
        print("This m/n is: "+str(m_over_n))
        FullStatsDF = Helper.ReadFullStatsCSV(DataDirectory,fname_prefix,m_over_n)

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

            # get the chi, elevation, and MLE for the tributaries
            TributariesX = list(ProfileDF_tribs[movern_key])
            TributariesElevation = list(ProfileDF_tribs['elevation'])
            TributariesMLE = list(ProfileDF_tribs['MLE'])

            # get the colourmap to colour channels by the MLE value
            #NUM_COLORS = len(MLE)
            MLE_array = np.asarray(TributariesMLE)
            this_cmap = plt.cm.Reds
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

            # label with the basin and m/n
            title_string = "Basin "+str(basin_key)+", $m/n$ = "+str(m_over_n)
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
            newFilename = DataDirectory+"MLE_profiles"+str(basin_key)+"_"+str(m_over_n)+".png"

            # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
            ax.tick_params(axis='both', width=1, pad = 2)
            for tick in ax.xaxis.get_major_ticks():
                tick.set_pad(2)

            plt.savefig(newFilename,format=FigFormat,dpi=300)
            ax.cla()
            ax2.cla()



def PlotProfilesRemovingOutliers(DataDirectory, fname_prefix, basin_list=[0], start_movern=0.2, d_movern=0.1, n_movern=7):
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
    FirstDF = Helper.ReadFullStatsCSV(DataDirectory,fname_prefix,start_movern)

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
    Outlier_counter, removed_sources_dict, best_fit_movern_dict, MLEs_dict = CheckMLEOutliers(DataDirectory, fname_prefix, basin_list, start_movern, d_movern, n_movern)

    # Now get the chi profiles of all the basins and channels
    # Load from file and put into a pandas data frame
    ProfileDF = Helper.ReadChiProfileCSV(DataDirectory,fname_prefix)

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
            FullStatsDF = Helper.ReadFullStatsCSV(DataDirectory,fname_prefix,best_fit_movern)

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
            title_string = "Basin "+str(basin_number)+", best fit $m/n$ = "+str(best_fit_movern)
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

def PlotMLEWithMOverN(DataDirectory, fname_prefix, basin_list = [0], size_format='ESURF', FigFormat='png', start_movern=0.2, d_movern = 0.1, n_movern = 7):
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
    FirstDF = Helper.ReadFullStatsCSV(DataDirectory,fname_prefix,start_movern)

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
    Outlier_counter, removed_sources_dict, best_fit_movern_dict, MLEs_dict = CheckMLEOutliers(DataDirectory, fname_prefix, basin_list, start_movern, d_movern, n_movern)

    # Now get the chi profiles of all the basins and channels
    # Load from file and put into a pandas data frame
    ProfileDF = Helper.ReadChiProfileCSV(DataDirectory,fname_prefix)

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
        print best_fit_moverns

        # loop through the number of removed tributaires and get the MLE for each m/n for each iteration
        for i in range(n_removed_sources+1):

            # get the MLE of the best fit m over n
            best_fit_movern = best_fit_moverns[i]
            # get the index in the MLE list
            idx = int(round((best_fit_movern - start_movern)/d_movern,0))
            print idx
            best_fit_MLE = basin_MLEs[idx][i]
            print ("The best fit MLE is: "+str(best_fit_MLE)+", where m/n = " +str(best_fit_movern))

            # get the MLEs for this iteration
            these_MLEs = basin_MLEs[:,i]

            # get the ratio of these MLEs to the best fit
            ratio_MLEs = [x/best_fit_MLE for x in these_MLEs]

            # no removed tributaries
            if i == 0:
                # plot the data
                ax.plot(m_over_n_values,ratio_MLEs, lw=1.5, label = str(i), c='k', linestyle = '-', zorder=100)

                # get the limits for the arrow
                max_MLE = max(ratio_MLEs)
                min_MLE = min(ratio_MLEs)
                dy = (max_MLE-min_MLE)/8
                spacing = 1.1
                # add arrow at best fit m/n
                ax.add_patch(
                    patches.Arrow(
                        best_fit_movern, #x
                        max_MLE+(dy*spacing), #y
                        0, #dx
                        -dy, #dy
                        width = 0.05,
                        facecolor = 'r',
                        edgecolor = 'r'
                    )
                )
            #remove tribs
            else:
                # plot the data
                ax.plot(m_over_n_values,ratio_MLEs, lw=1, label = str(i), c='0.5', linestyle = ls[i]) # different linestyle for each iteration?

        # set the axes labels
        ax.set_xlabel('$m/n$')
        ax.set_ylabel('$MLE$ ratio')

        # set the ylim
        ax.set_ylim(min_MLE,max_MLE+(dy*spacing))

        # add the legend
        ax.legend(loc='right', bbox_to_anchor=(1.25,0.5), title = 'Iterations', frameon=False)

        # some formatting of the figure
        ax.spines['top'].set_linewidth(1)
        ax.spines['left'].set_linewidth(1)
        ax.spines['right'].set_linewidth(1)
        ax.spines['bottom'].set_linewidth(1)

        # label with the basin and m/n
        best_fit_movern = best_fit_movern_dict[basin_number][0]
        title_string = "Basin "+str(basin_number)+"; Best fit $m/n$: "+str(best_fit_movern)
        #ax.set_title(title_string)
        ax.text(0, 1.1, title_string,
            verticalalignment='top', horizontalalignment='left',
            transform=ax.transAxes,
            color='red', fontsize=10)

        #save the plot
        newFilename = DataDirectory+"MLE_fxn_movern_"+str(basin_number)+"."+FigFormat

        # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
        ax.tick_params(axis='both', width=1, pad = 2)
        for tick in ax.xaxis.get_major_ticks():
            tick.set_pad(2)

        plt.savefig(newFilename,format=FigFormat,dpi=300)
        ax.cla()

#=============================================================================
# RASTER PLOTTING FUNCTIONS
# Functions that interface with LSDMapFigure to plot the m/n analysis with
# spatial data.
#=============================================================================
def MakeRasterPlotsBasins(DataDirectory, fname_prefix, size_format='ESURF', FigFormat='png'):
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
    #import modules
    from LSDMapFigure.PlottingRaster import MapFigure
    from LSDMapFigure.PlottingRaster import BaseRaster
    import LSDPlottingTools.LSDMap_BasicManipulation as LSDMap_BM
    import LSDPlottingTools.LSDMap_PointTools as LSDMap_PT

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
    BaseLevelDF = Helper.ReadBaselevelKeysCSV(DataDirectory, fname_prefix)

    baselevel_keys = list(BaseLevelDF['baselevel_key'])
    baselevel_keys = [float(x) for x in baselevel_keys]

    baselevel_junctions = list(BaseLevelDF['baselevel_junction'])
    baselevel_junctions = [float(x) for x in baselevel_junctions]

    print ('Basin keys are: ')
    print baselevel_keys

    # get a discrete colormap
    cmap = plt.cm.jet

    # going to make the basin plots - need to have bil extensions.
    print("I'm going to make the basin plots. Your topographic data must be in ENVI bil format or I'll break!!")

    # get the rasters
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext
    BasinsName = fname_prefix+'_AllBasins'+raster_ext
    print (BasinsName)

    # create the map figure
    MF = MapFigure(HillshadeName, DataDirectory,coord_type="UTM_km", colourbar_location='bottom')

    # add the basins drape
    MF.add_drape_image(BasinsName, DataDirectory, colourmap = cmap, alpha = 0.8, colorbarlabel='Basin ID', discrete_cmap=True, n_colours=len(baselevel_keys), show_colourbar = True, modify_raster_values=True, old_values=baselevel_junctions, new_values=baselevel_keys, cbar_type = int)

    # add the basin outlines
    Basins = LSDMap_BM.GetBasinOutlines(DataDirectory, fname_prefix)
    MF.plot_polygon_outlines(Basins, linewidth=0.8)

    # add the basin labelling
    # get basin keys as a point data object
    basin_points = LSDMap_PT.LSDMap_PointData(BaseLevelDF, data_type='pandas', PANDEX=True)
    MF.add_text_annotation_from_points(basin_points,column_for_plotting='baselevel_key', PANDEX=True, text_colour='k')

    # Save the figure
    ImageName = DataDirectory+fname_prefix+'_basin_keys.'+FigFormat
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 300)

def MakeRasterPlotsMOverN(DataDirectory, fname_prefix, n_movern=7, size_format='ESURF', FigFormat='png'):
    """
    This function makes a shaded relief plot of the DEM with the basins coloured
    by the best fit m/n

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        n_movern (int): The number of m/n values tested, default = 7.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.

    Returns:
        Shaded relief plot with the basins coloured by best fit m/n

    Author: FJC
    """
    #import modules
    from LSDMapFigure.PlottingRaster import MapFigure
    from LSDMapFigure.PlottingRaster import BaseRaster
    import LSDPlottingTools.LSDMap_BasicManipulation as LSDMap_BM

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
    BasinDF = Helper.ReadBasinStatsCSV(DataDirectory,fname_prefix)

    basin_junctions = list(BasinDF['outlet_jn'])
    basin_junctions = [float(x) for x in basin_junctions]

    # get the best fit m/n for each junction
    MOverNDict = SimpleMaxMLECheck(BasinDF)
    print MOverNDict

    m_over_ns = MOverNDict.values()

    # get a discrete colormap
    mn_cmap = plt.cm.Reds

    # going to make the basin plots - need to have bil extensions.
    print("I'm going to make the basin m/n plots. Your topographic data must be in ENVI bil format or I'll break!!")

    # get the rasters
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext
    BasinsName = fname_prefix+'_AllBasins'+raster_ext
    print (BasinsName)

    # create the map figure
    MF = MapFigure(HillshadeName, DataDirectory,coord_type="UTM_km", colourbar_location='bottom')
    # add the basins drape
    MF.add_drape_image(BasinsName, DataDirectory, colourmap = mn_cmap, alpha = 0.8, colorbarlabel='Best fit m/n', discrete_cmap=True, n_colours=n_movern, show_colourbar = True, modify_raster_values=True, old_values=basin_junctions, new_values=m_over_ns, cbar_type=float)
    Basins = LSDMap_BM.GetBasinOutlines(DataDirectory, fname_prefix)
    MF.plot_polygon_outlines(Basins, linewidth=0.8)

    ImageName = DataDirectory+fname_prefix+'_basins_movern.'+FigFormat
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 300) # Save the figure

if __name__ == "__main__":

    # Change these filenames and paths to suit your own files
    DataDirectory = '/home/s0923330/DEMs_for_analysis/kentucky_srtm/'
    #fname_prefix = 'Kentucky_DEM'
    #DataDirectory = 'T:\\analysis_for_papers\\movern_testing\\'
    #DataDirectory = 'C:\\VagrantBoxes\\LSDTopoTools\\Topographic_projects\\Irian_jaya\\'
    fname_prefix = 'Kentucky_DEM'

    # specify the figure size and format
    size_format='ESURF'
    FigFormat = 'png'

    # either specify a list of the basins, or set as empty to get all of them
    basin_list = []

    # specify the m/n values tested
    start_movern = 0.2
    d_movern = 0.1
    n_movern = 8

    # analyse the MLE
    #CheckMLEOutliers(DataDirectory, fname_prefix, basin_list, start_movern=0.2, d_movern=0.1, n_movern=7)
    #PlotProfilesRemovingOutliers(DataDirectory, fname_prefix, basin_list, start_movern, d_movern, n_movern)
    #PlotMLEWithMOverN(DataDirectory, fname_prefix, basin_list=basin_list, start_movern=start_movern, d_movern=d_movern, n_movern=n_movern)

    # run the plotting functions
    #MakePlotsWithMLEStats(DataDirectory, fname_prefix, basin_list, start_movern, d_movern, n_movern)
    # MakeChiPlotsMLE(DataDirectory, fname_prefix, basin_list, start_movern, d_movern, n_movern,
    #                   size_format=size_format, FigFormat=FigFormat)

    # run the raster plotting
    MakeRasterPlotsBasins(DataDirectory, fname_prefix, size_format, FigFormat)
    #MakeRasterPlotsMOverN(DataDirectory, fname_prefix, n_movern, size_format, FigFormat)
