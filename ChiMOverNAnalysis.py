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
import seaborn as sns
from matplotlib import rcParams
import pandas as pd
from matplotlib import colors
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
    pd_DF = pd.DataFrame.from_csv(DataDirectory+movern_basin_stats_file)
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
    profile_suffix = "_movern.csv"
    basin_stats_suffix = "_movernstats_basinstats.csv"
    movern_profile_file = fname_prefix+profile_suffix
    movern_basin_stats_file = fname_prefix+basin_stats_suffix

    # get the data frames
    ProfileDF = pd.read_csv(DataDirectory+movern_profile_file)
    BasinStatsDF = pd.read_csv(DataDirectory+movern_basin_stats_file)

    # get the number of basins
    basin_keys = list(BasinStatsDF['basin_key'])
    basin_keys = [float(x) for x in basin_keys]

    # get the list of basins
    if basin_list == []:
        print("You didn't give me a list of basins, so I'll just run the analysis on all of them!")
        basin_list = basin_keys

    # loop through each m over n value
    end_movern = start_movern+d_movern*(n_movern-1)
    m_over_n_values = np.linspace(start_movern,end_movern,n_movern)

    for m_over_n in m_over_n_values:

        print("This m/n is: "+str(m_over_n))

        # get the full stats data for this m over n
        fullstats_suffix = "_movernstats_%s_fullstats.csv" %(str(m_over_n))
        fullstats_profile_file = fname_prefix+fullstats_suffix
        print("The stats filename is "+ fullstats_profile_file)
        FullStatsDF = pd.read_csv(DataDirectory+fullstats_profile_file)

        # loop through all the basins in the basin list
        for basin_key in basin_list:
            print("This basin key is %s") %str(basin_key)

            # mask the data frames for this basin
            ProfileDF_basin = ProfileDF[ProfileDF['basin_key'] == basin_key]
            FullStatsDF_basin = FullStatsDF[FullStatsDF['basin_key'] == basin_key]

            # get the data frame for the main stem
            ProfileDF_MS = ProfileDF_basin[ProfileDF_basin['source_key'] == FullStatsDF_basin['reference_source_key'][0]]

            # get the data frame for the tributaries
            ProfileDF_basin = ProfileDF_basin[ProfileDF_basin['source_key'] != FullStatsDF_basin['reference_source_key'][0]]
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


def CheckMLEOutliers(DataDirectory, fname_prefix, basin_list=[0], start_movern=0.2, d_movern=0.1, n_movern=7):
    """
    This function uses the fullstas files to search for outliers in the 
    channels, and also give statistics of the distribution of MLE values. 

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



    # Get a vector of the m over n values
    end_movern = start_movern+d_movern*(n_movern-1)
    m_over_n_values = np.linspace(start_movern,end_movern,n_movern)
    
    
    # we open the first file just so that we can get a counter list
    full_filename = DataDirectory+fname_prefix+"_movernstats_"+str(m_over_n_values[0])+"_fullstats.csv"
    FirstDF = pd.read_csv(full_filename)

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
        
        print("The basin list is now: ")
        print(basin_list)
    
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
        
        
    
    # Now we loop through all the files, calulating the outliers
    for m_over_n in m_over_n_values:
        full_filename = DataDirectory+fname_prefix+"_movernstats_"+str(m_over_n)+"_fullstats.csv"
        #print("\n\n\nFilename is: ")
        #print(full_filename)
        
        #load the file
        FullStatsDF = pd.read_csv(full_filename)
        
        # loop through the basins
        for basin in basin_list:
            
            # mask the data so you only get the correct basin
            FullStatsDF_basin = FullStatsDF[FullStatsDF['basin_key'] == basin]

            # extract the relevant data
            MLE_values = list(FullStatsDF_basin['MLE'])
            RMSE_values = list(FullStatsDF_basin['RMSE'])
            trib_values = list(FullStatsDF_basin['test_source_key']) 
            #print("The MLE values are: ")
            #print(MLE_values)          
            
            # now get the outliers
            MLE_array = np.asarray(MLE_values)
            RMSE_array = np.asarray(RMSE_values)
            RMSE_outliers = LSDP.lsdstatsutilities.is_outlier(RMSE_array)
            MLE_outliers = LSDP.lsdstatsutilities.is_outlier(RMSE_array)
 
            #print("The RMSE outliers are:")
            #print(RMSE_outliers)            
            #import itertools
            #filtered = list(itertools.compress(RMSE_array, RMSE_outliers))
            #print("The filtered array is: ")
            #print(filtered)
           
            # now check each of the outlier arrays to see if e need to flip the array
            RMSE_index_min = np.argmin(RMSE_array)
            #RMSE_min = np.min(RMSE_array)
 
            # if the max MLE is an outlier, flip the outlier vector
            if (RMSE_outliers[RMSE_index_min]):
                RMSE_outliers = [not i for i in RMSE_outliers]
                
            #filtered = list(itertools.compress(RMSE_array, RMSE_outliers))
            #print("Now the filtered is: ")
            #print(filtered)              
            
            
            #print("\nThe MLE outliers are:")
            #print(MLE_outliers)            
            #filtered = list(itertools.compress(MLE_array, MLE_outliers))
            #print("The filtered array is: ")
            #print(filtered)
           
            MLE_index_max = np.argmax(MLE_array) 
            #MLE_max = np.argmax(MLE_array)
            
            # if the max MLE is an outlier, flip the outlier vector
            if (MLE_outliers[MLE_index_max]):
                MLE_outliers = [not i for i in MLE_outliers]
                
            #filtered = list(itertools.compress(MLE_array, MLE_outliers))
            #print("Now the filtered is: ")
            #print(filtered)            
            
            # turn the outliers vector into an integer
            int_Outlier = [int(i) for i in MLE_outliers]
            #print("Integer outliers are: ")
            #print(int_Outlier)
            
            Outlier_counter[basin] = Outlier_counter[basin]+int_Outlier
            
            

           
            # some formatting of the figure
            ax.spines['top'].set_linewidth(1)
            ax.spines['left'].set_linewidth(1)
            ax.spines['right'].set_linewidth(1)
            ax.spines['bottom'].set_linewidth(1)

            # make the lables
            ax.set_xlabel("$\chi$ (m)")
            ax.set_ylabel("Elevation (m)")
            
            # make a pretty seaborn plot
            sns.distplot(MLE_array, ax=ax, rug=True, hist=False)
            ax.plot(MLE_outliers, np.zeros_like(MLE_outliers), 'ro', clip_on=False)

            kwargs = dict(y=0.95, x=0.05, ha='left', va='top')
            ax.set_title('MAD-based Outliers', **kwargs)

            title_string = "Basin "+str(basin)+", $m/n$ = "+str(m_over_n)
            ax.text(0.05, 0.95, title_string,
                    verticalalignment='top', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='black', fontsize=10)            
        
            #save the plot
            newFilename = DataDirectory+"MLE_outliers"+str(basin)+"_"+str(m_over_n)+".png"

            # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
            ax.tick_params(axis='both', width=1, pad = 2)
            for tick in ax.xaxis.get_major_ticks():
                tick.set_pad(2)

            plt.savefig(newFilename,format=FigFormat,dpi=300)
            ax.cla()  
    
    # now show the outlier counter
    for basin in basin_list:
        print("The outlier counter in basin: "+str(basin)+" is: ")
        print(Outlier_counter[basin])
        
    

if __name__ == "__main__":

    # Change these filenames and paths to suit your own files
    #DataDirectory = '/home/s0923330/DEMs_for_analysis/kentucky_srtm/'
    #fname_prefix = 'Kentucky_chi'
    DataDirectory = 'T:\\analysis_for_papers\\movern_testing\\'
    #DataDirectory = 'C:\\VagrantBoxes\\LSDTopoTools\\Topographic_projects\\Irian_jaya\\'
    fname_prefix = 'Irian_Jaya_PP'    
    
    size_format='ESURF'
    FigFormat = 'png'

    # either specify a list of the basins, or set as empty to get all of them
    basin_list = []

    # specify the m/n values tested
    start_movern = 0.2
    d_movern = 0.1
    n_movern = 7
    
    CheckMLEOutliers(DataDirectory, fname_prefix, basin_list, start_movern=0.2, d_movern=0.1, n_movern=7)
    

    # run the plotting function
    #MakePlotsWithMLEStats(DataDirectory, fname_prefix, basin_list, start_movern, d_movern, n_movern)
    #MakeChiPlotsMLE(DataDirectory, fname_prefix, basin_list, start_movern, d_movern, n_movern,
    #                    size_format=size_format, FigFormat=FigFormat)
