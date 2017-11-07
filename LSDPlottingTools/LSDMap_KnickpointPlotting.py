## LSDMap_KnickpointPlotting.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## These functions are tools for analysing and plotting knickpoint data
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## FJC
## 05/06/2017
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import LSDPlottingTools.LSDMap_PointTools as LSDMap_PD
import matplotlib.pyplot as plt
import time as clock
from matplotlib import rcParams
import matplotlib.cm as cm
from LSDMapFigure.PlottingRaster import MapFigure
from LSDMapFigure.PlottingRaster import BaseRaster
from LSDMapFigure import PlottingHelpers as Helper
from LSDPlottingTools import colours as lsdcolours
from LSDPlottingTools import statsutilities as SUT
from LSDPlottingTools import init_plotting_DV
import LSDPlottingTools as LSDP
import sys
import os
import pandas as pd
from scipy.stats import norm


############################## Data (pre-)processing treatment ##########################################

def get_outliers_from_DF(df, method = ""):
    """
    proxy function to detect the outliers according to a specified method
    
    param:
        df (pandas Dataframe): the dataframe originally read from the _KsnKs.csv, potentially preprocessed
        method (str): method of outlier detection. Basin, source, general TOCOMPLETE AS WE ADD METHOD

    return:
        Dataframe containing the outliers
    
    Author: BG - 05/10/2017
    """
    if( method not in ["basin","river","general"]):
        print("ERROR: The method you are trying to use is not a valid one. \n I am aborting.")
        quit()

    if(method == "river"):
        print("I gonna detect your outliers using the river method: I'll try to detect the outliers in comparison of the river")
        df = SUT.extract_outliers_by_header(df,data_column_name = "rad_diff", header_for_group = "source_key", threshold = 2.5)

    elif(method == "basin"):
        print("I gonna detect your outliers using the river method: I'll try to detect the outliers in comparison of the basins")
        df = SUT.extract_outliers_by_header(df,data_column_name = "rad_diff", header_for_group = "basin_key", threshold = 2.5)

    elif(method == "basin"):
        print("I gonna detect your outliers using the river method: I'll try to detect the outliers in comparison of the basins")
        df["general"] = pd.Serie(np.ones(df.shape[0]),index = df.index)
        df = SUT.extract_outliers_by_header(df,data_column_name = "rad_diff", header_for_group = "general", threshold = 2.5)

    return df










###########################################################################################################

#################### Plotting function ####################################################################



def map_knickpoint_standard(DataDirectory, fname_prefix, size_format='ESURF', FigFormat='png', basin_list = [], mancut = 0, outlier_detection_method = "None"):
    
    """
    This creates a basic knickpoint map

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_list (list): List of the basin ID you want.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        mancut (float): manual cutoff for the data selection.
        outlier_detection_method (str): determine the outlier detection method to select the right knickpoints

    Returns:
        Shaded relief plot with the basins outlines and the knickpoints sized by intensity

    Author: BG - 05/10/2017
    """


    # check if a directory exists for the chi plots. If not then make it.
    raster_directory = DataDirectory+'raster_plots/'
    if not os.path.isdir(raster_directory):
        os.makedirs(raster_directory)

    # Set up fonts for plots
    basls = basin_list
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





    # get the rasters
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext
    BasinsName = fname_prefix+'_AllBasins'+raster_ext

    
    # create the map figure
    MF = MapFigure(HillshadeName, DataDirectory,coord_type="UTM_km", alpha = 0.7)
    
    # plot the basin outlines
    Basins = LSDP.GetBasinOutlines(DataDirectory, BasinsName)
    MF.plot_polygon_outlines(Basins, linewidth=0.5)

    # add the channel network without color
    ChannelDF = Helper.ReadMChiSegCSV(DataDirectory,fname_prefix)
    ChannelPoints = LSDP.LSDMap_PointData(ChannelDF, data_type = "pandas", PANDEX = True)
    MF.add_point_data(ChannelPoints,show_colourbar="False", scale_points=True,scaled_data_in_log= True, column_for_scaling='drainage_area',alpha=0.1,max_point_size = 0.5,min_point_size = 0.1,zorder=100)

    # add the knickpoints plots
    
    Kdf = Helper.ReadKnickpointCSV(DataDirectory,fname_prefix)
    # Selecting the basins in case you choose specific ones
    if(len(basls)>0):
        Kdf = Kdf[kdf["basin_key"].isin(basls)]

    # Sorting data in case of manual cut_off
    if(mancut>0):
        print("I am manually cutting off the data. If you need a automated outliers detection, switch mancut option off")
        Kdf = Kdf[Kdf["rad_diff"]> mancut]
    elif(outlier_detection_method != "None"):
        print("I will now select the outliers following the method " + outlier_detection_method)
        Kdf = get_outliers_from_DF(Kdf, method = outlier_detection_method)



    KdfPoints = LSDP.LSDMap_PointData(Kdf, data_type = "pandas", PANDEX = True)
    MF.add_point_data(KdfPoints,this_colourmap = "RdBu_r",column_for_plotting = "sign",show_colourbar="False", scale_points=True, scaled_data_in_log= False, column_for_scaling='rad_diff',alpha=1,max_point_size = 15,min_point_size = 1,zorder=200)

    #Saving and stuffs
    if(outlier_detection_method == "None"):
        outlier_detection_method = "raw"  
    ImageName = raster_directory+fname_prefix+"_knickpoint_"+ outlier_detection_method +"_map."+FigFormat
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 300) # Save the figure






def basic_hist(DataDirectory, fname_prefix ,basin_list = [] , size_format="ESURF", FigFormat=".png"):
    """
    Plot a simple histogram of the knickpoint repartition
    """
    print(" \n ########################## \n I am now going to print a basic histogram of your knickpoint in the area requested \n  ")

    from matplotlib.ticker import MaxNLocator

    # check if a directory exists for the sumarry plots. If not then make it.
    raster_directory = DataDirectory+'summary_plots/'
    if not os.path.isdir(raster_directory):
        os.makedirs(raster_directory)

    # loading the file:
    df = Helper.ReadKnickpointCSV(DataDirectory, fname_prefix)

    # Selecting the basin
    if(len(basin_list)>0):
        print("Selecting the basins %s" %(basin_list))
        df = df[df["basin_key"].isin(basin_list)]
    else:
        print("I am plotting for all the basins")

    # creating the figure

    plt.clf()
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))


    # Creating the axis
    
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.15,right=0.95,top=0.95)
    ax = fig.add_subplot(gs[0:100,0:100])

    # Plotting
    print("plotting ...")
    
    ls_baboty = []
    for i in df["basin_key"].unique():
        ls_baboty.append(df["rad_diff"][df["basin_key"] == i])

    n,bins, patch = ax.hist(ls_baboty,bins = 50, stacked  = True)
    n = np.array(n)

    # setting the yticks to be int

    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    

    #setting the labels
    ax.set_xlabel("Delta ksn")
    ax.set_ylabel("count")


    #Saving and stuffs   
    ImageName = raster_directory+fname_prefix+"_hist."+FigFormat
    plt.savefig(ImageName, dpi = 300) # Save the figure
    print(" done with saving your figure " + ImageName)
    plt.clf()





def chi_profile_knickpoint(DataDirectory, fname_prefix, size_format='ESURF', FigFormat='png', basin_list = [], mancut = 0, outlier_detection_method = "None", grouping = "basin_key", segments = True):
    
    """
    This creates a chi profiles with the knickpoint on top of the profile

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_list (list): List of the basin ID you want.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        mancut (float): manual cutoff for the data selection.
        outlier_detection_method (str): determine the outlier detection method to select the right knickpoints
        segments (bool): if segments is True, it plots the Mchi segmented elevation

    Returns:
        Shaded relief plot with the basins outlines and the knickpoints sized by intensity

    Author: BG - 05/10/2017
    """

    
    # check if a directory exists for the chi plots. If not then make it.
    raster_directory = DataDirectory+'chi_profile_knickpoint/'
    if(grouping == "source_key"):
        raster_directory = DataDirectory+'chi_profile_knickpoint_by_source/'
    if not os.path.isdir(raster_directory):
        os.makedirs(raster_directory)


    # Set up fonts for plots
    basls = basin_list
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    

    # add the channel network without color
    ChannelDF = Helper.ReadMChiSegCSV(DataDirectory,fname_prefix, type = "knickpoint")
    # add the knickpoints plots
    
    Kdf = Helper.ReadKnickpointCSV(DataDirectory,fname_prefix)

    segmented_elev = Helper.ReadMChiSegCSV(DataDirectory,fname_prefix, type = "knickpoint")
    # Selecting the basins in case you choose specific ones
    if(len(basls)>0):
        Kdf = Kdf[Kdf["basin_key"].isin(basls)]
        segmented_elev = segmented_elev[segmented_elev["basin_key"].isin(basls)]

    # Sorting data in case of manual cut_off
    if(mancut>0):
        print("I am manually cutting off the data. If you need a automated outliers detection, switch mancut option off")
        Kdf = Kdf[Kdf["rad_diff"]> mancut]
    elif(outlier_detection_method != "None"):
        print("I will now select the outliers following the method " + outlier_detection_method)
        Kdf = get_outliers_from_DF(Kdf, method = outlier_detection_method)


    #now plotting
    for hussard in Kdf[grouping].unique():
        print(hussard)
        # make a figure
        if size_format == "geomorphology":
            fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
            
        elif size_format == "big":
            fig = plt.figure(1, facecolor='white',figsize=(16,9))
            
        else:
            fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))

        # create the axis
        gs = plt.GridSpec(100,100,bottom=0.15,left=0.15,right=0.95,top=0.95)
        ax = fig.add_subplot(gs[0:100,0:100])

        # Selecting the data for this basin
        tcdf = ChannelDF[ChannelDF[grouping] == hussard] 
        tKdf = Kdf[Kdf[grouping] == hussard]
        tsegmented_elev = segmented_elev[segmented_elev[grouping] == hussard]
        # Setting the alpha, if segmented elevation is plotted, we want to see it under the chi plots
        if (segments):
            alo = 0.7
        else:
            alo = 1

        

        # Plotting the segmented elevation -  it can be computationally expensive depending on your number of segments
        if(segments):
            tcolseg = "#01DF01"
            for seg in tsegmented_elev["segment_number"].unique():
                
                # Setting the color of the segment and inverting it each turn
                if(tcolseg == "#01DF01"):
                    tcolseg = "#2E64FE"
                else:
                    tcolseg = "#01DF01"
                # selecting the unique segment I want
                teploseg = tsegmented_elev[tsegmented_elev["segment_number"] == seg]
                #plotting the segments
                ax.plot(teploseg["chi"],teploseg["segmented_elevation"], color = tcolseg, lw = 0.68)
        ## end with the segments

        # Plotting the chi profile
        ax.scatter(tcdf["chi"],tcdf["elevation"], c = tcdf["source_key"], cmap = "Accent", s = 1, lw = 0, alpha = alo)
        # Plotting the knickpoints
        sizel = abs(tKdf["rad_diff"]) * 200
        ax.scatter(tKdf["chi"],tKdf["elevation"], c = tKdf["sign"], s = sizel, alpha = 0.6, lw = 0.5,  cmap = "gnuplot", edgecolor = "k")

        # Details
        ax.set_xlabel("Chi")
        ax.set_ylabel("z")
        # saving details
        save_name = raster_directory + fname_prefix + "_" + grouping + str(hussard) + "_"+outlier_detection_method+"."+FigFormat

        plt.savefig(save_name, dpi = 600)
        plt.clf()


    print("done")


def chi_profile_knickzone(DataDirectory, fname_prefix, size_format='ESURF', FigFormat='png', basin_list = []):
    
    """
    This creates a chi profiles with the knickpoint on top of the profile

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_list (list): List of the basin ID you want.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        mancut (float): manual cutoff for the data selection.
        outlier_detection_method (str): determine the outlier detection method to select the right knickpoints
        segments (bool): if segments is True, it plots the Mchi segmented elevation

    Returns:
        Shaded relief plot with the basins outlines and the knickpoints sized by intensity

    Author: BG - 05/10/2017
    """

    
    # check if a directory exists for the chi plots. If not then make it.
    raster_directory = DataDirectory+'knickzone_river/'
    if not os.path.isdir(raster_directory):
        os.makedirs(raster_directory)


    # Set up fonts for plots
    basls = basin_list
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # read the knickpoint informations
    Kdf = Helper.ReadKnickpointCSV(DataDirectory,fname_prefix)
    Cdf = Helper.ReadMChiSegCSV(DataDirectory, fname_prefix, type = "knickpoint")

    # Selecting the basins in case you choose specific ones
    if(len(basls)>0):
        Kdf = Kdf[Kdf["basin_key"].isin(basls)]
        Cdf = Cdf[Cdf["basin_key"].isin(basls)]

    #now plotting
    print("I am plotting one figure per river, it can take a while. If you are processing a large area, I would recommend to select main channels")
    #for hussard in Kdf["source_key"].unique():#  TO KEEP!!!!!! TESTING ONE RIVER ATM
    for hussard in [0,19]:
        # make a figure
        if size_format == "geomorphology":
            fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
            
        elif size_format == "big":
            fig = plt.figure(1, facecolor='white',figsize=(16,9))
            
        else:
            fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))

        # create the axis
        gs = plt.GridSpec(100,100,bottom=0.15,left=0.15,right=0.85,top=0.95)
        ax = fig.add_subplot(gs[0:100,0:100])

        # Selecting the data for this river
        tKdf = Kdf[Kdf["source_key"] == hussard]
        tCdf = Cdf[Cdf["source_key"] == hussard]
        #Sorting by Chi values, not automatic since I am probably weridly using itrator to print the map in c++
        tKdf = tKdf.sort_values("chi")
        tCdf = tCdf.sort_values("chi")

        # Plotting the cumul ksn_variation
        ## shifting first an initial value at 0 for the variations
        tKdf.iloc[0, tKdf.columns.get_loc('chi')] = tCdf["chi"].min()
        ## then plotting
        ax.plot(tKdf["chi"],tKdf["cumul_ksn"], lw = 0.75, c = '#787878')
        ax.fill_between(tKdf["chi"],0,tKdf["cumul_ksn"], color = "k", alpha = 0.3)
        

        gs = plt.GridSpec(100,100,bottom=0.15,left=0.15,right=0.85,top=0.95)
        ax2 = fig.add_subplot(gs[0:100,0:100])
        ax2.patch.set_visible(False)
        ax.yaxis.set_ticks_position('right')
        ax2.yaxis.set_label_position('right')
        ax2.xaxis.set_ticks_position('none')
        ax2.yaxis.label.set_color('#787878') 
        ax2.plot(tCdf["chi"],tCdf["segmented_elevation"], lw = 1.2 , c ='#0089B9',zorder = 5)

        # plotting the knickpoint and setting a min-max size
        # min_size = 2
        # max_size = 700
        # sizepoint = (((tKdf["diff"]-tKdf["diff"].min())/tKdf["diff"].max())*max_size)+min_size
        #sizepoint = tKdf["diff"].copy()
        #sizepoint[sizepoint<20] = 20
        ax2.scatter(tKdf["chi"],tKdf["elevation"], c = tKdf["sign"],cmap = 'RdBu', s = tKdf["diff"].abs(), alpha = 0.7, lw = 0.5, edgecolor = "k", zorder = 10)

        ax.set_xlim(ax2.get_xlim())
        # Plotting the derivative
        #ax2.plot(tKdf["chi"],tKdf["deriv_cumul_ksn"].abs(), lw = 0.5, c = '#000000')

        #deriv_cumul_ksn

        # Details
        ax.set_xlabel(r'$\chi$')
        ax.set_ylabel(r'$k_(sn)$')
        ax2.set_ylabel(r'$\sum \Delta k_(sn)$')
        # saving details
        save_name = raster_directory + fname_prefix + "_KZ_Source" + str(hussard) + "."+FigFormat

        plt.savefig(save_name, dpi = 400)
        plt.clf()


    print("done")























#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~ THE FOLLOWING FUNCTIONS ARE TESTS AND EARLIER WORK, NOT USED WITH THE PlotKnickpointAnalysis ~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################



def plot_knickpoint_elevations(PointData, DataDirectory, basin_key=0, kp_threshold=0,
                               FigFileName='Image.pdf', FigFormat='pdf', size_format='ESURF', kp_type = "rad_diff"):
    """
    Function to create a plot of knickpoint elevation vs flow distance for each
    basin. Knickpoints are colour-coded by source node, and the marker size represents
    the magnitude of the knickpoint.

    Args:
        PointData: the LSDMap_PointData object with the knickpoint information
        DataDirectory (str): the data directory for the knickpoint file
        csv_name (str): name of the csv file with the knickpoint information
        basin_key (int): key to select the basin of interest
        kp_threshold (int): threshold knickpoint magnitude, any knickpoint below this will be removed (This option may be removed soon)
        kp_type (string): switch between diff and ratio data
        FigFileName (str): The name of the figure file
        FigFormat (str): format of output figure, can be 'pdf' (default), 'png', 'return', or 'show'
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).

    Returns:
        Plot of knickpoint elevations against flow distance

    Author: FJC
    """
    #PointData = LSDMap_PD.LSDMap_PointData(kp_csv_fname)
    # thin out small knickpoints
    KPData = PointData
    #KPData.ThinData(kp_type,kp_threshold)

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

    # get the data

    # Soert the dataset to the basin key
    KPData.selectValue('basin_key', value = basin_key, operator = "==")

    elevation = KPData.QueryData('elevation')

    flow_distance = KPData.QueryData('flow distance')

    magnitude = KPData.QueryData(kp_type)
    print("For the plotting, if you want to manage the scale, " +kp_type + " max is "+ str(magnitude.max()) +" and min is " + str(magnitude.min()))

    source = KPData.QueryData('source_key')




    #colour by source - this is the same as the script to colour channels over a raster,
    # (BasicChannelPlotGridPlotCategories) so that the colour scheme should match
    # make a color map of fixed colors
    NUM_COLORS = len(np.unique(source))
    this_cmap = plt.cm.Set1
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    plt.cm.ScalarMappable(norm=cNorm, cmap=this_cmap)
    channel_data = [x % NUM_COLORS for x in source]

    # normalise magnitude for plotting
    min_size = np.min(magnitude)
    max_size = np.max(magnitude)
    #normSize = [100*((x - min_size)/(max_size - min_size)) for x in magnitude]
    normSize = 100 * (magnitude - min_size)/(max_size - min_size)

    # now get the channel profiles that correspond to each knickpoint source and basin
    # PointData.ThinDataFromKey('basin_key',basin_key)
    # PointData.ThinDataSelection('source_key',maskSource)
    #
    # channel_elev = PointData.QueryData('elevation')
    # channel_elev = [float(x) for x in channel_elev]
    # channel_dist = PointData.QueryData('flow_distance')
    # channel_dist = [float(x) for x in channel_dist]

    # now plot the knickpoint elevations and flow distances
    #ax.scatter(channel_dist, channel_elev, s=0.1, c='k')
    ax.scatter(flow_distance, elevation, c = channel_data, cmap=this_cmap, s = normSize, lw=0.5,edgecolors='k',zorder=100)
    ax.set_xlabel('Flow distance (m)')
    ax.set_ylabel('Elevation (m)')

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


def plot_diff_ratio(PointData, DataDirectory, saveName = "Basic_diff_ratio", save_fmt = ".png", size_format = "ESURF", log_data = False):
    """
    Basic plot to have a general view of the knickpoints: diff against ratio colored by elevation

    Args:
        PointData: A PointData object
        DataDirectory: Where the data is saved
        saveName: save name

    returns:
        Nothing, sorry.
    Author: BG
    """
    plt.clf()
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

    diff = PointData.QueryData("diff")
    ratio = PointData.QueryData("ratio")

    if(log_data):
        print("I am logging the data")
        diff = np.log10(diff)
        ratio = np.log10(ratio)

    elevation =PointData.QueryData("elevation")
    ax.scatter(diff,ratio, s=0.5, lw = 0, c = elevation)
    ax.set_xlabel("Diff")
    ax.set_ylabel("Ratio")

    plt.savefig(DataDirectory+saveName+save_fmt,dpi=500)

def plot_basic_DA(PointData, DataDirectory, saveName = "Basic_DA", save_fmt = ".png", size_format = "ESURF", log_data = False):
    """
    Basic plot to have a general view of the knickpoints: drainage area against ratio and diff colored by elevation

    Args:
        PointData: A PointData object
        DataDirectory: Where the data is saved
        saveName: save name

    returns:
        Nothing, sorry.
    Author: BG
    """
    plt.clf()
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(2, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(2, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(2, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    diff = PointData.QueryData("diff")
    ratio = PointData.QueryData("ratio")
    DA = PointData.QueryData("drainage area")

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax1 = fig.add_subplot(gs[10:50,10:95])
    ax2 = fig.add_subplot(gs[50:100,10:95])



    if(log_data):
        print("I am logging the data")
        diff = np.log10(diff)
        ratio = np.log10(ratio)

    elevation = PointData.QueryData("elevation")
    DA = np.log10(DA)
    ax1.scatter(DA,ratio, s=0.7, lw = 0, c = elevation)
    ax1.set_ylabel("Ratio")
    ax1.tick_params(axis = 'x', length = 0, width = 0, labelsize = 0)
    ax1.spines['bottom'].set_visible(False)
    ax2.scatter(DA,diff,s=0.7, lw = 0, c = elevation)
    ax2.set_ylabel("Diff")
    ax2.set_xlabel("Drainage area")
    #ax2.set_xticks([1,2,3,4,5,6,7])
    ax2.tick_params(axis = 'x', labelsize = 6)
    ax1.set_xticks([4,5,6,7,8,9,10])
    ax2.set_xticks([4,5,6,7,8,9,10])
    ax2.set_xticklabels([ur"$10^{4}$",ur"$10^{5}$",ur"$10^{6}$",ur"$10^{7}$",ur"$10^{8}$",ur"$10^{9}$",ur"$10^{10}$"])

    plt.savefig(DataDirectory+saveName+save_fmt,dpi=500)

def plot_basic_FD(PointData, DataDirectory, saveName = "Basic_FD", save_fmt = ".png", size_format = "ESURF", log_data = False):
    """
    Basic plot to have a general view of the knickpoints: flow distance against ratio and diff colored by elevation

    Args:
        PointData: A PointData object
        DataDirectory: Where the data is saved
        saveName: save name

    returns:
        Nothing, sorry.
    Author: BG
    """
    plt.clf()
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(2, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(2, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(2, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    diff = PointData.QueryData("diff")
    ratio = PointData.QueryData("ratio")
    FD = PointData.QueryData("flow distance")

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax1 = fig.add_subplot(gs[10:50,10:95])
    ax2 = fig.add_subplot(gs[50:100,10:95])



    if(log_data):
        print("I am logging the data")
        diff = np.log10(diff)
        ratio = np.log10(ratio)

    elevation = PointData.QueryData("elevation")

    ax1.scatter(FD,ratio, s=0.7, lw = 0, c = elevation)
    ax1.set_ylabel("Ratio")
    ax1.tick_params(axis = 'x', length = 0, width = 0, labelsize = 0)
    ax1.spines['bottom'].set_visible(False)
    ax2.scatter(FD,diff,s=0.7, lw = 0, c = elevation)
    ax2.set_ylabel("Diff")
    ax2.set_xlabel("Flow distance")

    #ax2.tick_params(axis = 'x', labelsize = 6)
    #ax1.set_xticks([4,5,6,7,8,9,10])
    #ax2.set_xticks([4,5,6,7,8,9,10])
    #ax2.set_xticklabels([ur"$10^{4}$",ur"$10^{5}$",ur"$10^{6}$",ur"$10^{7}$",ur"$10^{8}$",ur"$10^{9}$",ur"$10^{10}$"])

    plt.savefig(DataDirectory+saveName+save_fmt,dpi=500)

def plot_basic_Z(PointData, DataDirectory, saveName = "Basic_Z", save_fmt = ".png", size_format = "ESURF", log_data = False):
    """
    Basic plot to have a general view of the knickpoints: flow distance against ratio and diff colored by elevation

    Args:
        PointData: A PointData object
        DataDirectory: Where the data is saved
        saveName: save name

    returns:
        Nothing, sorry.
    Author: BG
    """
    plt.clf()
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(2, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(2, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(2, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    diff = PointData.QueryData("diff")
    ratio = PointData.QueryData("ratio")
    Z = PointData.QueryData("elevation")

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax1 = fig.add_subplot(gs[10:50,10:95])
    ax2 = fig.add_subplot(gs[50:100,10:95])



    if(log_data):
        print("I am logging the data")
        diff = np.log10(diff)
        ratio = np.log10(ratio)

    sign = PointData.QueryData("sign")

    ax1.scatter(Z,ratio, s=0.7, lw = 0, c = sign)
    ax1.set_ylabel("Ratio")
    ax1.tick_params(axis = 'x', length = 0, width = 0, labelsize = 0)
    ax1.spines['bottom'].set_visible(False)
    ax2.scatter(Z,diff,s=0.7, lw = 0, c = sign)
    ax2.set_ylabel("Diff")
    ax2.set_xlabel("elevation")

    #ax2.tick_params(axis = 'x', labelsize = 6)
    #ax1.set_xticks([4,5,6,7,8,9,10])
    #ax2.set_xticks([4,5,6,7,8,9,10])
    #ax2.set_xticklabels([ur"$10^{4}$",ur"$10^{5}$",ur"$10^{6}$",ur"$10^{7}$",ur"$10^{8}$",ur"$10^{9}$",ur"$10^{10}$"])

    plt.savefig(DataDirectory+saveName+save_fmt,dpi=500)

def plot_outliers_x_vs_diff_ratio(PointData, DataDirectory,x_col = "elevation", saveName = "Outliers", save_fmt = ".png", size_format = "ESURF", log_data = False, ylim_ratio = [], ylim_diff = []):
    """
    Basic plot to have a general view of the knickpoints: flow distance against ratio and diff colored by elevation

    Args:
        PointData: A PointData object
        DataDirectory: Where the data is saved
        saveName: save name

    returns:
        Nothing, sorry.
    Author: BG
    """
    # Merging the dictionnary
    if(isinstance(PointData, dict)):
        print("Your data is a dictionnary of dataframes, let me create a PointTool object that contains all of these.")
        PointData = pd.concat(PointData)
        PointData = LSDMap_PD.LSDMap_PointData(PointData,data_type ="pandas", PANDEX = True)


    plt.clf()
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(2, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(2, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(2, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax1 = fig.add_subplot(gs[5:45,10:95])
    ax2 = fig.add_subplot(gs[55:100,10:95])


    diff = PointData.QueryData("diff")
    ratio = PointData.QueryData("ratio")
    Z = PointData.QueryData(x_col)
    PointData.selectValue('diff_outlier',value = True, operator = "==")
    PointData.selectValue('ratio_outlier',value = True, operator = "==")

    diffout = PointData.QueryData("diff")
    ratioout = PointData.QueryData("ratio")
    Z_out = PointData.QueryData(x_col)
    if(log_data):
        print("I am logging the data")

        diff = np.log10(diff)
        ratio = np.log10(ratio)
        if(isinstance(diffout, list) == False):
            diffout = [diffout]
            ratioout = [ ratioout]

        for i in range(len(diffout)):

            diffout[i]= np.log10(diffout[i])

            ratioout[i]= np.log10(ratioout[i])

    sign = PointData.QueryData("sign")

    ax1.scatter(Z,ratio, s=1.5, lw = 0, c = "gray")
    ax1.scatter(Z_out,ratioout, s=1.5, lw = 0, c = sign)
    ax1.set_ylabel("Ratio")
    ax1.tick_params(axis = 'x', length = 0, width = 0, labelsize = 0)
    ax1.spines['bottom'].set_visible(False)
    ax2.scatter(Z,diff,s=1.5, lw = 0, c = "gray")
    ax2.scatter(Z_out,diffout,s=1.5, lw = 0, c = sign)
    ax2.set_ylabel("Diff")
    ax2.set_xlabel(x_col)
    if(ylim_diff != []):
        ax2.set_ylim(ylim_diff[0],ylim_diff[1])
    if ylim_ratio != []:
        ax.set_ylim(ylim_ratio[0],ylim_ratio[1])

    #ax2.tick_params(axis = 'x', labelsize = 6)
    #ax1.set_xticks([4,5,6,7,8,9,10])
    #ax2.set_xticks([4,5,6,7,8,9,10])
    #ax2.set_xticklabels([ur"$10^{4}$",ur"$10^{5}$",ur"$10^{6}$",ur"$10^{7}$",ur"$10^{8}$",ur"$10^{9}$",ur"$10^{10}$"])

    plt.savefig(DataDirectory+saveName+save_fmt,dpi=500)

def plot_outliers_vs_others(PointData, DataDirectory, saveName = "Basic_diff_ratio", save_fmt = ".png", size_format = "ESURF", log_data = False):
    """
    Basic plot to have a general view of the knickpoints: diff against ratio colored by elevation

    Args:
        PointData: A PointData object or a dictionnary of dataframe containing outliers and none outliers
        DataDirectory: Where the data is saved
        saveName: save name

    returns:
        Nothing, sorry.
    Author: BG
    """
    # Merging the dictionnary
    if(isinstance(PointData, dict)):
        print("Your data is a dictionnary of dataframes, let me create a PointTool object that contains all of these.")
        PointData = pd.concat(PointData)
        PointData = LSDMap_PD.LSDMap_PointData(PointData,data_type ="pandas", PANDEX = True)


    plt.clf()
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


    diff = PointData.QueryData("diff")
    ratio = PointData.QueryData("ratio")

    PointData.selectValue('diff_outlier',value = True, operator = "==")
    PointData.selectValue('ratio_outlier',value = True, operator = "==")

    diffout = PointData.QueryData("diff")
    ratioout = PointData.QueryData("ratio")
    elevation = PointData.QueryData("elevation")
    if(log_data):
        print("I am logging the data")

        diff = np.log10(diff)
        ratio = np.log10(ratio)
        if(isinstance(diffout, list) == False):
            diffout = [diffout]
            ratioout = [ ratioout]

        for i in range(len(diffout)):

            diffout[i]= np.log10(diffout[i])

            ratioout[i]= np.log10(ratioout[i])
    print("Now plotting the outliers vs the non-outliers")
    ax.scatter(diff,ratio, s=0.5, lw = 0, c = 'gray')
    ax.scatter(diffout,ratioout, s = 0.5,c = elevation,lw = 0)
    ax.set_xlabel("Diff")
    ax.set_ylabel("Ratio")

    plt.savefig(DataDirectory+saveName+save_fmt,dpi=500)
    print("Your figure is " +DataDirectory+saveName+save_fmt)





############################ Mapping scripts ######################################




def map_custom():
    """
    Testing function to plot custom maps before creating real function for mapping routines

    Args:
        Yes.
    returns:
        No.
    Author:
        BG
    """
    ###### Parameters ######
    Directory = "/home/s1675537/PhD/DataStoreBoris/GIS/Data/Carpathian/knickpoint/" # reading directory
    wDirectory = Directory # writing directory
    Base_file = "Buzau" # It will be the cabkground raster. Each other raster you want to drap on it will be cropped to its extents including nodata
    csv_file = Directory + "test_sign_m_chi.csv" # Name of your point file, add a similar line with different name if you have more than one point file
    DrapeRasterName = "Buzau_hs.bil" # if you want to drap a raster on your background one. Just add a similar line in case you want another raster to drap and so on
    wname = "sign_test" # name of your output file
    dpi = 500 # Quality of your output image, don't exceed 900
    fig_size_inches = 7 # Figure size in Inches

    ##### Now we can load and plot the data

    BackgroundRasterName = Base_file + ".bil" # Ignore this line
    thisPointData = LSDP.LSDMap_PointData(csv_file, PANDEX = True) # Load the point file #1, add a similar line with different name if you have more than one point file.

    plt.clf() # Ignore this line

    MF = MapFigure(BackgroundRasterName, Directory,coord_type="UTM_km") # load the background raster

    MF.add_drape_image(DrapeRasterName,Directory, # Calling the function will add a drapped raster on the top of the background one
                        colourmap = "gray", # colormap used for this raster, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                        alpha=0.5, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                        show_colourbar = False, # Well, this one is explicit I think
                        colorbarlabel = "Colourbar") # Name of your Colourbar, it might bug though



    MF.add_point_data( thisPointData, # this function plot the requested point file using the lat/long column in the csv file
                       column_for_plotting = "m_chi_sign",  # Column used to color the data
                       this_colourmap = "cubehelix", # Colormap used, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                       colorbarlabel = "Colourbar", # Label
                       scale_points = False, # All the point will have the same size if False
                       column_for_scaling = "None", # If scale point True, you can scale the size of your points using one of the columns
                       scaled_data_in_log = False, # If scale point True, you can log the scaling
                       max_point_size = 5, # max size if scale point True again
                       min_point_size = 0.5, # You should be able to guess that one now
                       coulor_log = False, # do you want a log scale for your colorbar ?
                       coulor_manual_scale = [], #Do you want to manually limit the scale of your colorbar? if not let is false
                       manual_size = 0.5, # If none of above is choosen but you want to put another value than 0.5 to scale your point
                       alpha = 1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                       minimum_log_scale_cut_off = -10) # you probably won't need this

    ImageName = wDirectory+str(int(clock.time()))+wname+".png" # Ignore this
    ax_style = "Normal" # Ignore this
    MF.save_fig(fig_width_inches = fig_size_inches, FigFileName = ImageName, axis_style = ax_style, Fig_dpi = dpi) # Save the figure

def map_knickpoint_sign(PointData, DataDirectory, Raster_base_name, HS_name = "none",Time_in_name = False, river_network = "none", saveName = "none", size = 2, outliers = 'none'):
    """
    Will create a map of the knickpoint simply colored by sign.

    Args:
        PointData (PointTools object)
        DataDirectory (str): directory where the data will be saved and loaded.
        Raster_base_name (str): Base name of your files without the .bil
        HS_name (str): name of your Hillshade file, by default baseName + _hs.bil like LSDTT create it
        Time_in_name (bool): Option to add timing info in the nae of the figure. Can be useful if you test loads of parameters and you want to be sure that your files names are different (but awful).
    returns:
        No, but creates a map named map_knickpoint_sign.png
    Author:
        BG
    """
    ###### Parameters ######
    if(isinstance(PointData, dict)):
        print("Your data is a dictionnary of dataframes, let me create a PointTool object that contains all of these.")
        PointData = pd.concat(PointData)
        PointData = LSDMap_PD.LSDMap_PointData(PointData,data_type ="pandas", PANDEX = False)
    if(outliers != 'none' ):
        PointData.selectValue(outliers, operator = "==", value = True)
        print PointData

    Directory = DataDirectory # reading directory
    wDirectory = Directory # writing directory
    Base_file = Raster_base_name # It will be the cabkground raster. Each other raster you want to drap on it will be cropped to its extents including nodata
    if(saveName == "none"):
        wname = "map_knickpoint_sign" # name of your output file
    else:
        wname = saveName
    dpi = 500 # Quality of your output image, don't exceed 900
    fig_size_inches = 7 # Figure size in Inches
    if(HS_name == "none"):
        HS_name = Raster_base_name+("_hs.bil")
    DrapeRasterName = HS_name # if you want to drap a raster on your background one. Just add a similar line in case you want another raster to drap and so on

    ##### Now we can load and plot the data

    BackgroundRasterName = Base_file + ".bil" # Ignore this line
    plt.clf() # Ignore this line

    MF = MapFigure(BackgroundRasterName, Directory,coord_type="UTM_km", NFF_opti = True) # load the background raster

    MF.add_drape_image(DrapeRasterName,Directory, # Calling the function will add a drapped raster on the top of the background one
                        colourmap = "gray", # colormap used for this raster, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                        alpha=0.5, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                        show_colourbar = False, # Well, this one is explicit I think
                        colorbarlabel = "Colourbar", # Name of your Colourbar, it might bug though
                        NFF_opti = True)


    if(isinstance(river_network,LSDP.LSDMap_PointData)):
        MF.add_point_data( river_network, # this function plot the requested point file using the lat/long column in the csv file
                           column_for_plotting = "none",  # Column used to color the data
                           this_colourmap = "cubehelix", # Colormap used, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                           colorbarlabel = "Colourbar", # Label
                           scale_points = False, # All the point will have the same size if False
                           column_for_scaling = "None", # If scale point True, you can scale the size of your points using one of the columns
                           scaled_data_in_log = False, # If scale point True, you can log the scaling
                           max_point_size = 5, # max size if scale point True again
                           min_point_size = 0.5, # You should be able to guess that one now
                           coulor_log = False, # do you want a log scale for your colorbar ?
                           coulor_manual_scale = [], #Do you want to manually limit the scale of your colorbar? if not let is false
                           manual_size = 0.5, # If none of above is choosen but you want to put another value than 0.5 to scale your point
                           alpha = 1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                           minimum_log_scale_cut_off = -10) # you probably won't need this

    MF.add_point_data( PointData, # this function plot the requested point file using the lat/long column in the csv file
                       column_for_plotting = "sign",  # Column used to color the data
                       this_colourmap = "cubehelix", # Colormap used, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                       colorbarlabel = "Colourbar", # Label
                       scale_points = False, # All the point will have the same size if False
                       column_for_scaling = "None", # If scale point True, you can scale the size of your points using one of the columns
                       scaled_data_in_log = False, # If scale point True, you can log the scaling
                       max_point_size = 5, # max size if scale point True again
                       min_point_size = 0.5, # You should be able to guess that one now
                       coulor_log = False, # do you want a log scale for your colorbar ?
                       coulor_manual_scale = [], #Do you want to manually limit the scale of your colorbar? if not let is false
                       manual_size = size, # If none of above is choosen but you want to put another value than 0.5 to scale your point
                       alpha = 1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                       minimum_log_scale_cut_off = -10 # you probably won't need this
                      )
    if(Time_in_name):
        ImageName = wDirectory+str(int(clock.time()))+wname+".png" # Ignore this
    else:
        ImageName = wDirectory+wname+".png" # Ignore this
    ax_style = "Normal" # Ignore this
    MF.save_fig(fig_width_inches = fig_size_inches, FigFileName = ImageName, axis_style = ax_style, Fig_dpi = dpi) # Save the figure

def map_knickpoint_diff_sized_colored_ratio(PointData, DataDirectory, Raster_base_name, HS_name = "none",Time_in_name = False, river_network = "none", saveName = "none", log = False):
    """
    Will create a map of the knickpoint simply colored by sign.

    Args:
        PointData (PointTools object)
        DataDirectory (str): directory where the data will be saved and loaded.
        Raster_base_name (str): Base name of your files without the .bil
        HS_name (str): name of your Hillshade file, by default baseName + _hs.bil like LSDTT create it
        Time_in_name (bool): Option to add timing info in the nae of the figure. Can be useful if you test loads of parameters and you want to be sure that your files names are different (but awful).
    returns:
        No, but creates a map named map_knickpoint_sign.png
    Author:
        BG
    """
    ###### Parameters ######
    Directory = DataDirectory # reading directory
    wDirectory = Directory # writing directory
    Base_file = Raster_base_name # It will be the cabkground raster. Each other raster you want to drap on it will be cropped to its extents including nodata
    if(saveName == "none"):
        wname = "map_knickpoint_diff_sized_colored_ratio" # name of your output file
    else:
        wname = saveName
    dpi = 500 # Quality of your output image, don't exceed 900
    fig_size_inches = 7 # Figure size in Inches
    if(HS_name == "none"):
        HS_name = Raster_base_name+("_hs.bil")
    DrapeRasterName = HS_name # if you want to drap a raster on your background one. Just add a similar line in case you want another raster to drap and so on

    ##### Now we can load and plot the data

    BackgroundRasterName = Base_file + ".bil" # Ignore this line
    plt.clf() # Ignore this line

    MF = MapFigure(BackgroundRasterName, Directory,coord_type="UTM_km") # load the background raster

    MF.add_drape_image(DrapeRasterName,Directory, # Calling the function will add a drapped raster on the top of the background one
                        colourmap = "gray", # colormap used for this raster, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                        alpha=0.5, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                        show_colourbar = True, # Well, this one is explicit I think
                        colorbarlabel = "Colourbar") # Name of your Colourbar, it might bug though



    if(isinstance(river_network,LSDP.LSDMap_PointData)):
        MF.add_point_data( river_network, # this function plot the requested point file using the lat/long column in the csv file
                           column_for_plotting = "none",  # Column used to color the data
                           this_colourmap = "cubehelix", # Colormap used, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                           colorbarlabel = "Colourbar", # Label
                           scale_points = False, # All the point will have the same size if False
                           column_for_scaling = "None", # If scale point True, you can scale the size of your points using one of the columns
                           scaled_data_in_log = False, # If scale point True, you can log the scaling
                           max_point_size = 5, # max size if scale point True again
                           min_point_size = 0.5, # You should be able to guess that one now
                           coulor_log = False, # do you want a log scale for your colorbar ?
                           coulor_manual_scale = [], #Do you want to manually limit the scale of your colorbar? if not let is false
                           manual_size = 0.1, # If none of above is choosen but you want to put another value than 0.5 to scale your point
                           alpha = 1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                           minimum_log_scale_cut_off = -10) # you probably won't need this

    MF.add_point_data( PointData, # this function plot the requested point file using the lat/long column in the csv file
                       column_for_plotting = "ratio",  # Column used to color the data
                       this_colourmap = "cubehelix", # Colormap used, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                       colorbarlabel = "Ratio", # Label
                       scale_points = True, # All the point will have the same size if False
                       column_for_scaling = "diff", # If scale point True, you can scale the size of your points using one of the columns
                       scaled_data_in_log = log, # If scale point True, you can log the scaling
                       max_point_size = 20, # max size if scale point True again
                       min_point_size = 0.5, # You should be able to guess that one now
                       coulor_log = log, # do you want a log scale for your colorbar ?
                       coulor_manual_scale = [], #Do you want to manually limit the scale of your colorbar? if not let is false
                       manual_size = 0.5, # If none of above is choosen but you want to put another value than 0.5 to scale your point
                       alpha = 1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                       minimum_log_scale_cut_off = -10) # you probably won't need this

    if(Time_in_name):
        ImageName = wDirectory+str(int(clock.time()))+wname+".png" # Ignore this
    else:
        ImageName = wDirectory+wname+".png" # Ignore this
    ax_style = "Normal" # Ignore this
    MF.save_fig(fig_width_inches = fig_size_inches, FigFileName = ImageName, axis_style = ax_style, Fig_dpi = dpi) # Save the figure






def DEPRECATED_map_knickpoint_standard(PointData, DataDirectory, Raster_base_name, HS_name = "none",Time_in_name = False, river_network = "none", saveName = "none", log = False):
    """
    Will create a map of the knickpoint sized by their delta and colored by size.

    Args:
        PointData (PointTools object)
        DataDirectory (str): directory where the data will be saved and loaded.
        Raster_base_name (str): Base name of your files without the .bil
        HS_name (str): name of your Hillshade file, by default baseName + _hs.bil like LSDTT create it
        Time_in_name (bool): Option to add timing info in the nae of the figure. Can be useful if you test loads of parameters and you want to be sure that your files names are different (but awful).
    returns:
        No, but creates a map named map_knickpoint_sign.png
    Author:
        BG
    """
    ###### Parameters ######
    Directory = DataDirectory # reading directory
    wDirectory = Directory # writing directory
    Base_file = Raster_base_name # It will be the cabkground raster. Each other raster you want to drap on it will be cropped to its extents including nodata
    if(saveName == "none"):
        wname = "map_knickpoint_std" # name of your output file
    else:
        wname = saveName
    dpi = 500 # Quality of your output image, don't exceed 900
    fig_size_inches = 7 # Figure size in Inches
    if(HS_name == "none"):
        HS_name = Raster_base_name+("_hs.bil")
    DrapeRasterName = HS_name # if you want to drap a raster on your background one. Just add a similar line in case you want another raster to drap and so on

    ##### Now we can load and plot the data

    BackgroundRasterName = Base_file + ".bil" # Ignore this line
    plt.clf() # Ignore this line

    MF = MapFigure(BackgroundRasterName, Directory,coord_type="UTM_km") # load the background raster

    MF.add_drape_image(DrapeRasterName,Directory, # Calling the function will add a drapped raster on the top of the background one
                        colourmap = "gray", # colormap used for this raster, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                        alpha=0.5, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                        show_colourbar = True, # Well, this one is explicit I think
                        colorbarlabel = "Colourbar") # Name of your Colourbar, it might bug though



    if(isinstance(river_network,LSDP.LSDMap_PointData)):
        MF.add_point_data( river_network, # this function plot the requested point file using the lat/long column in the csv file
                           column_for_plotting = "none",  # Column used to color the data
                           this_colourmap = "cubehelix", # Colormap used, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                           colorbarlabel = "Colourbar", # Label
                           scale_points = False, # All the point will have the same size if False
                           column_for_scaling = "None", # If scale point True, you can scale the size of your points using one of the columns
                           scaled_data_in_log = False, # If scale point True, you can log the scaling
                           max_point_size = 5, # max size if scale point True again
                           min_point_size = 0.5, # You should be able to guess that one now
                           coulor_log = False, # do you want a log scale for your colorbar ?
                           coulor_manual_scale = [], #Do you want to manually limit the scale of your colorbar? if not let is false
                           manual_size = 0.1, # If none of above is choosen but you want to put another value than 0.5 to scale your point
                           alpha = 1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                           minimum_log_scale_cut_off = -10) # you probably won't need this

    MF.add_point_data( PointData, # this function plot the requested point file using the lat/long column in the csv file
                       column_for_plotting = "sign",  # Column used to color the data
                       this_colourmap = "cubehelix", # Colormap used, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                       colorbarlabel = "diff", # Label
                       scale_points = True, # All the point will have the same size if False
                       column_for_scaling = "diff", # If scale point True, you can scale the size of your points using one of the columns
                       scaled_data_in_log = log, # If scale point True, you can log the scaling
                       max_point_size = 20, # max size if scale point True again
                       min_point_size = 5, # You should be able to guess that one now
                       coulor_log = log, # do you want a log scale for your colorbar ?
                       coulor_manual_scale = [], #Do you want to manually limit the scale of your colorbar? if not let is false
                       manual_size = 0.5, # If none of above is choosen but you want to put another value than 0.5 to scale your point
                       alpha = 1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                       minimum_log_scale_cut_off = -10) # you probably won't need this

    if(Time_in_name):
        ImageName = wDirectory+str(int(clock.time()))+wname+".png" # Ignore this
    else:
        ImageName = wDirectory+wname+".png" # Ignore this
    ax_style = "Normal" # Ignore this
    MF.save_fig(fig_width_inches = fig_size_inches, FigFileName = ImageName, axis_style = ax_style, Fig_dpi = dpi) # Save the figure


def plot_pdf_diff_ratio(df, DataDirectory, saveName = "pdf_diff_ratio", save_fmt = ".png", size_format = "ESURF",  xlim =[]):
    """
    Basic plot to have a general view of the knickpoints: flow distance against ratio and diff colored by elevation

    Args:
        PointData: A PointData object
        DataDirectory: Where the data is saved
        saveName: save name

    returns:
        Nothing, sorry.
    Author: BG
    """
    plt.clf()
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(2, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(2, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(2, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35


    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax1 = fig.add_subplot(gs[10:50,10:95])
    ax2 = fig.add_subplot(gs[50:100,10:95])


    ax1.scatter(df["ratio"],norm.pdf(df["ratio"]),lw =0, s = 1, c = "red")
    ax1.set_ylabel("Ratio")
    ax1.tick_params(axis = 'x', length = 0, width = 0, labelsize = 0)
    ax1.spines['bottom'].set_visible(False)
    ax2.scatter(df["diff"],norm.pdf(df["diff"]),lw =0, s = 1, c = "red")
    ax2.set_ylabel("Diff")
    ax2.set_xlabel("PDF")


    #'###### Setting the limits
    if(xlim != []):
        ax2.set_xlim(xlim[0],xlim[1])
        ax1.set_xlim(xlim[0],xlim[1])



    #ax2.tick_params(axis = 'x', labelsize = 6)
    #ax1.set_xticks([4,5,6,7,8,9,10])
    #ax2.set_xticks([4,5,6,7,8,9,10])
    #ax2.set_xticklabels([ur"$10^{4}$",ur"$10^{5}$",ur"$10^{6}$",ur"$10^{7}$",ur"$10^{8}$",ur"$10^{9}$",ur"$10^{10}$"])

    plt.savefig(DataDirectory+saveName+save_fmt,dpi=500)

def violin_by_bin(ldf, DataDirectory, saveName = "Violin", column = "elevation", size_format = "ESURF"):

    """
    Will plot violin from a list of bins. NOT READY YET.

    Author: BG

    matplotlib description:
        Violin plots are similar to histograms and box plots in that they show
    an abstract representation of the probability distribution of the
    sample. Rather than showing counts of data points that fall into bins
    or order statistics, violin plots use kernel density estimation (KDE) to
    compute an empirical distribution of the sample. That computation
    is controlled by several parameters. This example demonstrates how to
    modify the number of points at which the KDE is evaluated (``points``)
    and how to modify the band-width of the KDE (``bw_method``).

    For more information on violin plots and KDE, the scikit-learn docs
    have a great section: http://scikit-learn.org/stable/modules/density.html
    """

    plt.clf()
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(2, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(2, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(2, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35


    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax1 = fig.add_subplot(gs[10:50,10:95])
    ax2 = fig.add_subplot(gs[50:100,10:95])


    ax2.set_ylabel("Ratio")
    ax1.set_ylabel("Diff")
    ax2.set_xlabel(column)
    plt.savefig(DataDirectory+saveName+save_fmt,dpi=500)


def pdf_from_bin(ldf, DataDirectory, saveName = "BasicPDF_", column = "elevation", size_format = "ESURF" ):

    """
    Produce some simple pdf plots from a list of pandas dataframe.

    Arg:

    Returns: nothing, but produce a plot.

    Author: BG
    """

    for inch in ldf:
        plt.clf()
        label_size = 10
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['arial']
        rcParams['font.size'] = label_size

        # make a figure
        if size_format == "geomorphology":
            fig = plt.figure(2, facecolor='white',figsize=(6.25,3.5))
            l_pad = -40
        elif size_format == "big":
            fig = plt.figure(2, facecolor='white',figsize=(16,9))
            l_pad = -50
        else:
            fig = plt.figure(2, facecolor='white',figsize=(4.92126,3.5))
            l_pad = -35



        gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
        ax1 = fig.add_subplot(gs[10:50,10:95])
        ax2 = fig.add_subplot(gs[50:100,10:95])

        ax2.scatter(ldf[inch]["diff"],norm.pdf(ldf[inch]["diff"]), s = 1.5, lw = 0)
        ax1.scatter(ldf[inch]["ratio"],norm.pdf(ldf[inch]["ratio"]), s = 1.5, lw = 0)

        ax2.set_ylabel("PDF (Diff)")
        ax1.set_ylabel("PDF (Ratio)")
        ax2.set_xlabel("Diff/ratio binned by " + column + "_" + inch)
        plt.savefig(DataDirectory+saveName+inch+"_"+column+".png",dpi=500)


def pdf_from_bin_one_col(ldf, DataDirectory, saveName = "BasicPDF_", column = "elevation", size_format = "ESURF", pdf_col = "diff", combine_diff_sign = False, argsort = False ):

    """
    Produce some simple pdf plots from a dict of pandas dataframe.

    Arg:

    Returns: nothing, but produce a plot.

    Author: BG
    """




    for inch in ldf:
        plt.clf()
        label_size = 10
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['arial']
        rcParams['font.size'] = label_size



        if(combine_diff_sign):
            ldf[inch]["diff"][ldf[inch]["sign"] == -1] = -ldf[inch]["diff"][ldf[inch]["sign"] == -1]

        data = np.array(ldf[inch][pdf_col].values)


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
        ax1 = fig.add_subplot(gs[10:100,5:95])

        print(data.shape)
        if(data.shape[0]>0):
            ax1.hist(data, 100, normed=1, facecolor='green', alpha=0.75)






        ax1.set_ylabel("PDF")
        ax1.set_xlabel("elevation by " + pdf_col)
        ax1.set_xlim(-100,100)
        plt.savefig(DataDirectory+saveName+inch+"_"+column+".png",dpi=500)

def plot_2d_density_map(dataframe, DataDirectory, columns = ["drainage area", "diff"], bin = 50,   saveName = "BasicPDF_", size_format = "ESURF",):

    """
    Plots a 2d histogram or density plot or heatmap depending how you name it of two variables.

    Args:
        dataframe: a Pandas dataframe
        columns (list of str): The x,y columns to plot
        bin (int): number of bins

    returns:
        Nothing yet, plot a figure.
    """

    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    X_data = dataframe[columns[0]]
    Y_data = dataframe[columns[1]]









#
