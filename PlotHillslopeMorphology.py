#=============================================================================
# Script to plot various hillslope metrics that are generated from the
# <<not sure which on, maybe make_spaghetti.cpp>> driver file
#
# Authors:
#     Fiona J. Clubb
#     Martin D Hurst
#     Simon M. Mudd
#=============================================================================
#=============================================================================
# IMPORT MODULES
#=============================================================================
# set backend to run on server
import matplotlib
matplotlib.use('Agg')

#from __future__ import print_function
import sys
import os
from LSDPlottingTools import LSDMap_HillslopeMorphology as HS
from LSDMapFigure import PlottingHelpers as Helper


#=============================================================================
# This parses a list of lists separated string. Each list is separated by a colon
#=============================================================================
def parse_list_of_list_from_string(a_string):
    """
    This parses a list of lists separated string. Each list is separated by a colon

    Args:
        a_string (str): This creates a list of lists. Each sub list is separated by colons and the sub list items are separated by commas. So `1,2,3:4,5` would produce [ [1,2,3],[4,5]]

    Returns:
        list_of_list (list): A list of lists

    Author: SMM

    Date: 11/01/2018
    """

    if len(a_string) == 0:
        print("No list of list found. I will return an empty list.")
        list_of_list = []
    else:
        listified_entry = [item for item in a_string.split(':')]
        list_of_list = []

        # now loop through these creating a dict
        for entry in listified_entry:
            split_entry = [int(item) for item in entry.split(',')]
            list_of_list.append(split_entry)

    print("This list of lists is: ")
    print(list_of_list)

    return list_of_list


#=============================================================================
# This is just a welcome screen that is displayed if no arguments are provided.
#=============================================================================
def print_welcome():

    print("\n\n=======================================================================")
    print("Hello! I'm going to plot the hillslope-channel results for you.")
    print("You will need to tell me which directory to look in.")
    print("Use the -dir flag to define the working directory.")
    print("If you don't do this I will assume the data is in the same directory as this script.")
    print("For help type:")
    print("   python PlotHillslopeMorphologyAnalysis.py -h\n")
    print("=======================================================================\n\n ")

#=============================================================================
# This is the main function that runs the whole thing
#=============================================================================
def main(argv):

    # print("On some windows systems you need to set an environment variable GDAL_DATA")
    # print("If the code crashes here it means the environment variable is not set")
    # print("Let me check gdal enviroment for you. Currently is is:")
    # print(os.environ['GDAL_DATA'])
    #os.environ['GDAL_DATA'] = os.popen('gdal-config --datadir').read().rstrip()
    #print("Now I am going to get the updated version:")
    #print(os.environ['GDAL_DATA'])

    # If there are no arguments, send to the welcome screen
    if not len(sys.argv) > 1:
        full_paramfile = print_welcome()
        sys.exit()

    # Get the arguments
    import argparse
    parser = argparse.ArgumentParser()

    # The location of the data files
    parser.add_argument("-dir", "--base_directory", type=str, help="The base directory that contains your data files. If this isn't defined I'll assume it's the same as the current directory.")
    parser.add_argument("-fname", "--fname_prefix", type=str, help="The prefix of your DEM WITHOUT EXTENSION!!! This must be supplied or you will get an error (unless you're running the parallel plotting).")
    parser.add_argument("-plotdir","--plot_directory",type=str, help="The directory to place output plots in. If it isn't defined I'll create a new folder inside --base-directory called \"plots\" and put them there.")
    # The basins you want to plot
    parser.add_argument("-basin_keys", "--basin_keys",type=str,default = "", help = "This is a comma delimited string that gets the list of basins you want for the plotting. Default = no basins")
    parser.add_argument("-basin_lists", "--basin_lists",type=str,default = "", help = "This is a string that initiates a list of a list for grouping basins. The object becomes a list of a list but the syntax is comma seperated lists, and each one is separated by a colon. Default = no dict")
    parser.add_argument("-basin_labels", "--basin_labels",type=str,default = "", help = "This is a list of labels if you want to change the naming of the basins. Must be the same length as the basin keys")
    # The analysis you want to do
    parser.add_argument("-Mchi", "--plot_mchi", type=bool, default=False, help="If this is true, I'll make some plots of the hillslope-channel data against Mchi")
    parser.add_argument("-CHT", "--plot_CHT", type=bool, default=False, help="If this is true, I'll make some plots of hilltop curvature against data from the channel segments.")
    parser.add_argument("-segments", "--plot_segments", type=bool, default=False, help="If this is true, I'll make some long profile plots of the channel segments.")
    parser.add_argument("-in_basin", "--plot_data_within_basin", type=bool, default=False, help="If this is true, I'll make plots of the hillslope data vs distance upstream for each basin")
    parser.add_argument("-means", "--plot_mean_basin_data", type=bool, default=False, help="If this is true I'll make plots of the mean hillslope data vs basin ID")
    parser.add_argument("-means_uplift", "--plot_means_with_uplift", type=bool, default=False, help="If this is true I'l make plots of the mean hillslope data vs uplift rate")
    parser.add_argument("-traces", "--plot_hillslope_traces",type=bool, default=False, help="if this is true I'll plot a hillshade with hillslope traces overlain")
    parser.add_argument("-angles", "--plot_junction_angles", type=bool, default=False, help="If this is true I'll make plots of the basin junction angles vs basin ID")
    parser.add_argument("-determine_sc", "--determine_sc", type=bool, default=False, help="If this is true I will create some plots to determine what the correct critical slope value should be, based on Grieve et al. (2016, ESPL)")
    parser.add_argument("-profile_plots", "--profile_plots", type=bool, default=False, help="If this is true I will plot E*, R*, and either Ksn or elevation against chi for each basin.")
    parser.add_argument("-EsRs_basin", "--plot_Es_Rs_by_basin", type=bool, default=False, help="This plots E* R* data for each basin.")
    parser.add_argument("-Ksn_vs_EsRs_basin", "--plot_Ksn_vs_Es_Rs_by_basin", type=bool, default=False, help="This plots K_sn vs E* R* data for each basin.")
    parser.add_argument("-plot_stacked_Es_Rs_by_chi", "--plot_stacked_Es_Rs_by_chi", type=bool, default=False, help="This plots stacked E* or R* as function of Chi.")
    parser.add_argument("-plot_clustered_Es_Rs_by_chi", "--plot_clustered_Es_Rs_by_chi", type=bool, default=False, help="This plots clustered E* or R* as function of Chi.")
    parser.add_argument("-plot_clustered_Es_Rs_plus_theoretical", "--plot_clustered_Es_Rs_plus_theoretical", type=bool, default=False, help="This plots clustered E* or R* plots that have been clustered.")
    parser.add_argument("-plot_Lh_hist", "--plot_Lh_hist", type=bool, default=False, help="This plots the full distribution of hillslope lengths for specified basins")

    # Parameters that are used within plotting functions
    parser.add_argument("-sc", "--sc", type=float, default=0.8, help="Critical slope for E*R* calculations. Default = 0.8")
    parser.add_argument("-extent", "--custom_plot_extent", type=float, nargs=4, default=None, help="four values required to define the [xmin, xmax, ymin, ymax] extent to plot map data")
    parser.add_argument("-plot_Ksn", "--plot_Ksn", type=bool, default=False, help="If this is true I will plot Ksn in the profile plots.")
    parser.add_argument("-plot_chi_by_basin", "--plot_chi_by_basin", type=bool, default=False, help="If this is true I will plot with chi coordinate in main stem basin plots. If false I will plot by flow distance.")
    parser.add_argument("-min_traces", "--minimum_traces", type=int, default=50, help="The minimium number of traces to be used for plotting segment data.")
    parser.add_argument("-EsRs_colour_by", "--EsRs_colour_by", type=str, default="ksn", help="What to colour E* R* plots by.")
    parser.add_argument("-mainstem_only", "--mainstem_only", type=bool, default=False, help="If true, use only mainstem data for E* R* plots.")




    parser.add_argument("-common_max_Ksn", "--common_max_Ksn", type=float, default=-99, help="If positive, use as common maximum Ksn for all plots")
    parser.add_argument("-common_max_Es", "--common_max_Es", type=float, default=-99, help="If positive, use as common maximum E star for all plots")

    #parse the argments
    args = parser.parse_args()

    if not args.fname_prefix:
        print("WARNING! You haven't supplied your DEM name. Please specify this with the flag '-fname'")
        sys.exit()

    # get the base directory
    if args.base_directory:
        this_dir = args.base_directory
        # check if you remembered a / at the end of your path_name
        if not this_dir.endswith("/"):
            print("You forgot the '/' at the end of the directory, appending...")
            this_dir = this_dir+"/"
    else:
        this_dir = os.getcwd()


    # check the basins
    print("You told me that the basin keys are: ")
    print(args.basin_keys)

    if len(args.basin_keys) == 0:
        print("No basins found, I will plot all the basins")
        df = HS.ReadChannelData(this_dir, args.fname_prefix)
        these_basin_keys = list(df['basin_key'].unique())
        print(these_basin_keys)
    else:
        these_basin_keys = [int(item) for item in args.basin_keys.split(',')]
        print("The basins I will plot are:")
        print(these_basin_keys)

    if len(args.basin_labels) != 0:
        these_basin_labels = [item for item in args.basin_labels.split(',')]
        print("I will use these basin labels:")
        print(these_basin_labels)


    # This is the basin stack
    basin_stack_list = parse_list_of_list_from_string(args.basin_lists)

    # separate the data into basin csvs
    #HS.SaveHillslopeDataByBasin(this_dir, args.fname_prefix)
    #HS.SaveChannelDataByBasin(this_dir, args.fname_prefix)

    # make the directory for saving the plots
    if args.plot_directory:
        PlotDirectory = args.plot_directory
        # check if you remembered a / at the end of your path_name
        if not this_dir.endswith("/"):
            print("You forgot the '/' at the end of the plot directory, appending...")
            this_dir = this_dir+"/"

    else:
        PlotDirectory = this_dir+'hillslope_plots/'

    if not os.path.isdir(PlotDirectory):
        os.makedirs(PlotDirectory)

    # SMM: This has been tested
    if args.plot_mchi:
        print("Hello. I am going to plot chi--elevation and distance--elevation plots by basin.")
        for basin_key in these_basin_keys:
            HS.PlotChiElevationMChi(this_dir, args.fname_prefix, PlotDirectory, basin_key)
            HS.PlotLongProfileMChi(this_dir, args.fname_prefix, PlotDirectory, basin_key)

    # SMM: This has been tested
    if args.plot_segments:
        print("Hello. I am going to plot the segments in chi--elevation and distance--elevation plots by basin.")
        print("Note these are the segments linked to ridgetops and not the ones from the segementation algorithm")
        for basin_key in these_basin_keys:
            HS.PlotLongProfileSegments(this_dir, args.fname_prefix, PlotDirectory, basin_key)
            HS.PlotChiElevationSegments(this_dir, args.fname_prefix, PlotDirectory, basin_key)

    # SMM this isn't working
    if args.plot_CHT:
        for basin_key in these_basin_keys:
            HS.PlotCHTAgainstChannelData(this_dir, args.fname_prefix, PlotDirectory, basin_key)


    # SMM: This has been tested. Plots data along prfiles (i.e., L_H as a fxn of chi or flow distance)
    if args.plot_data_within_basin:
        print("Hello, I will print plots of things like L_H and other stuff as a function of position. Only the mainstem is plotted.")
        print("This plots by segment.")
        for basin_key in these_basin_keys:
            print("This basin key is: ", basin_key)
            HS.PlotHillslopeDataVsDistance(this_dir, args.fname_prefix, PlotDirectory, basin_key, args.plot_chi_by_basin, args.minimum_traces)

    # SMM Not tested since it is hard coded for mendocino
    if args.plot_mean_basin_data:
        # NOTE this requires uplift data!!
        HS.PlotHillslopeDataWithBasins(this_dir, args.fname_prefix, PlotDirectory)
        HS.PlotEStarRStarSubPlots(this_dir, args.fname_prefix, PlotDirectory, args.sc)
        #HS.PlotDataAgainstErosionRate(this_dir, args.fname_prefix, PlotDirectory)
        #HS.Make3DHillslopePlot(this_dir, args.fname_prefix, PlotDirectory)

    # SMM: again, this is hard coded for mendocino at the moment
    if args.plot_means_with_uplift:
        HS.PlotBasinDataAgainstUplift(this_dir, args.fname_prefix, PlotDirectory)

    if args.plot_hillslope_traces:
      if args.custom_plot_extent:
        HS.PlotHillslopeTraces(this_dir, args.fname_prefix, PlotDirectory, args.custom_plot_extent)
      else:
        HS.PlotHillslopeTraces(this_dir, args.fname_prefix, PlotDirectory)

    # make a plot of hillslope length distribution for specific basins (FJC)
    if args.plot_Lh_hist:
        HS.PlotHillslopeLengthDistribution(this_dir, args.fname_prefix, PlotDirectory, these_basin_keys, these_basin_labels)


    # This seems to be working but not extensively tested
    if args.determine_sc:
        HS.DetermineSc(this_dir, args.fname_prefix, PlotDirectory)

    # SMM: This has been tested
    if args.profile_plots:
        print("Hello, I will print plots of things like L_H and other stuff as a function of position. This plots things in profile.")
        if args.mainstem_only:
            print("I am only going to print data from the main stem.")
        else:
            print("All tributaries can be plotted.")
        HS.PlotChiProfileHillslopeData(this_dir, args.fname_prefix, PlotDirectory, these_basin_keys, args.plot_Ksn, args.sc,
                                       args.mainstem_only, args.minimum_traces,
                                       args.common_max_Es, args.common_max_Ksn)

    # SMM
    if args.plot_Ksn_vs_Es_Rs_by_basin:
        print("Hello, I will print plots k_sn vs E* and R*.")
        if args.mainstem_only:
            print("I am only going to print data from the main stem.")
        else:
            print("All tributaries can be plotted.")
        HS.PlotCatchmentKsnEsRs(this_dir, args.fname_prefix, PlotDirectory, these_basin_keys, args.sc, args.mainstem_only, args.minimum_traces)


    # SMM: This has been tested
    if args.plot_Es_Rs_by_basin:
        print("Let me print the E* R* plots")
        print("It makes E* R* plots for each basin and colours the data points by things like k_sn or chi.")
        for basin_key in these_basin_keys:
            HS.PlotEStarRStarWithinBasin(this_dir, args.fname_prefix, PlotDirectory, basin_key,args.minimum_traces, args.sc, args.mainstem_only, args.EsRs_colour_by)

    # SMM: Vaguely works 14-6-2018 but needs quite a bit of editing and options to be complete
    if args.plot_stacked_Es_Rs_by_chi:
        print("Let me print a stacked plot of the median E* and R* values")
        HS.PlotStackedEsRsFxnChi(this_dir, args.fname_prefix, PlotDirectory, these_basin_keys, args.sc, args.mainstem_only)

    # SMM: Working on this as of 14-6-2018
    if args.plot_clustered_Es_Rs_by_chi:
        print("Let me print a cluster plots of E* and R* values as function of chi")
        HS.PlotClusteredEsRsFxnChi(this_dir, args.fname_prefix, PlotDirectory, args.sc, args.mainstem_only, basin_stack_list)

    # SMM: Working on this as of 14-6-2018
    if args.plot_clustered_Es_Rs_plus_theoretical:
        print("Let me print a cluster plots of E* and R* values as function of chi")
        HS.PlotClusteredEsRs(this_dir, args.fname_prefix, PlotDirectory, args.sc, args.mainstem_only, basin_stack_list, args.EsRs_colour_by)

#=============================================================================
if __name__ == "__main__":
    main(sys.argv[1:])
