#=============================================================================
# Script to plot the m/n analysis using the MLE collinearity
#
# Authors:
#     Fiona J. Clubb
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
from decimal import Decimal
from LSDPlottingTools import LSDMap_MOverNPlotting as MN
from LSDPlottingTools import LSDMap_SAPlotting as SA
from LSDMapFigure import PlottingHelpers as Helper

#=============================================================================
# This is just a welcome screen that is displayed if no arguments are provided.
#=============================================================================
def print_welcome():

    print("\n\n=======================================================================")
    print("Hello! I'm going to plot the m/n analysis results for you.")
    print("You will need to tell me which directory to look in.")
    print("Use the -dir flag to define the working directory.")
    print("If you don't do this I will assume the data is in the same directory as this script.")
    print("You also need the -fname flag which will give the prefix of the raster files.")
    print("See our documentation for computing the data needed for these visualisation scripts:")
    print("https://lsdtopotools.github.io/LSDTT_documentation/LSDTT_chi_analysis.html#_calculating_concavity")
    print("For help type:")
    print("   python PlotMOverNAnalysis.py -h\n")
    print("=======================================================================\n\n ")


#=============================================================================
# This parses a comma separated string
#=============================================================================
def parse_list_from_string(a_string):
    """
    This just parses a comma separated string and returns either a float or an int. 
    Will return a float if there is a decimal in any of the data entries

    Args:
        a_string (str): The string to be parsed

    Returns:
        A list of integers or floats

    Author: SMM

    Date: 10/01/2018
    """
    print("Hey pardner, I'm a gonna parse a string into a list. Yeehaw.")
    if len(a_string) == 0:
        print("No items found, I am returning and empty list.")
        return_list = []
    elif "." in a_string:
        print("I found a decimal in your string. I am going to assume these are floats")
        return_list = [float(item) for item in a_string.split(',')]
        print("The parsed string is:")
        print(return_list)
    else:
        return_list = [int(item) for item in a_string.split(',')]
        print("The parsed string is:")
        print(return_list)

    return return_list    
 
    
#=============================================================================
# This parses a comma separated string into strings
#=============================================================================
def parse_string_list_from_string(a_string):
    """
    This just parses a comma separated string and returns either a float or an int. 
    Will return a float if there is a decimal in any of the data entries

    Args:
        a_string (str): The string to be parsed

    Returns:
        A list of integers

    Author: SMM

    Date: 18/11/2019
    """
    print("Hey pardner, I'm a gonna parse a string into a list. Yeehaw.")
    if len(a_string) == 0:
        print("No items found, I am returning and empty list.")
        return_list = []
    else:
        return_list = [item for item in a_string.split(',')]
        print("The parsed string is:")
        print(return_list)

    return return_list       
    
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
# Returns true if the data only uses the disorder metric
#============================================================================= 
def check_if_disorder_metric_only(this_dir, this_fname_prefix):
    Chi_disorder_only = False
    from pathlib import Path
    config = Path(this_dir+this_fname_prefix+'_movernstats_0.5_fullstats.csv')
    if config.is_file():
        print("The bootstrap files exist")
    else:
        print("I didn't find the bootstrap files.")
        Chi_disorder_only = True
    return Chi_disorder_only
    
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

    # What sort of analyses you want
    parser.add_argument("-PR", "--plot_rasters", type=bool, default=False, help="If this is true, I'll make raster plots of the m/n value and basin keys")
    parser.add_argument("-chi", "--plot_basic_chi", type=bool, default=False, help="If this is true I'll make basin chi plots for each basin coloured by elevation.")
    parser.add_argument("-PC", "--plot_chi_profiles", type=bool, default=False, help="If this is true, I'll make chi-elevation plots for each basin coloured by the MLE")
    parser.add_argument("-K", "--plot_chi_by_K", type=bool, default=False, help="If this is true, I'll make chi-elevation plots for each basin coloured by K. NOTE - you MUST have a column in your chi csv with the K value or this will break!")
    parser.add_argument("-pcbl", "--plot_chi_by_lith", type=bool, default=False, help="If this is true, I'll make chi-elevation plots for each basin coloured by litho. NOTE - you MUST have a column in your chi csv with the K value or this will break!")
    parser.add_argument("-PO", "--plot_outliers", type=bool, default=False, help="If this is true, I'll make chi-elevation plots with the outliers removed")
    parser.add_argument("-MLE", "--plot_MLE_movern", type=bool, default=False, help="If this is true, I'll make a plot of the MLE values for each m/n showing how the MLE values change as you remove the tributaries")
    parser.add_argument("-SA", "--plot_SA_data", type=bool, default=False, help="If this is true, I'll make a plot of the MLE values for each m/n showing how the MLE values change as you remove the tributaries")
    parser.add_argument("-MCMC", "--plot_MCMC", type=bool, default=False, help="If this is true, I'll make a plot of the MCMC analysis. Specify which basins you want with the -basin_keys flag.")
    parser.add_argument("-pts", "--point_uncertainty", type=bool, default=False, help="If this is true, I'll make a plot of the range in m/n from the MC points analysis")
    parser.add_argument("-hist", "--plot_histogram", type=bool, default=False, help="If this is true, I'll make plots of the pdfs of m/n values for each method.")
    parser.add_argument("-disorder", "--plot_disorder", type=bool, default=False, help="If this is true, I'll make plots of the chi disorder analysis.")
    parser.add_argument("-SUM", "--plot_summary", type=bool, default=False, help="If this is true, I'll make the summary CSV file and plot of the best fit concavity from each of the methods.")
    parser.add_argument("-ALL", "--all_movern_estimates", type=bool, default=False, help="If this is true, I'll make all the plots")
    parser.add_argument("-DisFxnDist", "--disorder_function_of_distance", type=bool, default=False, help="If this is true, I'll make a plot of the disorder metric as a function of position")


    # Plotting options
    parser.add_argument("-points", "--point_analysis", type=bool, default=False, help="If this is true then I'll assume that you're running the MLE analysis using the point method. Default = False")
    parser.add_argument("-show_SA_raw", "--show_SA_raw", type=bool, default=True, help="Show the raw S-A data in background of SA plot. Default = True")
    parser.add_argument("-show_SA_segments", "--show_SA_segments", type=bool, default=False, help="Show the segmented S-A data in SA plot. Default = False")
    parser.add_argument("-test_SA_regression", "--test_SA_regression", type=bool, default=False, help="If this is true I'll print the regression stats for the slope area plots.")
    parser.add_argument("-show_legend", "--show_legend", type=bool, default=True, help="If this is true, I'll display the legend for plots.")
    parser.add_argument("-no_legend", "--no_legend", dest="show_legend", action="store_false", help="Flag to not display legends, I'll not display the legend for plots, default is for legend to be displayed. Note taht setting show_legend False does not achieve this due to bool issues with python parsing")

    # Options about basin selection
    parser.add_argument("-basin_keys", "--basin_keys",type=str,default = "", help = "This is a comma delimited string that gets the list of basins you want for the plotting. Default = no basins")
    parser.add_argument("-basin_lists", "--basin_lists",type=str,default = "", help = "This is a string that initiates a list of a list for grouping basins. The object becomes a list of a list but the syntax is comma seperated lists, and each one is separated by a colon. Default = no dict")
    parser.add_argument("-group_names", "--group_names",type=str,default = "", help = "Names of the groups provided by basin_lists. Used in legends")
    
    
    # These control the format of your figures
    parser.add_argument("-fmt", "--FigFormat", type=str, default='png', help="Set the figure format for the plots. Default is png")
    parser.add_argument("-size", "--size_format", type=str, default='ESURF', help="Set the size format for the figure. Can be 'big' (16 inches wide), 'geomorphology' (6.25 inches wide), or 'ESURF' (4.92 inches wide) (defualt esurf).")
    parser.add_argument("-animate", "--animate", type=bool, default=True, help="If this is true I will create an animation of the chi plots. Must be used with the -PC flag set to True.")
    parser.add_argument("-keep_pngs", "--keep_pngs", type=bool, default=False, help="If this is true I will delete the png files when I animate the figures. Must be used with the -animate flag set to True.")
    parser.add_argument("-parallel", "--parallel", type=bool, default=False, help="If this is true I'll assume you ran the code in parallel and append all your CSVs together before plotting.")

    args = parser.parse_args()

    print(argv)
    print(args)

    if not args.fname_prefix:
        if not args.parallel:
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

 
    # See if a basin info file exists and if so get the basin list
    print("Let me check if there is a basins info csv file.")
    BasinInfoPrefix = args.fname_prefix+"_AllBasinsInfo.csv"
    BasinInfoFileName = this_dir+BasinInfoPrefix
    existing_basin_keys = []
    if os.path.isfile(BasinInfoFileName):
        print("There is a basins info csv file")
        BasinInfoDF = Helper.ReadBasinInfoCSV(this_dir, args.fname_prefix)
        existing_basin_keys = list(BasinInfoDF['basin_key'])
        existing_basin_keys = [int(x) for x in existing_basin_keys]
    else:
        print("I didn't find a basins info csv file. Check directory or filename.")

    # Parse any lists, dicts, or list of lists from the arguments
    these_basin_keys = parse_list_from_string(args.basin_keys)
    basin_stack_list = parse_list_of_list_from_string(args.basin_lists)
    basin_stack_names = parse_string_list_from_string(args.group_names)


    # If the basin keys are not supplied then assume all basins are used.
    if these_basin_keys == []:
        these_basin_keys = existing_basin_keys

    # Python is so amazing. Look at the line below.
    Mask_basin_keys = [i for i in existing_basin_keys if i not in these_basin_keys]
    print("All basins are: ")
    print(existing_basin_keys)
    print("The basins to keep are:")
    print(these_basin_keys)
    print("The basins to mask are:")
    print(Mask_basin_keys)

    # This is an old version. It passes empty strings to the plotting functions. 
    #if len(args.basin_keys) == 0:
    #    print("No basins found, I will plot all of them")
    #    # Note that if you pass an empty list to the plotting functions, they will plot all the basins
    #    these_basin_keys = []
    #else:
    #    these_basin_keys = [int(item) for item in args.basin_keys.split(',')]
    #    print("The basins I will plot are:")
    #    print(these_basin_keys)     
        
        

        
    # This checks to see if chi points method is being used. 
    # If not, assumes only the disorder metric has been calculated
    Using_disorder_metric_only = check_if_disorder_metric_only(this_dir, args.fname_prefix)

    
    if not args.parallel:
        BasinDF = Helper.ReadBasinStatsCSV(this_dir, args.fname_prefix)
    else:
        BasinDF = Helper.AppendBasinCSVs(this_dir, args.fname_prefix)

        # if parallel, get the fname from the data directory. This assumes that your directory is called
        # something sensible that relates to the DEM name.
        split_fname = this_dir.split("/")
        split_fname = split_fname[len(split_fname)-2]


    # get the range of moverns, needed for plotting
    # we need the column headers
    columns = BasinDF.columns[BasinDF.columns.str.contains('m_over_n')].tolist()
    moverns = [float(x.split("=")[-1]) for x in columns]
    start_movern = moverns[0]
    n_movern = len(moverns)
    x = Decimal((moverns[-1] - moverns[0])/(n_movern-1))
    d_movern = round(x,2)
    print ('Start movern, n_movern, d_movern: ')
    print (start_movern, n_movern, d_movern)

    # some formatting for the figures
    if args.FigFormat == "manuscipt_svg":
        print("You chose the manuscript svg option. This only works with the -ALL flag. For other flags it will default to simple svg")
        simple_format = "svg"
    elif args.FigFormat == "manuscript_png":
        print("You chose the manuscript png option. This only works with the -ALL flag. For other flags it will default to simple png")
        simple_format = "png"
    else:
        simple_format = args.FigFormat


    # make the plots depending on your choices
    if args.plot_rasters:
        MN.MakeRasterPlotsBasins(this_dir, args.fname_prefix, args.size_format, simple_format, parallel=args.parallel)
    if args.plot_basic_chi:
        MN.MakePlotsWithMLEStats(this_dir, args.fname_prefix, basin_list=these_basin_keys, start_movern=start_movern, d_movern=d_movern, n_movern=n_movern,parallel=args.parallel)
 
    if args.plot_chi_profiles:    
        if Using_disorder_metric_only:
            MN.MakeChiPlotsChi(this_dir, args.fname_prefix, basin_list=these_basin_keys, 
                           start_movern=start_movern, d_movern=d_movern, n_movern=n_movern,
                           size_format=args.size_format, FigFormat = args.FigFormat, animate=True, keep_pngs=True,parallel=args.parallel)
        else:
            MN.MakeChiPlotsMLE(this_dir, args.fname_prefix, basin_list=these_basin_keys, 
                           start_movern=start_movern, d_movern=d_movern, n_movern=n_movern,
                           size_format=args.size_format, FigFormat = args.FigFormat, animate=True, keep_pngs=True,parallel=args.parallel)  
    
    if args.plot_chi_by_K:
        MN.MakeChiPlotsColouredByK(this_dir, args.fname_prefix, basin_list=these_basin_keys, start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, size_format=args.size_format, FigFormat=simple_format, animate=args.animate, keep_pngs=args.keep_pngs, parallel=args.parallel)
    if args.plot_chi_by_lith:
        MN.MakeChiPlotsColouredByLith(this_dir, args.fname_prefix, basin_list=these_basin_keys, start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, size_format=args.size_format, FigFormat=simple_format, animate=args.animate, keep_pngs=args.keep_pngs,parallel=args.parallel)
    if args.plot_outliers:
        MN.PlotProfilesRemovingOutliers(this_dir, args.fname_prefix, basin_list=these_basin_keys, start_movern=start_movern, d_movern=d_movern, n_movern=n_movern,parallel=args.parallel)
    if args.plot_MLE_movern:
        MN.PlotMLEWithMOverN(this_dir, args.fname_prefix,basin_list=these_basin_keys, start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, size_format=args.size_format, FigFormat =simple_format,parallel=args.parallel)
    if args.plot_SA_data:
        SA.SAPlotDriver(this_dir, args.fname_prefix, FigFormat = simple_format,size_format=args.size_format,
                        show_raw = args.show_SA_raw, show_segments = args.show_SA_segments,basin_keys = these_basin_keys, parallel=args.parallel)
    if args.test_SA_regression:
        #SA.TestSARegression(this_dir, args.fname_prefix)
        SA.LinearRegressionRawDataByChannel(this_dir,args.fname_prefix, basin_list=these_basin_keys, parallel=args.parallel)
        #SA.LinearRegressionSegmentedData(this_dir, args.fname_prefix, basin_list=these_basin_keys)
    if args.plot_MCMC:
        MN.plot_MCMC_analysis(this_dir, args.fname_prefix,basin_list=these_basin_keys, FigFormat= simple_format, size_format=args.size_format,parallel=args.parallel)
    if args.point_uncertainty:
        MN.PlotMCPointsUncertainty(this_dir, args.fname_prefix,basin_list=these_basin_keys, FigFormat=simple_format, size_format=args.size_format,start_movern=start_movern, d_movern=d_movern, n_movern=n_movern,parallel=args.parallel)
    if args.plot_histogram:
        MN.MakeMOverNSummaryHistogram(this_dir, args.fname_prefix,basin_list=these_basin_keys,start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, FigFormat=simple_format, size_format=args.size_format, show_legend=args.show_legend,Chi_disorder=True)

    if args.plot_summary:
        # This function creates a csv that has the concavity statistics in it
        MN.CompareMOverNEstimatesAllMethods(this_dir, args.fname_prefix, basin_list=these_basin_keys,
                                            start_movern=start_movern, d_movern=d_movern,
                                            n_movern=n_movern, parallel=args.parallel, Chi_disorder=True)
        
        MN.MakeMOverNSummaryPlot(this_dir, args.fname_prefix, basin_list=these_basin_keys,
                                 start_movern=start_movern, d_movern=d_movern,
                                 n_movern=n_movern, FigFormat = simple_format,size_format=args.size_format, show_legend=args.show_legend,parallel=args.parallel, Chi_disorder=True)

        # This only prints the summary plots for bootstrap and disorder metrics
        MN.MakeMOverNSummaryPlot(this_dir, args.fname_prefix, basin_list=these_basin_keys,
                                 start_movern=start_movern, d_movern=d_movern,
                                 n_movern=n_movern, FigFormat = simple_format,size_format=args.size_format,
                                 show_legend=args.show_legend,parallel=args.parallel,
                                 Chi_all = False, SA_raw = False, SA_segmented = False,
                                 SA_channels = False, Chi_bootstrap = True, Chi_disorder=True)

        MN.MakeMOverNSummaryHistogram(this_dir, args.fname_prefix,basin_list=these_basin_keys,
                                      start_movern=start_movern, d_movern=d_movern,
                                      n_movern=n_movern, FigFormat=args.FigFormat, size_format=args.size_format, show_legend=args.show_legend, Chi_disorder=True)
    if args.plot_disorder:
        MN.MakeRasterPlotsMOverN(this_dir, args.fname_prefix, start_movern, n_movern, d_movern, movern_method="Chi_disorder", size_format=args.size_format, FigFormat=args.FigFormat,parallel=args.parallel)
        
        # This function creates a csv that has the concavity statistics in it
        MN.CompareMOverNEstimatesAllMethods(this_dir, args.fname_prefix, basin_list=these_basin_keys, start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, parallel=args.parallel, Chi_disorder=True)
        
        MN.MakeMOverNSummaryPlot(this_dir, args.fname_prefix, basin_list=these_basin_keys,start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, FigFormat = simple_format,size_format=args.size_format, show_legend=args.show_legend,parallel=args.parallel, Chi_disorder=True)
        
        MN.MakeMOverNSummaryHistogram(this_dir, args.fname_prefix,basin_list=these_basin_keys,start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, FigFormat=args.FigFormat, size_format=args.size_format, show_legend=args.show_legend, Chi_disorder=True)

    if args.all_movern_estimates:
        print("I am going to print out loads and loads of figures for you.")
        # plot the rasters
        MN.MakeRasterPlotsBasins(this_dir, args.fname_prefix, args.size_format, args.FigFormat,parallel=args.parallel)
        
        if not Using_disorder_metric_only:
            MN.MakeRasterPlotsMOverN(this_dir, args.fname_prefix, start_movern, n_movern, d_movern, 
                                 movern_method="Chi_full", size_format=args.size_format,
                                 FigFormat=args.FigFormat,parallel=args.parallel)
            MN.MakeRasterPlotsMOverN(this_dir, args.fname_prefix, start_movern, n_movern, d_movern,
                                 movern_method="Chi_points", size_format=args.size_format,
                                 FigFormat=args.FigFormat,parallel=args.parallel)
            
        MN.MakeRasterPlotsMOverN(this_dir, args.fname_prefix, start_movern, n_movern, d_movern, 
                                 movern_method="SA", size_format=args.size_format,
                                 FigFormat=args.FigFormat,parallel=args.parallel)
        MN.MakeRasterPlotsMOverN(this_dir, args.fname_prefix, start_movern, n_movern, 
                                 d_movern, movern_method="Chi_disorder",
                                 size_format=args.size_format, FigFormat=args.FigFormat,parallel=args.parallel)

        # make the chi plots
        if Using_disorder_metric_only:
            MN.MakeChiPlotsChi(this_dir, args.fname_prefix, basin_list=these_basin_keys, 
                           start_movern=start_movern, d_movern=d_movern, n_movern=n_movern,
                           size_format=args.size_format, FigFormat = args.FigFormat, animate=True, keep_pngs=True,parallel=args.parallel)
        else:
            MN.MakeChiPlotsMLE(this_dir, args.fname_prefix, basin_list=these_basin_keys, 
                           start_movern=start_movern, d_movern=d_movern, n_movern=n_movern,
                           size_format=args.size_format, FigFormat = args.FigFormat, animate=True, keep_pngs=True,parallel=args.parallel)

        # make the SA plots
        SA.SAPlotDriver(this_dir, args.fname_prefix, FigFormat = args.FigFormat,size_format=args.size_format,
                        show_raw = args.show_SA_raw, show_segments = True, basin_keys = these_basin_keys, parallel=args.parallel)
        SA.SAPlotDriver(this_dir, args.fname_prefix, FigFormat = args.FigFormat,size_format=args.size_format,
                        show_raw = args.show_SA_raw, show_segments = False, basin_keys = these_basin_keys, parallel=args.parallel)

        #summary plots
        # This function creates a csv that has the concavity statistics in it
        MN.CompareMOverNEstimatesAllMethods(this_dir, args.fname_prefix, basin_list=these_basin_keys,
                                            start_movern=start_movern, d_movern=d_movern,
                                            n_movern=n_movern, parallel=args.parallel, Chi_disorder=True)
        
        MN.MakeMOverNSummaryPlot(this_dir, args.fname_prefix, basin_list=these_basin_keys,start_movern=start_movern, 
                                 d_movern=d_movern,
                                 n_movern=n_movern, FigFormat = simple_format,size_format=args.size_format, show_legend=args.show_legend,parallel=args.parallel, Chi_disorder=True)

        # This only prints the summary plots for bootstrap and disorder metrics
        MN.MakeMOverNSummaryPlot(this_dir, args.fname_prefix, basin_list=these_basin_keys,
                                 start_movern=start_movern, d_movern=d_movern,
                                 n_movern=n_movern, FigFormat = simple_format,size_format=args.size_format,
                                 show_legend=args.show_legend,parallel=args.parallel,
                                 Chi_all = False, SA_raw = False, SA_segmented = False,
                                 SA_channels = False, Chi_bootstrap = True, Chi_disorder=True)

        MN.MakeMOverNSummaryHistogram(this_dir, args.fname_prefix,basin_list=these_basin_keys,
                                      start_movern=start_movern, d_movern=d_movern,
                                      n_movern=n_movern, FigFormat=args.FigFormat, size_format=args.size_format, show_legend=args.show_legend, Chi_disorder=True)


        
    if args.disorder_function_of_distance:
        # This function creates a csv that has the concavity statistics in it
        print("=====================================================")
        print("=====================================================")
        print("\n\n\n\nI am going to get the summary information.")
        
        
        # See if the summary already exists
        
        
        print("Let me check if there is a concavity summary csv file.")
        SummaryPrefix = args.fname_prefix+"_movern_summary.csv"
        SummaryFileName = this_dir+"summary_plots/"+SummaryPrefix
        print("The summary filename is: "+SummaryFileName)
        if os.path.isfile(SummaryFileName):
            print("There is already a summary file")
        else:
            print("No summray csv found. I will calculate a new one.")
            MN.CompareMOverNEstimatesAllMethods(this_dir, args.fname_prefix, basin_list=these_basin_keys,
                                            start_movern=start_movern, d_movern=d_movern,
                                            n_movern=n_movern, parallel=args.parallel, Chi_disorder=True)
        
        # Okay, now we plot the metrics as a function of distance
        print("I am going to print the following lists of basins: ")
        print(basin_stack_list)
        
        MN.MakeMOverNDisorderDistancePlot(this_dir, args.fname_prefix, basin_list_list=basin_stack_list,
                                 start_movern=start_movern, d_movern=d_movern,
                                 n_movern=n_movern, FigFormat = simple_format,size_format=args.size_format,
                                 show_legend=args.show_legend,parallel=args.parallel,group_names=basin_stack_names)
        
        
        
        
#=============================================================================
if __name__ == "__main__":
    main(sys.argv[1:])
