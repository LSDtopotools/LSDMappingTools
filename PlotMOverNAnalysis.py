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
    print("Use the -wd flag to define the working directory.")
    print("If you don't do this I will assume the data is in the same directory as this script.")
    print("For help type:")
    print("   python PlotMOverNAnalysis.py -h\n")
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
    parser.add_argument("-dir", "--base_directory", type=str, help="The base directory with the m/n analysis. If this isn't defined I'll assume it's the same as the current directory.")
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
    # parser.add_argument("-basin_joyplot", "--basin_joyplot", type=bool, default=False, help="If this is true, I'll make a joyplot showing m/n for each basin from the chi points")
    parser.add_argument("-SUM", "--plot_summary", type=bool, default=False, help="If this is true, I'll make the summary CSV file and plot of the best fit m/n from each of the methods.")
    parser.add_argument("-ALL", "--all_movern_estimates", type=bool, default=False, help="If this is true, I'll make all the plots")

    # Plotting options
    parser.add_argument("-points", "--point_analysis", type=bool, default=False, help="If this is true then I'll assume that you're running the MLE analysis using the point method. Default = False")
    parser.add_argument("-show_SA_raw", "--show_SA_raw", type=bool, default=True, help="Show the raw S-A data in background of SA plot. Default = True")
    parser.add_argument("-show_SA_segments", "--show_SA_segments", type=bool, default=False, help="Show the segmented S-A data in SA plot. Default = False")
    parser.add_argument("-test_SA_regression", "--test_SA_regression", type=bool, default=False, help="If this is true I'll print the regression stats for the slope area plots.")
    parser.add_argument("-show_legend", "--show_legend", type=bool, default=True, help="If this is true, I'll display the legend for the SA plots.")

    parser.add_argument("-basin_keys", "--basin_keys",type=str,default = "", help = "This is a comma delimited string that gets the list of basins you want for the plotting. Default = no basins")

    # These control the format of your figures
    parser.add_argument("-fmt", "--FigFormat", type=str, default='png', help="Set the figure format for the plots. Default is png")
    parser.add_argument("-size", "--size_format", type=str, default='ESURF', help="Set the size format for the figure. Can be 'big' (16 inches wide), 'geomorphology' (6.25 inches wide), or 'ESURF' (4.92 inches wide) (defualt esurf).")
    parser.add_argument("-animate", "--animate", type=bool, default=True, help="If this is true I will create an animation of the chi plots. Must be used with the -PC flag set to True.")
    parser.add_argument("-keep_pngs", "--keep_pngs", type=bool, default=False, help="If this is true I will delete the png files when I animate the figures. Must be used with the -animate flag set to True.")
    parser.add_argument("-parallel", "--parallel", type=bool, default=False, help="If this is true I'll assume you ran the code in parallel and append all your CSVs together before plotting.")

    args = parser.parse_args()

    if not args.fname_prefix:
        if not args.parallel:
            print("WARNING! You haven't supplied your DEM name. Please specify this with the flag '-fname'")
            sys.exit()

    # check the basins
    print("You told me that the basin keys are: ")
    print(args.basin_keys)

    if len(args.basin_keys) == 0:
        print("No basins found, I will plot all of them")
        these_basin_keys = []
    else:
        these_basin_keys = [int(item) for item in args.basin_keys.split(',')]
        print("The basins I will plot are:")
        print(these_basin_keys)

    # get the base directory
    if args.base_directory:
        this_dir = args.base_directory
    else:
        this_dir = os.getcwd()

    # get the range of moverns, needed for plotting
    if not args.parallel:
        BasinDF = Helper.ReadBasinStatsCSV(this_dir, args.fname_prefix)
    else:
        BasinDF = Helper.AppendBasinCSVs(this_dir)
        
        # if parallel, get the fname from the data directory. This assumes that your directory is called
        # something sensible that relates to the DEM name.
        split_fname = this_dir.split("/")
        split_fname = split_fname[len(split_fname)-2]
        #args.fname_prefix = split_fname # commented out for now since base fname given, basins will always have basinX fname_prefix
        

    # we need the column headers
    columns = BasinDF.columns[BasinDF.columns.str.contains('m_over_n')].tolist()
    moverns = [float(x.split("=")[-1]) for x in columns]
    start_movern = moverns[0]
    n_movern = len(moverns)
    d_movern = (moverns[-1] - moverns[0])/(n_movern-1)

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
        MN.MakeRasterPlotsMOverN(this_dir, args.fname_prefix, start_movern, n_movern, d_movern, size_format=args.size_format, FigFormat=simple_format, parallel=args.parallel)
        MN.MakeRasterPlotsMOverN(this_dir, args.fname_prefix, start_movern, n_movern, d_movern, movern_method="Chi_points", size_format=args.size_format, FigFormat=simple_format,parallel=args.parallel)
        MN.MakeRasterPlotsMOverN(this_dir, args.fname_prefix, start_movern, n_movern, d_movern, movern_method="SA", size_format=args.size_format, FigFormat=simple_format,parallel=args.parallel)
    if args.plot_basic_chi:
        MN.MakePlotsWithMLEStats(this_dir, args.fname_prefix, basin_list=these_basin_keys, start_movern=start_movern, d_movern=d_movern, n_movern=n_movern,parallel=args.parallel)
    if args.plot_chi_profiles:
        MN.MakeChiPlotsMLE(this_dir, args.fname_prefix, basin_list=these_basin_keys, start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, size_format=args.size_format, FigFormat = simple_format, animate=args.animate, keep_pngs=args.keep_pngs, parallel=args.parallel)
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
                        show_raw = args.show_SA_raw, show_segments = args.show_SA_segments,basin_keys = these_basin_keys)
    if args.test_SA_regression:
        #SA.TestSARegression(this_dir, args.fname_prefix)
        SA.LinearRegressionRawDataByChannel(this_dir,args.fname_prefix, basin_list=these_basin_keys)
        #SA.LinearRegressionSegmentedData(this_dir, args.fname_prefix, basin_list=these_basin_keys)
    if args.plot_MCMC:
        MN.plot_MCMC_analysis(this_dir, args.fname_prefix,basin_list=these_basin_keys, FigFormat= simple_format, size_format=args.size_format,parallel=args.parallel)
    if args.point_uncertainty:
        MN.PlotMCPointsUncertainty(this_dir, args.fname_prefix,basin_list=these_basin_keys, FigFormat=simple_format, size_format=args.size_format,start_movern=start_movern, d_movern=d_movern, n_movern=n_movern,parallel=args.parallel)
    if args.plot_histogram:
        MN.MakeMOverNSummaryHistogram(this_dir, args.fname_prefix,basin_list=these_basin_keys,start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, FigFormat=simple_format, size_format=args.size_format, show_legend=args.show_legend, parallel=args.parallel)
    if args.plot_summary:
        MN.CompareMOverNEstimatesAllMethods(this_dir, args.fname_prefix, basin_list=these_basin_keys, start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, parallel=args.parallel)
        MN.MakeMOverNSummaryPlot(this_dir, args.fname_prefix, basin_list=these_basin_keys,start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, FigFormat = simple_format,size_format=args.size_format, show_legend=args.show_legend,parallel=args.parallel)
        MN.MakeMOverNPlotOneMethod(this_dir,args.fname_prefix,basin_list=these_basin_keys,start_movern=start_movern,d_movern=d_movern,n_movern=n_movern,FigFormat=args.FigFormat,size_format=args.size_format)
    # if args.basin_joyplot:
    #     MN.CompareMOverNEstimatesAllMethods(this_dir, args.fname_prefix, basin_list=these_basin_keys, start_movern=start_movern, d_movern=d_movern, n_movern=n_movern)
    #     MN.MakeBasinJoyplot(this_dir, args.fname_prefix, basin_list=these_basin_keys, FigFormat=simple_format, size_format=args.size_format)
    if args.all_movern_estimates:
        # plot the rasters
        MN.MakeRasterPlotsBasins(this_dir, args.fname_prefix, args.size_format, args.FigFormat,parallel=args.parallel)
        MN.MakeRasterPlotsMOverN(this_dir, args.fname_prefix, start_movern, n_movern, d_movern, movern_method="Chi_full", size_format=args.size_format, FigFormat=args.FigFormat,parallel=args.parallel)
        MN.MakeRasterPlotsMOverN(this_dir, args.fname_prefix, start_movern, n_movern, d_movern, movern_method="Chi_points", size_format=args.size_format, FigFormat=args.FigFormat,parallel=args.parallel)
        MN.MakeRasterPlotsMOverN(this_dir, args.fname_prefix, start_movern, n_movern, d_movern, movern_method="SA", size_format=args.size_format, FigFormat=args.FigFormat,parallel=args.parallel)

        # make the chi plots
        MN.MakeChiPlotsMLE(this_dir, args.fname_prefix, basin_list=these_basin_keys, start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, size_format=args.size_format, FigFormat = args.FigFormat, animate=True, keep_pngs=True,parallel=args.parallel)

        # make the SA plots
        SA.SAPlotDriver(this_dir, args.fname_prefix, FigFormat = args.FigFormat,size_format=args.size_format,
                        show_raw = args.show_SA_raw, show_segments = True, basin_keys = these_basin_keys)
        SA.SAPlotDriver(this_dir, args.fname_prefix, FigFormat = args.FigFormat,size_format=args.size_format,
                        show_raw = args.show_SA_raw, show_segments = False, basin_keys = these_basin_keys)

        #summary plots
        MN.CompareMOverNEstimatesAllMethods(this_dir, args.fname_prefix, basin_list=these_basin_keys, start_movern=start_movern, d_movern=d_movern, n_movern=n_movern,parallel=args.parallel)
        MN.MakeMOverNSummaryPlot(this_dir, args.fname_prefix, basin_list=these_basin_keys,start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, FigFormat = args.FigFormat,size_format=args.size_format, show_legend=args.show_legend)
        #joyplot
        MN.MakeMOverNPlotOneMethod(this_dir,args.fname_prefix,basin_list=these_basin_keys,start_movern=start_movern,d_movern=d_movern,n_movern=n_movern,FigFormat=args.FigFormat,size_format=args.size_format)
        MN.MakeMOverNSummaryHistogram(this_dir, args.fname_prefix,basin_list=these_basin_keys,start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, FigFormat=args.FigFormat, size_format=args.size_format, show_legend=args.show_legend)


#=============================================================================
if __name__ == "__main__":
    main(sys.argv[1:])
