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

    # If there are no arguments, send to the welcome screen
    if not len(sys.argv) > 1:
        full_paramfile = print_welcome()
        sys.exit()

    # Get the arguments
    import argparse
    parser = argparse.ArgumentParser()
    # The location of the data files
    parser.add_argument("-bd", "--base_directory", type=str, help="The base directory with the m/n analysis. If this isn't defined I'll assume it's the same as the current directory.")
    parser.add_argument("-fname", "--fname_prefix", type=str, help="The prefix of your DEM WITHOUT EXTENSION!!! This must be supplied or you will get an error.")

    # What sort of analyses you want
    parser.add_argument("-PR", "--plot_rasters", type=bool, default=False, help="If this is true, I'll make raster plots of the m/n value and basin keys")
    parser.add_argument("-PC", "--plot_chi_profiles", type=bool, default=False, help="If this is true, I'll make chi-elevation plots for each basin coloured by the MLE")
    parser.add_argument("-PO", "--plot_outliers", type=bool, default=False, help="If this is true, I'll make chi-elevation plots with the outliers removed")
    parser.add_argument("-MLE", "--plot_MLE_movern", type=bool, default=False, help="If this is true, I'll make a plot of the MLE values for each m/n showing how the MLE values change as you remove the tributaries")
    parser.add_argument("-SA", "--plot_SA_data", type=bool, default=False, help="If this is true, I'll make a plot of the MLE values for each m/n showing how the MLE values change as you remove the tributaries")
    parser.add_argument("-MCMC", "--plot_MCMC", type=bool, default=False, help="If this is true, I'll make a plot of the MCMC analysis. Specify which basins you want with the -basin_keys flag.")
    parser.add_argument("-SUM", "--plot_summary", type=bool, default=False, help="If this is true, I'll make the summary CSV file and plot of the best fit m/n from each of the methods.")
    parser.add_argument("-ALL", "--all_movern_estimates", type=bool, default=False, help="If this is true, I'll make all the plots")

    # Plotting options
    parser.add_argument("-points", "--point_analysis", type=bool, default=False, help="If this is true then I'll assume that you're running the MLE analysis using the point method. Default = False")
    parser.add_argument("-show_SA_raw", "--show_SA_raw", type=bool, default=True, help="Show the raw S-A data in background of SA plot. Default = True")
    parser.add_argument("-show_SA_segments", "--show_SA_segments", type=bool, default=False, help="Show the segmented S-A data in SA plot. Default = False")
    parser.add_argument("-test_SA_regression", "--test_SA_regression", type=bool, default=False, help="If this is true I'll print the regression stats for the slope area plots.")
    parser.add_argument("-show_legend", "--show_legend", type=bool, default=False, help="If this is true, I'll display the legend for the SA plots.")

    parser.add_argument("-basin_keys", "--basin_keys",type=str,default = "", help = "This is a comma delimited string that gets the list of basins you want for the plotting. Default = no basins")

    # These control the format of your figures
    parser.add_argument("-fmt", "--FigFormat", type=str, default='png', help="Set the figure format for the plots. Default is png")
    parser.add_argument("-size", "--size_format", type=str, default='ESURF', help="Set the size format for the figure. Can be 'big' (16 inches wide), 'geomorphology' (6.25 inches wide), or 'ESURF' (4.92 inches wide) (defualt esurf).")
    parser.add_argument("-animate", "--animate", type=bool, default=True, help="If this is true I will create an animation of the chi plots. Must be used with the -PC flag set to True.")
    parser.add_argument("-keep_pngs", "--keep_pngs", type=bool, default=False, help="If this is true I will delete the png files when I animate the figures. Must be used with the -animate flag set to True.")

    args = parser.parse_args()

    if not args.fname_prefix:
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
    BasinDF = Helper.ReadBasinStatsCSV(this_dir, args.fname_prefix)
    # we need the column headers
    columns = BasinDF.columns[BasinDF.columns.str.contains('m_over_n')].tolist()
    moverns = [float(x.split("=")[-1]) for x in columns]
    start_movern = moverns[0]
    n_movern = len(moverns)
    d_movern = (moverns[-1] - moverns[0])/(n_movern-1)

    # make the plots depending on your choices
    if args.plot_rasters:
        MN.MakeRasterPlotsBasins(this_dir, args.fname_prefix, args.size_format, args.FigFormat)
        MN.MakeRasterPlotsMOverN(this_dir, args.fname_prefix, start_movern, n_movern, d_movern, size_format=args.size_format, FigFormat=args.FigFormat)
        MN.MakeRasterPlotsMOverN(this_dir, args.fname_prefix, start_movern, n_movern, d_movern, movern_method="Chi_points", size_format=args.size_format, FigFormat=args.FigFormat)
        MN.MakeRasterPlotsMOverN(this_dir, args.fname_prefix, start_movern, n_movern, d_movern, movern_method="SA", size_format=args.size_format, FigFormat=args.FigFormat)
    if args.plot_chi_profiles:
        MN.MakeChiPlotsMLE(this_dir, args.fname_prefix, basin_list=these_basin_keys, start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, size_format=args.size_format, FigFormat = args.FigFormat, animate=args.animate, keep_pngs=args.keep_pngs)
    if args.plot_outliers:
        MN.PlotProfilesRemovingOutliers(this_dir, args.fname_prefix, basin_list=these_basin_keys, start_movern=start_movern, d_movern=d_movern, n_movern=n_movern)
    if args.plot_MLE_movern:
        MN.PlotMLEWithMOverN(this_dir, args.fname_prefix,basin_list=these_basin_keys, start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, size_format=args.size_format, FigFormat = args.FigFormat)
    if args.plot_SA_data:
        SA.SAPlotDriver(this_dir, args.fname_prefix, FigFormat = args.FigFormat,size_format=args.size_format,
                        show_raw = args.show_SA_raw, show_segments = args.show_SA_segments,basin_keys = these_basin_keys)
    if args.test_SA_regression:
        #SA.TestSARegression(this_dir, args.fname_prefix)
        SA.LinearRegressionRawDataByChannel(this_dir,args.fname_prefix, basin_list=these_basin_keys)
        #SA.LinearRegressionSegmentedData(this_dir, args.fname_prefix, basin_list=these_basin_keys)
    if args.plot_MCMC:
        MN.plot_MCMC_analysis(this_dir, args.fname_prefix,basin_list=these_basin_keys, FigFormat= args.FigFormat, size_format=args.size_format)
    if args.plot_summary:
        MN.CompareMOverNEstimatesAllMethods(this_dir, args.fname_prefix, basin_list=these_basin_keys, start_movern=start_movern, d_movern=d_movern, n_movern=n_movern)
        # with SA channels
        SA_channels = True
        print ("SHOULD I SHOW THE LEGEND???? "+str(args.show_legend))
        MN.MakeMOverNSummaryPlot(this_dir, args.fname_prefix, basin_list=these_basin_keys,start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, FigFormat = args.FigFormat,size_format=args.size_format, SA_channels=SA_channels,show_legend=args.show_legend)
        # now without
        SA_channels = False
        MN.MakeMOverNSummaryPlot(this_dir, args.fname_prefix, basin_list=these_basin_keys,start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, FigFormat = args.FigFormat,size_format=args.size_format, SA_channels=SA_channels, show_legend=args.show_legend)
    if args.all_movern_estimates:
        # plot the rasters
        MN.MakeRasterPlotsBasins(this_dir, args.fname_prefix, args.size_format, args.FigFormat)
        MN.MakeRasterPlotsMOverN(this_dir, args.fname_prefix, start_movern, n_movern, d_movern, movern_method="Chi_full", size_format=args.size_format, FigFormat=args.FigFormat)
        MN.MakeRasterPlotsMOverN(this_dir, args.fname_prefix, start_movern, n_movern, d_movern, movern_method="Chi_points", size_format=args.size_format, FigFormat=args.FigFormat)
        MN.MakeRasterPlotsMOverN(this_dir, args.fname_prefix, start_movern, n_movern, d_movern, movern_method="SA", size_format=args.size_format, FigFormat=args.FigFormat)

        # make the chi plots
        MN.MakeChiPlotsMLE(this_dir, args.fname_prefix, basin_list=these_basin_keys, start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, size_format=args.size_format, FigFormat = args.FigFormat, animate=True, keep_pngs=True)

        # make the SA plots
        SA.SAPlotDriver(this_dir, args.fname_prefix, FigFormat = args.FigFormat,size_format=args.size_format,
                        show_raw = args.show_SA_raw, show_segments = True, basin_keys = these_basin_keys)
        SA.SAPlotDriver(this_dir, args.fname_prefix, FigFormat = args.FigFormat,size_format=args.size_format,
                        show_raw = args.show_SA_raw, show_segments = False, basin_keys = these_basin_keys)

        #summary plots
        MN.CompareMOverNEstimatesAllMethods(this_dir, args.fname_prefix, basin_list=these_basin_keys, start_movern=start_movern, d_movern=d_movern, n_movern=n_movern)
        # with SA channels
        SA_channels = True
        MN.MakeMOverNSummaryPlot(this_dir, args.fname_prefix, basin_list=these_basin_keys,start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, FigFormat = args.FigFormat,size_format=args.size_format, SA_channels=SA_channels, show_legend=args.show_legend)
        # now without
        SA_channels = False
        MN.MakeMOverNSummaryPlot(this_dir, args.fname_prefix, basin_list=these_basin_keys,start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, FigFormat = args.FigFormat,size_format=args.size_format, SA_channels=SA_channels, show_legend=args.show_legend)


#=============================================================================
if __name__ == "__main__":
    main(sys.argv[1:])
