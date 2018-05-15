#=============================================================================
# Script to plot the knickpoint data produced with the LSDTT knickpoint picking algorithm 
#
# Authors:
#     Boris Gailleton, Fiona J. Clubb
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
from LSDPlottingTools import LSDMap_BasicMaps as BP
from LSDMapFigure import PlottingHelpers as Helper
from LSDPlottingTools import LSDMap_KnickpointPlotting as KP
from LSDPlottingTools import LSDMap_ChiPlotting as CP
from LSDPlottingTools import LSDMap_LithoPlotting as LP


#=============================================================================
# This is just a welcome screen that is displayed if no arguments are provided.
#=============================================================================
def print_welcome():

    print("\n\n=======================================================================")
    print("Hello! I'm going to plot some knickpoints stuffs for you.")
    print("You will need to tell me which directory to look in.")
    print("Use the -wd flag to define the working directory.")
    print("If you don't do this I will assume the data is in the same directory as this script.")
    print("I also need to know the common prefix of all your files generated whith LSDTopoTool")
    print("For help type:")
    print("   python PlotknickpointAnalysis.py -h\n")
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
    parser.add_argument("-dir", "--base_directory", type=str, help="The base directory with the m/n analysis. If this isn't defined I'll assume it's the same as the current directory.")
    parser.add_argument("-fname", "--fname_prefix", type=str, help="The prefix of your DEM WITHOUT EXTENSION!!! This must be supplied or you will get an error.")
    # Basin and source selection
    # Basins selection
    parser.add_argument("-basin_keys", "--basin_keys",type=str,default = "", help = "This is a comma delimited string that gets the list of basins you want for the plotting. Default = all basins")
    # Sources selection
    parser.add_argument("-source_keys", "--source_keys",type=str,default = "", help = "This is a comma delimited string that gets the list of sources you want for the plotting. Default = all sources")
    parser.add_argument("-min_sl", "--min_source_length", type=float , default = 0, help = "This is a minimum length for the river to plot, if you want to only plot the river profile of the main rivers for example. Default = 0 (no restrictions)")
    parser.add_argument("-main_stem", "--isolate_main_stem", type=bool , default =False, help = "set to True to only keep the main stem of each basins")

    # These control the format of your figures
    parser.add_argument("-fmt", "--FigFormat", type=str, default='png', help="Set the figure format for the plots. Default is png")
    parser.add_argument("-size", "--size_format", type=str, default='ESURF', help="Set the size format for the figure. Can be 'big' (16 inches wide), 'geomorphology' (6.25 inches wide), or 'ESURF' (4.92 inches wide) (defualt esurf).")

    # ALL
    parser.add_argument("-all", "--AllAnalysis", type=bool, default = False, help="Turn on to have fun")
    parser.add_argument("-allD", "--AllAnalysisDebug", type=bool, default = False, help="Turn on to have even more fun")


    # Mchi_related
    parser.add_argument("-mcstd", "--mchi_map_std", type=bool, default = False, help="Turn to True to plot a standart M_chi map on an HS. Small reminder, Mchi = Ksn if calculated with A0 = 1.")
    parser.add_argument("-mcbk", "--mchi_map_black", type=bool, default = False, help="Turn to True to plot a standart M_chi map on Black background. Small reminder, Mchi = Ksn if calculated with A0 = 1.")
    parser.add_argument("-minmc", "--min_mchi_map", type=int, default = 0, help="mininum value for the scale of your m_chi maps, default 0")
    parser.add_argument("-maxmc", "--max_mchi_map", type=int, default = 0, help="maximum value for the scale of your m_chi maps, default auto")
    parser.add_argument("-almc", "--alpha_mchi", type=float, default = 0, help="alpha_value for the background mchi raster")


    #knickpint related

    parser.add_argument("-ksnPs", "--ksn_per_source", type=bool, default = False, help="Print one figure per source key selected, with ksn -> f(chi & flow_distance) in the folder .../river_plots/. it displays the ksn out of Mudd et al., 2014 method, and the TVD one out of the *insert algorithm name*")
    parser.add_argument("-rivplot", "--river_profile", type=bool, default = False, help="Print one figure per source key selected, with elevation -> f(chi & flow_distance) in the folder .../river_plots/. it displays river profiles in a chi and distance spaces")
    parser.add_argument("-basplot", "--basin_plot", type=bool, default = False, help="Print one figure per basins key selected, with elevation -> f(chi & flow_distance) in the folder .../river_plots/. it displays river profiles in a chi and distance spaces")
    parser.add_argument("-rasplot", "--raster_plots", type = bool, default = False, help="Print raster plots with knickpoints on top of ksn in the folder .../raster_plots/")
    parser.add_argument("-rasplot_ld", "--raster_plots_large_dataset", type = bool, default = False, help="Print raster plots with knickpoints on top of ksn in the folder .../raster_plots/")
    parser.add_argument("-statplot", "--statistical_plots", type = bool, default = False, help="Print a bunch of statistics about the knickpoints in the folder .../raster_plots/")
    parser.add_argument("-stradivarius", "--multi_violin_plots", type = bool, default = False, help="Print a bunch of statistical distribution against elevation, chi, ... for the selected knickpoints")

    # Others
    parser.add_argument("-nbh", "--n_bin_hist", type = int, default = 0, help = "Customization of the number of bin you want for the general histogram. Default is an automatic in-built selection from numpy")
    parser.add_argument("-cov", "--cut_off_val", type = str, default = "0,0,0,0", help = "Cutoff value for the knickpoint magnitude (the drop/increase of ksn). Default is 0 (no cut)")
    parser.add_argument("-kal", "--kalib", type = bool, default = False, help = "Don't use that.")
    parser.add_argument("-segelev", "--print_segmented_elevation", type = bool, default = False, help = "This print the segmented elevation on the top of the river profiles, in transparent black. Useful to check segment boundaries and adjust target_nodes parameter. Default False.")
    parser.add_argument("-extent_rast_cmap", "--manual_extent_colormap_knickpoint_raster", type = str, default = "", help = "This print the segmented elevation on the top of the river profiles, in transparent black. Useful to check segment boundaries and adjust target_nodes parameter. Default False.")
    parser.add_argument("-size_kp_map", "--size_kp_map", type = bool, default = True, help = "This print the segmented elevation on the top of the river profiles, in transparent black. Useful to check segment boundaries and adjust target_nodes parameter. Default False.")
    parser.add_argument("-max_hist", "--maximum_extent_for_histogram", type = int, default = 0, help = "This print the segmented elevation on the top of the river profiles, in transparent black. Useful to check segment boundaries and adjust target_nodes parameter. Default False.")
    parser.add_argument("-lith_rast","--lithologic_raster", type = bool, default = False, help = "switch on if you have a _LITHRAST raster, it will plot a hillshade colored by lithologic unit")
    parser.add_argument("-save","--save_output", type = bool, default = False, help = "switch on if you have the willingness to save your selected knickpoints in a new csv file named prefix_output.csv")
    parser.add_argument("-coeff_size_kp", "--kp_coeff_size", type = float, default = 1, help = "qualitative size of knickpoints on river profile. Default 10, increase or decrease to adapt the size of the triangles")
    parser.add_argument("-fixed_size_kp_min_max_river", "--fixed_size_kp_min_max_river", type = str, default = "", help = "qualitative size of knickpoints on river profile. Default 50, increase or decrease to adapt the size of the triangles")
    parser.add_argument("-normelev_rel", "--normalise_elevation_to_outlet_relative", type = bool, default = False, help = "Normalize the elevation for each watershed, considering the outlet points as 0, but without changing the relief.")
    parser.add_argument("-normelev_abs", "--normalise_elevation_to_outlet_absolute", type = bool, default = False, help = "Normalize the elevation for each watershed, considering the outlet points as 0 and the maximum to 1.")
    parser.add_argument("-unicol_kp", "--unicolor_for_knickpoint_map", type = str, default = None, help = "use this flag to force a uni color for knickpoints on map. It has to be a matplotlib color code (simplest way is to use html colorcode like #AABB47)")
    parser.add_argument("-GS", "--general_size",type=str,default = "", help = "This is a comma delimited string used to scale your knickpoints relatively: min_value,min_size,max_value,max_size. see documentation for explanations")

    args = parser.parse_args()

    # Processing the basin/source keys selection
    print("I am reading the basin/source key selection and checking your parameters...")
    
    if len(args.basin_keys) == 0:
        print("No basins found, I will plot all of them")
        these_basin_keys = []
    else:
        these_basin_keys = [int(item) for item in args.basin_keys.split(',')]
        print("The basins I will plot are:")
        print(these_basin_keys)

    if len(args.source_keys) == 0:
        print("No sources found, I will plot all of them")
        these_source_keys = []
    else:
        these_source_keys = [int(item) for item in args.source_keys.split(',')]
        print("The sources preselected are:")
        print(these_source_keys)

    if len(args.fixed_size_kp_min_max_river) != 2:
        min_max_kp_river = []
    else:
        min_max_kp_river = [float(item) for item in args.min_max_kp_river.split(',')]
        print("You are recasting the size for all your knickpoints triangles plotted on the river_plots:")
        print(min_max_kp_river)

    if len(args.manual_extent_colormap_knickpoint_raster) > 0:
        manual_cmap_extent_raster_plot = [float(item) for item in args.manual_extent_colormap_knickpoint_raster.split(',')]
        print("You choose a manual colorbar for plotting:")
        print(manual_cmap_extent_raster_plot)
    else:
        manual_cmap_extent_raster_plot = []
        
    if(args.maximum_extent_for_histogram>0):
        ext_dseg_hist = [0,args.maximum_extent_for_histogram]
        ext_dksn_hist = [-args.maximum_extent_for_histogram,args.maximum_extent_for_histogram]
    else:
        ext_dseg_hist = []
        ext_dksn_hist = []

    if len(args.general_size) == 0:
        print("Sizing the knickpoints automatically")
        general_size_kp = []
    else:
        general_size_kp = [float(item) for item in args.general_size.split(',')]
        print(general_size_kp)


    print("Getting your cut off values...")
    if(args.cut_off_val != "auto"):
        try:
            covfefe = [float(item) for item in args.cut_off_val.replace(" ", "").split(',')]
            print("ok.")
            covfefe_t = [-covfefe[0],covfefe[1],-10000,covfefe[2]]
            covfefe = covfefe_t
        except ValueError:
            print("Something went wrong - I am defaulting the values")
            covfefe = [0,0,-10000,0]
        print("cut off values:")
        print(covfefe)
    else:
        covfefe = "auto"
        print("I will choose automatically the cut-off values, based on the quartiles")

    # Processing the size choice
    try:
        size = [int(item) for item in args.size_format.split(',')]
    except ValueError:
        size = args.size_format

    if not args.fname_prefix:
        print("WARNING! You haven't supplied your DEM name. Please specify this with the flag '-fname'")
        sys.exit()
    print("Done")

    # Normalisation???
    if(args.normalise_elevation_to_outlet_absolute and args.normalise_elevation_to_outlet_relative):
        print("Erm, you cannot normaliSe your elevation in relatively absolute or absolutely relative way! you have to choose.")
        quit()
    elif(args.normalise_elevation_to_outlet_absolute):

        normalisation = "absolute"
    elif(args.normalise_elevation_to_outlet_relative):
        normalisation = "relative"
    else:
        normalisation = None



    print("Loading the dataset:")

    KI = KP.KP_plotting(args.base_directory,args.fname_prefix, basin_key = these_basin_keys, source_key = these_source_keys, min_length = args.min_source_length, 
        cut_off_val = covfefe, main_stem = args.isolate_main_stem, normalisation = normalisation, size_kp = general_size_kp)

    
    if(args.AllAnalysisDebug):
        args.AllAnalysis = True
        args.ksn_per_source = True


    if(args.AllAnalysis):
        args.statistical_plots = True
        args.river_profile = True
        args.raster_plots = True
        # args.raster_plots_large_dataset = True
        args.basin_plot = True


    # Plotting hte knickpoints
    if(args.statistical_plots):
        
        if(int(args.n_bin_hist) == 0):
            n_b = "auto"
        else:
            n_b = int(args.n_bin_hist)

        KI.print_histogram(size = size, format = args.FigFormat, n_bin = n_b,x_extents = ext_dksn_hist)
        KI.print_histogram(size = size, format = args.FigFormat, n_bin = n_b, data = "delta_segelev",x_extents = ext_dseg_hist)

        KI.print_box_and_whisker(size = size, format = args.FigFormat, label_size = 8, binning = 'source_key', facecolor = "white", grid = True)
        KI.print_box_and_whisker(size = size, format = args.FigFormat, label_size = 8, binning = 'basin_key', facecolor = "white", grid = True)

    if(args.ksn_per_source):
        print("Printing a set of ksn values with the knickpoints and their magnitude in a Chi distance")
        KI.print_ksn_profile(size = size, format = args.FigFormat, x_axis = "chi", knickpoint = True, title = "auto", label_size = 8, facecolor = 'white', legend = True)
        # KI.print_ksn_profile(size = size, format = args.FigFormat, x_axis = "chi",y_axis = "segmented_elevation", knickpoint = True, title = "auto", label_size = 8, facecolor = 'white', legend = True)


        # print("Printing a set of ksn values with the knickpoints and their magnitude in a Flow distance")
        # KI.print_ksn_profile(size = size, format = args.FigFormat, x_axis = "flow_distance", knickpoint = True, title = "auto", label_size = 8, facecolor = 'white', legend = True)

    if(args.river_profile):
        print("Printing river profiles in chi spaces")
        KI.print_river_profile(size = size, format = args.FigFormat, x_axis = "chi", knickpoint = True, title = "auto", label_size = 8, facecolor = 'white', kalib = args.kalib,
         print_seg_elev = args.print_segmented_elevation, coeff_size = args.kp_coeff_size, size_recasting = min_max_kp_river)
        print("Printing river profiles in flow distance")
        KI.print_river_profile(size = size, format = args.FigFormat, x_axis = "flow_distance", knickpoint = True, title = "auto", label_size = 8, facecolor = 'white', kalib = args.kalib,
         print_seg_elev = args.print_segmented_elevation, coeff_size = args.kp_coeff_size, size_recasting = min_max_kp_river)
        print("Printing river profiles for the entire basins")

    if (args.basin_plot):
        KI.print_river_profile(size = size, format = args.FigFormat, x_axis = "flow_distance", knickpoint = True, title = "auto", label_size = 8, facecolor = 'white', binning = "basin_key", 
            kalib = args.kalib, print_seg_elev = args.print_segmented_elevation)
        KI.print_river_profile(size = size, format = args.FigFormat, x_axis = "chi", knickpoint = True, title = "auto", label_size = 8, facecolor = 'white', binning = "basin_key", 
            kalib = args.kalib, print_seg_elev = args.print_segmented_elevation)



    if(args.raster_plots):
        KI.print_map_topo(size = size, format = args.FigFormat,label_size = 8, return_fig = False, extent_cmap = [], kalib = False)
        KI.print_map_of_kp(unicolor_kp = args.unicolor_for_knickpoint_map,size = size, format = args.FigFormat, black_bg = False, scale_points = True, label_size = 6,
            size_kp = args.size_kp_map, extent_cmap = manual_cmap_extent_raster_plot, kalib = args.kalib)
        KI.print_map_of_kp(unicolor_kp = args.unicolor_for_knickpoint_map,size = size, format = args.FigFormat, black_bg = True, scale_points = True, label_size = 6,
            size_kp = args.size_kp_map, extent_cmap = manual_cmap_extent_raster_plot, kalib = args.kalib)

    if(args.raster_plots_large_dataset):
        KI.print_map_of_kp(unicolor_kp = args.unicolor_for_knickpoint_map,size = size, format = args.FigFormat, black_bg = False, scale_points = False, label_size = 6, 
            size_kp = args.size_kp_map, extent_cmap = manual_cmap_extent_raster_plot, kalib = args.kalib)
        KI.print_map_of_kp(unicolor_kp = args.unicolor_for_knickpoint_map,size = size, format = args.FigFormat, black_bg = True, scale_points = False, label_size = 6, 
            size_kp = args.size_kp_map, extent_cmap = manual_cmap_extent_raster_plot, kalib = args.kalib)

    if(args.multi_violin_plots):
        KI.stradivarius_analysis(size = size, format = args.FigFormat)


    # Preparing the min_max color for mchi maps
    if(args.max_mchi_map <= args.min_mchi_map):
        colo = []
    else:
        colo = [args.min_mchi_map,args.max_mchi_map]
    
    if(args.lithologic_raster):
        dict_file = LP.litho_pre_check(args.base_directory,"", fname = args.fname_prefix)
        LP.MakeRasterLithoBasinMap(args.base_directory, args.fname_prefix, args.fname_prefix+"_LITHRAST", dict_file["lithodict"], size_format= size, FigFormat=args.FigFormat, 
            m_chi = True, mancol = colo, log_scale_river = True, minmax_m_chi = [0.01,1])
        cml = LP.getLithoColorMap(args.fname_prefix, args.base_directory)
        KI.print_map_of_kp(unicolor_kp = args.unicolor_for_knickpoint_map,size = size, format = args.FigFormat, black_bg = False, scale_points = False, label_size = 6,
            size_kp = args.size_kp_map, extent_cmap = manual_cmap_extent_raster_plot, kalib = args.kalib, cml = cml, lith_raster = True)


    if(args.save_output):
        KI.save_output_csv()

    # Preparing the min_max color for mchi maps
    if(args.max_mchi_map <= args.min_mchi_map):
        colo = []
    else:
        colo = [args.min_mchi_map,args.max_mchi_map]

    if args.mchi_map_std:
        
        CP.map_Mchi_standard(args.base_directory, args.fname_prefix, size_format=args.size_format, FigFormat=args.FigFormat, basin_list = these_basin_keys, log = False, 
            colmanscal = colo, knickpoint = True, alpha_background = args.alpha_mchi)

    if args.mchi_map_black:
        
        CP.map_Mchi_standard(args.base_directory, args.fname_prefix, size_format=args.size_format, FigFormat=args.FigFormat, basin_list = these_basin_keys, log = False, 
            colmanscal = colo, bkbg = True, knickpoint = True, alpha_background = args.alpha_mchi)


#=============================================================================
if __name__ == "__main__":
    main(sys.argv[1:])

