#=============================================================================
# Script to plot chi analysis.
#
# Authors:
#   Simon M. Mudd
#   Fiona J. Clubb
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
import LSDMapWrappers as LSDMW
from LSDMapFigure import PlottingHelpers as phelp
import LSDPlottingTools as LSDP
from osgeo import ogr



#=============================================================================
# This is just a welcome screen that is displayed if no arguments are provided.
#=============================================================================
def print_welcome():

    print("\n\n=======================================================================")
    print("Hello! I'm going to plot chi analysis results for you.")
    print("You will need to tell me which directory to look in.")
    print("Use the -dir flag to define the working directory.")
    print("If you don't do this I will assume the data is in the same directory as this script.")
    print("For help type:")
    print("   python PlotChiAnalysis.py -h\n")
    print("=======================================================================\n\n ")

    
#=============================================================================
# Some functions for managing directories
#=============================================================================     
def MakeBasemapDirectory(this_dir):
    # check if a raster directory exists. If not then make it.
    basemap_directory = this_dir+'basemap_plots/'
    print("I am printing to a raster directory:")
    print(basemap_directory)
    if not os.path.isdir(basemap_directory):
        os.makedirs(basemap_directory)        
    
#=============================================================================
# This parses a comma separated string
#=============================================================================
def parse_list_from_string(a_string):
    """
    This just parses a comma separated string and returns an INTEGER list

    Args:
        a_string (str): The string to be parsed

    Returns:
        A list of integers

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
# This parses a dict separated string
#=============================================================================
def parse_dict_from_string(a_string):
    """
    This takes a string that is formatted to create a dict. The format is that each key/value pair is separated by a "," and each key and value are separated with a ":"

    Args:
        a_string (int): The input string

    Returns:
        A dictionary with the functions

    Author: SMM

    Date: 10/01/2018
    """
    if len(a_string) == 0:
        print("No rename dictionary found. I will return and empty dict.")
        this_rename_dict = {}
    else:
        listified_entry = [item for item in a_string.split(',')]
        this_rename_dict = {}

        # now loop through these creating a dict
        for entry in listified_entry:
            split_entry = entry.split(":")
            this_rename_dict[int(split_entry[0])]=split_entry[1]

    print("The parsed dict is: ")
    print(this_rename_dict)
    return this_rename_dict

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
# This takes the basin stack list and then gives each basin in a stack layer
# a constant value. Used for plotting.
#=============================================================================
def convert_basin_stack_to_value_dict(basin_stack_list):
    """
    This takes the basin stack list and then gives each basin in a stack layer a constant value. Used for plotting. So if there are several basin stacks each one gets a different value.

    Args:
        basin_stack_list (list of int lists): The basins that will be stacked. Each item in the list is a collection of basins that will be used in each indivdual stack plot. So, for example, if this is [ [1,2,3],[4,5]] then there will be two stacked plot, the first with basins 1,2,3 and the second with basins 4 and 5.

    Returns:
        this_value_dict (dict): A dictionary assigning a single value to each basin. Basins in the same stack will have the same value.

    Author: SMM

    Date: 11/01/2018

    """

    N_stacks = len(basin_stack_list)
    print("The number of stacks are: "+ str(N_stacks))
    if len(basin_stack_list) == 0:
        this_value_dict = {}
    else:
        this_value_dict = {}
        for idx,stack in enumerate(basin_stack_list):
            value = float(idx)/float(N_stacks)
            for item in stack:
                this_value_dict[item] = value
    return this_value_dict


#=============================================================================
# This pads an offset list so it is the same size as the basin list
#=============================================================================
def pad_offset_lists(basin_stack_list,offset_list):
    """
    This pads an offset list so it is the same size as the basin list. The offsets are the coordinate distances between the starting node of adjacent profile plots.

    Args:
        basin_stack_list (list of int lists): The basins that will be stacked. Each item in the list is a collection of basins that will be used in each indivdual stack plot. So, for example, if this is [ [1,2,3],[4,5]] then there will be two stacked plot, the first with basins 1,2,3 and the second with basins 4 and 5.
        offset_list (float list): A list of of the offset spacings for each basin stack.

    Return:
        final_offsets (float list): The locations of the offsets.

    Author: SMM

    Date: 09/01/2018
    """
    # I need to check chi the offsets
    n_basin_stacks = len(basin_stack_list)
    if len(offset_list) == 0:
        const_offset = 5
    else:
        const_offset = offset_list[-1]
    final_offsets = offset_list
    if len(offset_list) < n_basin_stacks:
        final_offsets = offset_list + [const_offset]*(n_basin_stacks - len(offset_list))
    else:
        final_offsets = offset_list

    print("Initial offsets are: ")
    print(offset_list)
    print("And const offset is: "+str(const_offset))
    print("Final offset is: ")
    print(final_offsets)

    return final_offsets



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
    parser.add_argument("-dir", "--base_directory", type=str, help="The base directory that contains your data files. If this isn't defined I'll assume it's the same as the current directory.")
    parser.add_argument("-fname", "--fname_prefix", type=str, help="The prefix of your DEM WITHOUT EXTENSION!!! This must be supplied or you will get an error (unless you're running the parallel plotting).")
    parser.add_argument("-out_fname", "--out_fname_prefix", type=str, help="The prefix of the figures WITHOUT EXTENSION!!! If not supplied the fname prefix will be used.")

    # Selecting and renaming basins
    parser.add_argument("-basin_keys", "--basin_keys",type=str,default = "", help = "This is a comma delimited string that gets the list of basins you want for the plotting. Default = no basins")
    parser.add_argument("-rename_dict", "--rename_dict",type=str,default = "", help = "This is a string that initiates a dictionary for renaming basins. The different dict entries should be comma separated, and key and value should be separated by a colon. Default = no dict")
    parser.add_argument("-basin_lists", "--basin_lists",type=str,default = "", help = "This is a string that initiates a list of a list for grouping basins. The object becomes a list of a list but the syntax is comma seperated lists, and each one is separated by a colon. Default = no dict")
    parser.add_argument("-chi_offsets", "--chi_offsets",type=str,default = "", help = "This is a string that initiates a list of chi offsets for each of the basin lists. Default = no list")
    parser.add_argument("-fd_offsets", "--flow_distance_offsets",type=str,default = "", help = "This is a string that initiates a list of flow distance offsets for each of the basin lists. Default = no list")


    # What sort of analyses you want
    parser.add_argument("-PB", "--plot_basins", type=bool, default=False, help="If this is true, I'll make a simple basin plot.")
    parser.add_argument("-PC", "--plot_chi_coord", type=bool, default=False, help="If this is true, I'll make a chi coordinate plot.")
    parser.add_argument("-all", "--all_chi_plots", type=bool, default=False, help="If this is true, I'll make all the plots including raster and chi profile plots.")
    parser.add_argument("-all_rasters", "--all_raster_plots", type=bool, default=False, help="If this is true, I'll make all the raster plots.")
    parser.add_argument("-all_stacks", "--all_stacked_plots", type=bool, default=False, help="If this is true, I'll make all the stacked plots.")

    # Some simple geographic functions that can aid in plotting regional maps. They do things like create shapefile that
    # can then be used with basemap. We don't include the basemap functions since that is not in the LSDTT toolchain (but
    # might get included later)
    parser.add_argument("-RF", "--create_raster_footprint_shapefile",type=bool, default=False, help="If true, create a shapefile from the raster. Can be used with basemap to make regional maps")
    parser.add_argument("-BM", "--create_basemap_figure",type=bool, default=False, help="If true, create a basemap file")

    # These control the format of your figures
    parser.add_argument("-fmt", "--FigFormat", type=str, default='png', help="Set the figure format for the plots. Default is png")
    parser.add_argument("-size", "--size_format", type=str, default='ESURF', help="Set the size format for the figure. Can be 'big' (16 inches wide), 'geomorphology' (6.25 inches wide), or 'ESURF' (4.92 inches wide) (defualt esurf).")
    parser.add_argument("-ar", "--figure_aspect_ratio", type=float, default=2, help="The aspect ratio of profile plots. Doesn't affect maps, whose aspect ratio is set by the size of the DEM.")
    parser.add_argument("-parallel", "--parallel", type=bool, default=False, help="If this is true I'll assume you ran the code in parallel and append all your CSVs together before plotting.")
    parser.add_argument("-dpi", "--dpi", type=int, default=250, help="The dots per inch of your figure.")
    parser.add_argument("-rotate_labels", "--rotate_labels", type=bool, default=False, help='If true I will rotate the labels of the basins on the stacked chi plots')
    parser.add_argument("-cmap", "--cmap", type=str, default="viridis", help = "The colourmap for the chi plots. Default = viridis")
    parser.add_argument("-data_fmt", "--plotting_data_format", type=str, default="log", help = "Plotting data format for chi plots. Default is log.")

    # These control the appearance of the basemap
    parser.add_argument("-bmpsm", "--basemap_parallel_spacing_multiplier", type=float, default=0.5, help="Basemap parallel spacing multiplier. Increase if parallels are too close on your basemap.")
    parser.add_argument("-bmrem", "--basemap_regional_extent_multiplier", type=float, default=4, help="Basemap regional extent multiplier. The multiple of the size of the raster to make the basemap extent")
    parser.add_argument("-bmortho", "--basemap_orthographic", type=bool, default=False, help="If this is true the basemap creates an orthographic map, that is a globe.")   
    parser.add_argument("-bmwidth", "--basemap_width_inches", type=float, default=4, help="Basemap width in inches (since matplotlib is written by yanks).")
    parser.add_argument("-bmar", "--basemap_aspect_ratio", type=float, default=1, help="Basemap aspect ratio.")    
    
    args = parser.parse_args()

    if not args.fname_prefix:
        if not args.parallel:
            print("WARNING! You haven't supplied your DEM name. Please specify this with the flag '-fname'")
            sys.exit()

    # get the base directory
    if args.base_directory:
        this_dir = args.base_directory
        # check if you remembered a / at the end of your path_name
        if not this_dir.endswith(os.sep):
            print("You forgot the separator at the end of the directory, appending...")
            this_dir = this_dir+os.sep
    else:
        this_dir = os.getcwd()

    # See if you should create a shapefile of the raster footprint
    if args.create_raster_footprint_shapefile:
        print("Let me create a shapefile of the raster footprint")

        driver_name = "ESRI shapefile"
        driver = ogr.GetDriverByName(driver_name)

        print("Driver is: ")
        print(driver)

        print("Now I'll try it from LSDPlottingTools")

        RasterFile = args.fname_prefix+".bil"
        LSDP.CreateShapefileOfRasterFootprint(this_dir, RasterFile)

    # See if you should create a basemap
    if args.create_basemap_figure:
        import LSDBasemapTools as LSDBM
        
        MakeBasemapDirectory(this_dir)
        RasterFile = args.fname_prefix+".bil"
        basemap_out_prefix = "/basemap_plots/"+out_fname_prefix
        
        # This gets the positioning
        centre_lat, centre_long, extent_lat, extent_long, xproj_extent, yproj_extent = LSDP.GetCentreAndExtentOfRaster(this_dir, RasterFile)
        
        FWI = args.basemap_width_inches
        FHI = FWI/(args.basemap_aspect_ratio)
        
        print("The basemap centrepoint is: "+str(centre_lat)+"," +str(centre_long))
        LSDBM.GenerateBasemapImageAutomated(this_dir, RasterFile, FigWidthInches = FWI, FigHeightInches = FHI, regional_extent_multiplier = args.basemap_regional_extent_multiplier, label_spacing_multiplier = args.basemap_parallel_spacing_multiplier, out_fname_prefix = basemap_out_prefix, fig_dpi = args.dpi, is_orthographic = args.basemap_orthographic)

    # See if a basin info file exists and if so get the basin list
    print("Let me check if there is a basins info csv file.")
    BasinInfoPrefix = args.fname_prefix+"_AllBasinsInfo.csv"
    BasinInfoFileName = this_dir+BasinInfoPrefix
    existing_basin_keys = []
    if os.path.isfile(BasinInfoFileName):
        print("There is a basins info csv file")
        BasinInfoDF = phelp.ReadBasinInfoCSV(this_dir, args.fname_prefix)
        existing_basin_keys = list(BasinInfoDF['basin_key'])
        existing_basin_keys = [int(x) for x in existing_basin_keys]
    else:
        print("I didn't find a basins info csv file. Check directory or filename.")

    # Parse any lists, dicts, or list of lists from the arguments
    these_basin_keys = parse_list_from_string(args.basin_keys)
    this_rename_dict = parse_dict_from_string(args.rename_dict)
    basin_stack_list = parse_list_of_list_from_string(args.basin_lists)
    chi_offset_list = parse_list_from_string(args.chi_offsets)
    fd_offset_list = parse_list_from_string(args.flow_distance_offsets)

    # If the basin keys are not supplited then assume all basins are used.
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

    # Look to see if there is a basin stack list. If there is, organise it so that we have different values in the
    # value dict
    if len(basin_stack_list) == 0:
        temp_stack = []
        temp_stack.append(these_basin_keys)
        this_value_dict = convert_basin_stack_to_value_dict(temp_stack)
    else:
        this_value_dict = convert_basin_stack_to_value_dict(basin_stack_list)

    # Now make sure all basins have a value dict value
    value_dict_single_basin = {}
    for basin in these_basin_keys:
        value_dict_single_basin[basin] = 1
        if basin not in this_value_dict:
            this_value_dict[basin] = 1

    #print("The value dict is:")
    #print(this_value_dict)

    # Now if there is a rename dict, replace the value dict values with the rename keys
    if len(this_rename_dict) != 0:
        #print("There is a rename dict. Let me adjust some values.")
        rename_value_dict = {}
        for key in this_value_dict:
            #print("Key is: "+str(key))
            if key in this_rename_dict:
                #print("I found a rename key in the value dict, changing to :"+ str(this_rename_dict[key]))
                rename_value_dict[this_rename_dict[key]] = this_value_dict[key]
            else:
                rename_value_dict[key] = this_value_dict[key]
        this_value_dict = rename_value_dict
    #print("The new value dict is: ")
    print(this_value_dict)


    # Set default offsets
    if len(chi_offset_list) == 0:
        chi_offset_list.append(5)
    if len(fd_offset_list) == 0:
        fd_offset_list.append(20000)

    #print("I am matching the offest list lengths to the number of basin stacks")
    final_chi_offsets = pad_offset_lists(basin_stack_list,chi_offset_list)
    final_fd_offsets = pad_offset_lists(basin_stack_list,fd_offset_list)



    # some formatting for the figures
    if args.FigFormat == "manuscipt_svg":
        print("You chose the manuscript svg option. This only works with the -ALL flag. For other flags it will default to simple svg")
        simple_format = "svg"
    elif args.FigFormat == "manuscript_png":
        print("You chose the manuscript png option. This only works with the -ALL flag. For other flags it will default to simple png")
        simple_format = "png"
    else:
        simple_format = args.FigFormat


    # This just plots the basins. Useful for checking on basin selection
    if args.plot_basins:
        print("I am only going to print basins.")

        # check if a raster directory exists. If not then make it.
        raster_directory = this_dir+'raster_plots/'
        if not os.path.isdir(raster_directory):
            os.makedirs(raster_directory)

        raster_out_prefix = "/raster_plots/"+args.fname_prefix
        # Now for raster plots
        # First the basins, labeled:
        LSDMW.PrintBasins_Complex(this_dir,args.fname_prefix,use_keys_not_junctions = True, show_colourbar = False,Remove_Basins = Mask_basin_keys, Rename_Basins = this_rename_dict,cmap = "jet", size_format = args.size_format,fig_format = simple_format, dpi = args.dpi, out_fname_prefix = raster_out_prefix+"_basins")

    # This plots the chi coordinate. It plots three different versions.
    # extension _CC_basins are the absins used in the chi plot
    # extension _CC_raster plots the chi raster
    # extension _CC_channels plots the channels
    if args.plot_chi_coord:
        print("I am only going to print basins.")

        # check if a raster directory exists. If not then make it.
        raster_directory = this_dir+'raster_plots/'
        if not os.path.isdir(raster_directory):
            os.makedirs(raster_directory)

        # Get the names of the relevant files
        ChannelFname = args.fname_prefix+"_chi_data_map.csv"

        raster_out_prefix = "/raster_plots/"+args.fname_prefix
        # Now for raster plots
        # First the basins, labeled:
        LSDMW.PrintBasins_Complex(this_dir,args.fname_prefix,use_keys_not_junctions = True, show_colourbar = False,Remove_Basins = Mask_basin_keys, Rename_Basins = this_rename_dict,cmap = "jet", size_format = args.size_format,fig_format = simple_format, dpi = args.dpi, out_fname_prefix = raster_out_prefix+"_CC_basins")

        # Then the chi plot for the rasters. Only call this if the masked raster exists
        masked_fname = this_dir+args.fname_prefix+"_MaskedChi.bil"
        print("\n\n\nThe filename of the chi raster is: "+masked_fname+ " I am checking if it exists.")
        import os.path as osp
        if osp.isfile(masked_fname):
            print("The chi raster exists. I'll drape the channels over the chi raster")
            LSDMW.PrintChiCoordChannelsAndBasins(this_dir,args.fname_prefix, ChannelFileName = ChannelFname, add_basin_labels = False, cmap = "cubehelix", cbar_loc = "top", size_format = args.size_format, fig_format = simple_format, dpi = args.dpi,plotting_column = "chi", colour_log = False, colorbarlabel = "$\chi$", Basin_remove_list = Mask_basin_keys, Basin_rename_dict = this_rename_dict , value_dict = this_value_dict, out_fname_prefix = raster_out_prefix+"_CC_raster", plot_chi_raster = True)
        else:
            print("The chi raster doesn't exist, I am skpping to the channel chi plots.")

        LSDMW.PrintChiCoordChannelsAndBasins(this_dir,args.fname_prefix, ChannelFileName = ChannelFname, add_basin_labels = False, cmap = "cubehelix", cbar_loc = "top", size_format = args.size_format, fig_format = simple_format, dpi = args.dpi,plotting_column = "chi", colour_log = False, colorbarlabel = "$\chi$", Basin_remove_list = Mask_basin_keys, Basin_rename_dict = this_rename_dict , value_dict = this_value_dict, out_fname_prefix = raster_out_prefix+"_CC_channels", plot_chi_raster = False)

    # This bundles a number of different analyses
    if args.all_chi_plots:
        print("You have chosen to plot all raster and stacked plots.")
        args.all_raster_plots = True
        args.all_stacked_plots = True

    # make the plots depending on your choices
    if args.all_raster_plots:
        print("I am goint to print some raster plots for you.")

        # check if a raster directory exists. If not then make it.
        raster_directory = this_dir+'raster_plots/'
        if not os.path.isdir(raster_directory):
            os.makedirs(raster_directory)

        # Get the names of the relevant files
        ChannelFname = args.fname_prefix+"_MChiSegmented.csv"

        raster_out_prefix = "/raster_plots/"+args.fname_prefix

        # Now for raster plots
        # First the basins, labeled:
        LSDMW.PrintBasins_Complex(this_dir,args.fname_prefix,use_keys_not_junctions = True, show_colourbar = False,Remove_Basins = Mask_basin_keys, Rename_Basins = this_rename_dict,cmap = "jet", size_format = args.size_format,fig_format = simple_format, dpi = args.dpi, out_fname_prefix = raster_out_prefix+"_basins")

        # Basins colour coded
        print("The value dict is: ")
        print(this_value_dict)
        LSDMW.PrintBasins_Complex(this_dir,args.fname_prefix,use_keys_not_junctions = True, show_colourbar = False,Remove_Basins = Mask_basin_keys, Rename_Basins = this_rename_dict, Value_dict = this_value_dict, cmap = "gray", size_format = args.size_format,fig_format = simple_format, dpi = args.dpi, out_fname_prefix = raster_out_prefix+"_stack_basins")

        # Now the chi steepness
        LSDMW.PrintChiChannelsAndBasins(this_dir, args.fname_prefix, ChannelFileName = ChannelFname, add_basin_labels = False, cmap = "viridis", cbar_loc = "right", size_format = args.size_format, fig_format = simple_format, dpi = args.dpi,plotting_column="m_chi",colorbarlabel = "$\mathrm{log}_{10} \; \mathrm{of} \; k_{sn}$", Basin_remove_list = Mask_basin_keys, Basin_rename_dict = this_rename_dict, value_dict = value_dict_single_basin, out_fname_prefix = raster_out_prefix+"_ksn")

        # Now plot the channels coloured by the source number
        LSDMW.PrintChiChannelsAndBasins(this_dir, args.fname_prefix, ChannelFileName = ChannelFname, add_basin_labels = False, cmap = "tab20b", cbar_loc = "None", size_format = args.size_format, fig_format = simple_format, dpi = args.dpi,plotting_column="source_key", Basin_remove_list = Mask_basin_keys, Basin_rename_dict = this_rename_dict, value_dict = this_value_dict, out_fname_prefix = raster_out_prefix+"sources", discrete_colours = True, NColours = 20, colour_log = False)

    if args.all_stacked_plots:

        # check if a chi profile directory exists. If not then make it.
        chi_profile_directory = this_dir+'chi_profile_plots/'
        if not os.path.isdir(chi_profile_directory):
            os.makedirs(chi_profile_directory)

         # Get the names of the relevant files
        ChannelFname = args.fname_prefix+"_MChiSegmented.csv"

        raster_out_prefix = "/raster_plots/"+args.fname_prefix

        print("I am going to plot some chi stacks for you.")
        cbl = "$\mathrm{log}_{10} \; \mathrm{of} \; k_{sn}$"
        i = 0
        print(basin_stack_list)
        for little_list in basin_stack_list:
            i = i+1
            this_prefix = "chi_profile_plots/Stacked_"+str(i)

            # This prints the chi profiles coloured by k_sn
            LSDMW.PrintChiStacked(this_dir, args.fname_prefix, ChannelFname, cmap = args.cmap, size_format = args.size_format, fig_format = simple_format, dpi = args.dpi,axis_data_name="chi",plot_data_name = "m_chi",colorbarlabel = cbl, plotting_data_format=args.plotting_data_format, cbar_loc = "bottom", Basin_select_list = little_list, Basin_rename_dict = this_rename_dict, out_fname_prefix = this_prefix+"_chi",X_offset = final_chi_offsets[i-1], figure_aspect_ratio = args.figure_aspect_ratio, rotate_labels = args.rotate_labels)

            # This prints channel profiles coloured by k_sn
            LSDMW.PrintChiStacked(this_dir, args.fname_prefix, ChannelFname, cmap = args.cmap, size_format = args.size_format, fig_format = simple_format, dpi = args.dpi,axis_data_name="flow_distance",plot_data_name = "m_chi", plotting_data_format =args.plotting_data_format, colorbarlabel = cbl, Basin_select_list = little_list, Basin_rename_dict = this_rename_dict, out_fname_prefix = this_prefix+"_FD", X_offset = final_fd_offsets[i-1], figure_aspect_ratio = args.figure_aspect_ratio)

            # This prints the channel profiles coloured by source number
            LSDMW.PrintChiStacked(this_dir, args.fname_prefix, ChannelFname, cmap = "tab20b", size_format = args.size_format, fig_format = simple_format, dpi = args.dpi,axis_data_name="flow_distance",plot_data_name = "source_key", plotting_data_format = 'normal', colorbarlabel = cbl, cbar_loc = "None", discrete_colours = True, NColours = 20, Basin_select_list = little_list, Basin_rename_dict = this_rename_dict, out_fname_prefix = this_prefix+"_Sources", X_offset = final_fd_offsets[i-1], figure_aspect_ratio = args.figure_aspect_ratio)


#=============================================================================
if __name__ == "__main__":
    main(sys.argv[1:])
