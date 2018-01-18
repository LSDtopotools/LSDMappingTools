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
# This parses a comma separated string
#=============================================================================    
def parse_list_from_string(a_string):
    if len(a_string) == 0:
        print("No items found, I am returning and empty list.")
        return_list = []
    else:
        return_list = [int(item) for item in a_string.split(',')]
        print("The parsed string is:")
        print(return_list)
        
    return return_list

#=============================================================================
# This parses a dict separated string
#=============================================================================    
def parse_dict_from_string(a_string):
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
# This pads an offset list so it is the same size as the basin list
#=============================================================================     
def pad_offset_lists(basin_stack_list,offset_list):
    
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
    parser.add_argument("-dir", "--base_directory", type=str, help="The base directory with the m/n analysis. If this isn't defined I'll assume it's the same as the current directory.")
    parser.add_argument("-fname", "--fname_prefix", type=str, help="The prefix of your DEM WITHOUT EXTENSION!!! This must be supplied or you will get an error (unless you're running the parallel plotting).")
    parser.add_argument("-out_fname", "--out_fname_prefix", type=str, help="The prefix of the figures WITHOUT EXTENSION!!! If not supplied the fname prefix will be used.")

    # Selecting and renaming basins
    parser.add_argument("-basin_keys", "--basin_keys",type=str,default = "", help = "This is a comma delimited string that gets the list of basins you want for the plotting. Default = no basins")  
    parser.add_argument("-rename_dict", "--rename_dict",type=str,default = "", help = "This is a string that initiates a dictionary for renaming basins. The different dict entries should be comma separated, and key and value should be separated by a colon. Default = no dict")   
    parser.add_argument("-basin_lists", "--basin_lists",type=str,default = "", help = "This is a string that initiates a list of a list for grouping basins. The object becomes a list of a list but the syntax is comma seperated lists, and each one is separated by a colon. Default = no dict")
    parser.add_argument("-chi_offsets", "--chi_offsets",type=str,default = "", help = "This is a string that initiates a list of chi offsets for each of the basin lists. Default = no list")
    parser.add_argument("-fd_offsets", "--flow_distance_offsets",type=str,default = "", help = "This is a string that initiates a list of flow distance offsets for each of the basin lists. Default = no list")
    
    
    # What sort of analyses you want
    parser.add_argument("-PR", "--plot_rasters", type=bool, default=False, help="If this is true, I'll make raster plots of the m/n value and basin keys")
    parser.add_argument("-PCS", "--plot_chi_stack", type=bool, default=False, help="If this is true, I'll make a stack of chi plots.")

    parser.add_argument("-all", "--all_chi_plots", type=bool, default=False, help="If this is true, I'll make all the plots including raster and chi profile plots.")
    parser.add_argument("-all_rasters", "--all_raster_plots", type=bool, default=False, help="If this is true, I'll make all the raster plots.")
    parser.add_argument("-all_stacks", "--all_stacked_plots", type=bool, default=False, help="If this is true, I'll make all the stacked plots.")

    # These control the format of your figures
    parser.add_argument("-fmt", "--FigFormat", type=str, default='png', help="Set the figure format for the plots. Default is png")
    parser.add_argument("-size", "--size_format", type=str, default='ESURF', help="Set the size format for the figure. Can be 'big' (16 inches wide), 'geomorphology' (6.25 inches wide), or 'ESURF' (4.92 inches wide) (defualt esurf).")
    parser.add_argument("-parallel", "--parallel", type=bool, default=False, help="If this is true I'll assume you ran the code in parallel and append all your CSVs together before plotting.")
    parser.add_argument("-dpi", "--dpi", type=int, default=250, help="The dots per inch of your figure.")
    
    args = parser.parse_args()

    if not args.fname_prefix:
        if not args.parallel:
            print("WARNING! You haven't supplied your DEM name. Please specify this with the flag '-fname'")
            sys.exit()

            
    
    
    # Parse any lists, dicts, or list of lists from the arguments   
    these_basin_keys = parse_list_from_string(args.basin_keys)
    this_rename_dict = parse_dict_from_string(args.rename_dict)
    basin_stack_list = parse_list_of_list_from_string(args.basin_lists)
    chi_offset_list = parse_list_from_string(args.chi_offsets)
    fd_offset_list = parse_list_from_string(args.flow_distance_offsets)

    # Set default offsets
    if len(chi_offset_list) == 0:
        chi_offset_list.append(5)
    if len(fd_offset_list) == 0:
        fd_offset_list.append(20000)
    
    print("I am matching the offest list lengths to the number of basin stacks")    
    final_chi_offsets = pad_offset_lists(basin_stack_list,chi_offset_list)
    final_fd_offsets = pad_offset_lists(basin_stack_list,fd_offset_list)

    # get the base directory
    if args.base_directory:
        this_dir = args.base_directory
        # check if you remembered a / at the end of your path_name
        if not this_dir.endswith("/"):
            print("You forgot the '/' at the end of the directory, appending...")
            this_dir = this_dir+"/"
    else:
        this_dir = os.getcwd()

       
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
    
    print("Let me get a list of the masked basins")
       
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
    
    # get a value dict for a single colour basin
    value_dict_single_basin = {}
    for basin in these_basin_keys:
        value_dict_single_basin[basin] = 1

    # some formatting for the figures
    if args.FigFormat == "manuscipt_svg":
        print("You chose the manuscript svg option. This only works with the -ALL flag. For other flags it will default to simple svg")
        simple_format = "svg"
    elif args.FigFormat == "manuscript_png":
        print("You chose the manuscript png option. This only works with the -ALL flag. For other flags it will default to simple png")
        simple_format = "png"
    else:
        simple_format = args.FigFormat

    ChannelFname = args.fname_prefix+"_MChiSegmented.csv"
 


        
    # This bundles a number of different analyses    
    if args.all_stacked_plots:
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
        LSDMW.PrintBasins_Complex(this_dir,args.fname_prefix,use_keys_not_junctions = True, show_colourbar = False,Remove_Basins = Mask_basin_keys, Rename_Basins = this_rename_dict,cmap = "jet", cbar_loc = "right", size_format = args.size_format,fig_format = simple_format, dpi = args.dpi, out_fname_prefix = raster_out_prefix+"_basins")
        
        # Now the chi steepness
        LSDMW.PrintChiChannelsAndBasins(this_dir, args.fname_prefix, ChannelFileName = ChannelFname, add_basin_labels = False, cmap = "viridis", cbar_loc = "right", size_format = args.size_format, fig_format = simple_format, dpi = args.dpi,plotting_column="m_chi",colorbarlabel = "$\mathrm{log}_{10} \; \mathrm{of} \; k_{sn}$", Basin_remove_list = Mask_basin_keys, Basin_rename_dict = this_rename_dict, value_dict = value_dict_single_basin, out_fname_prefix = raster_out_prefix+"ksn")
        
        # Now plot the channels coloured by the source number
        LSDMW.PrintChiChannelsAndBasins(this_dir, args.fname_prefix, ChannelFileName = ChannelFname, add_basin_labels = False, cmap = "tab20b", cbar_loc = "None", size_format = args.size_format, fig_format = simple_format, dpi = args.dpi,plotting_column="source_key", Basin_remove_list = Mask_basin_keys, Basin_rename_dict = this_rename_dict, value_dict = value_dict_single_basin, out_fname_prefix = raster_out_prefix+"sources", discrete_colours = True, NColours = 20, colour_log = False)
        
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
        this_rename_dict = {}
        for little_list in basin_stack_list:
            i = i+1
            this_prefix = "chi_profile_plots/Stacked_"+str(i)+"_" 
            
            LSDMW.PrintChiStacked(this_dir, args.fname_prefix, ChannelFname, cmap = "viridis", size_format = "ESURF", fig_format = "png", dpi = 250,axis_data_name="chi",plot_data_name = "m_chi",colorbarlabel = cbl, cbar_loc = "bottom", Basin_select_list = little_list, Basin_rename_dict = this_rename_dict, out_fname_prefix = this_prefix+"_chi")
        
            LSDMW.PrintChiStacked(this_dir, args.fname_prefix, ChannelFname, cmap = "viridis", size_format = "ESURF", fig_format = "png", dpi = 250,axis_data_name="flow_distance",plot_data_name = "m_chi", plotting_data_format = 'log', colorbarlabel = cbl, Basin_select_list = little_list, Basin_rename_dict = this_rename_dict, out_fname_prefix = this_prefix+"_FD")    

            LSDMW.PrintChiStacked(this_dir, args.fname_prefix, ChannelFname, cmap = "viridis", size_format = "ESURF", fig_format = "png", dpi = 250,axis_data_name="flow_distance",plot_data_name = "source_key", plotting_data_format = 'normal', colorbarlabel = cbl, cbar_loc = "None", discrete_colours = True, NColours = 15, Basin_select_list = little_list, Basin_rename_dict = this_rename_dict, out_fname_prefix = this_prefix+"_Sources")        


#=============================================================================
if __name__ == "__main__":
    main(sys.argv[1:])
