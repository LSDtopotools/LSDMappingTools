#=============================================================================
# Script to plot the litho information acquired from 
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
from LSDPlottingTools import LSDMap_MOverNPlotting as MN
from LSDPlottingTools import LSDMap_LithoPlotting as LP
from LSDPlottingTools import LSDMap_SAPlotting as SA
from LSDMapFigure import PlottingHelpers as Helper

#=============================================================================
# This is just a welcome screen that is displayed if no arguments are provided.
#=============================================================================
def print_welcome():

    print("\n\n=======================================================================")
    print("Hello! I'm going to plot some litho stuffs for you.")
    print("You will need to tell me which directory to look in.")
    print("Use the -wd flag to define the working directory.")
    print("If you don't do this I will assume the data is in the same directory as this script.")
    print("I also need to know some basic info like the prefix of your base dem and the name of your lithologic DEM")
    print("For help type:")
    print("   python PlotLithoInfo.py -h\n")
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
    parser.add_argument("-fname", "--fname_prefix", type=str, help="The prefix of your DEM WITHOUT EXTENSION!!! This must be supplied or you will get an error.")

    # Basin selection stuffs
    parser.add_argument("-basin_keys", "--basin_keys",type=str,default = "", help = "This is a comma delimited string that gets the list of basins you want for the plotting. Default = no basins")

    # What sort of analyses you want
    parser.add_argument("-c", "--check", type=bool, default=True, help="Turn to false if you really know what you are doing, like specify all the non-automatic file names for example")
    parser.add_argument("-lk", "--lithokey_file", type=str, default="", help="This is in case you wanna manually specify the lithokey file")
    parser.add_argument("-leg", "--legend", type=str, default="", help="Turn to True to print a separate legend file with the geology lithokey equivalent")
    parser.add_argument("-LM", "--LithoMap", type=str, default="", help="Turn to True to print a Lithologic map with the basins keys")
    parser.add_argument("-BLM", "--Basic_LithoMap", type=str, default="", help="Turn to True to print a basic Lithologic map won the top of a hillshade")
    parser.add_argument("-mn", "--movern", type=bool, default = False, help="Turn to True to turn to True and provide the movern plots related to lithology")
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
    


    # some formatting for the figures
    if args.FigFormat == "manuscipt_svg":
        print("You chose the manuscript svg option. This only works with the -ALL flag. For other flags it will default to simple svg")
        simple_format = "svg"
    elif args.FigFormat == "manuscript_png":
        print("You chose the manuscript png option. This only works with the -ALL flag. For other flags it will default to simple png")
        simple_format = "png"
    else:
        simple_format = args.FigFormat

    if args.check:
        print("I am now checking your files from the rasterization.")
        dict_file = LP.litho_pre_check(this_dir,args.lithokey_file, fname = args.fname_prefix)

    if args.LithoMap:
        LP.MakeRasterLithoBasinMap(this_dir, args.fname_prefix, args.fname_prefix+"_LITHRAST", dict_file["lithodict"], size_format='ESURF', FigFormat=args.FigFormat, basins = True)
    if args.Basic_LithoMap:
        LP.MakeRasterLithoBasinMap(this_dir, args.fname_prefix, args.fname_prefix+"_LITHRAST", dict_file["lithodict"], size_format='ESURF', FigFormat=args.FigFormat, basins = False)
    
    if args.legend:
        print ("ongoing work on the legend")

    if args.movern:
        # get the range of moverns, needed for plotting
        BasinDF = Helper.ReadBasinStatsCSV(this_dir, args.fname_prefix)
        # we need the column headers
        columns = BasinDF.columns[BasinDF.columns.str.contains('m_over_n')].tolist()
        moverns = [float(x.split("=")[-1]) for x in columns]
        start_movern = moverns[0]
        n_movern = len(moverns)
        d_movern = (moverns[-1] - moverns[0])/(n_movern-1)
        LP.movern_two_litho(args.fname_prefix,this_dir, litho = [9,10], basin_list=these_basin_keys,start_movern=start_movern, d_movern=d_movern, n_movern=n_movern)
        LP.MakeRasterLithoBasinMap(this_dir, args.fname_prefix, args.fname_prefix+"_LITHRAST", dict_file["lithodict"], size_format='ESURF', FigFormat=args.FigFormat)
        #quit()
        #MN.MakeRasterPlotsMOverN(this_dir, args.fname_prefix, start_movern, n_movern, d_movern, size_format=args.size_format, FigFormat=simple_format,lith = True)
        #MN.MakeRasterPlotsMOverN(this_dir, args.fname_prefix, start_movern, n_movern, d_movern, movern_method="Chi_points", size_format=args.size_format, FigFormat=simple_format,lith = True)
        #MN.MakeRasterPlotsMOverN(this_dir, args.fname_prefix, start_movern, n_movern, d_movern, movern_method="SA", size_format=args.size_format, FigFormat=simple_format,lith = True)
        LP.MakeChiPlotsByLith(this_dir, args.fname_prefix, basin_list=these_basin_keys,plot_colorbar = False, start_movern=start_movern, d_movern=d_movern, n_movern=n_movern, size_format=args.size_format, FigFormat=simple_format, animate=args.animate, keep_pngs=True)


    

#=============================================================================
if __name__ == "__main__":
    main(sys.argv[1:])
