# Driver to make the terrace profile plots
# FJC 01/11/17

# import modules
import sys
import os

from LSDPlottingTools import LSDMap_TerracePlotting as TerracePlotter

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
    parser.add_argument("-dir", "--base_directory", type=str, help="The base directory with the terrace analysis. If this isn't defined I'll assume it's the same as the current directory.")
    parser.add_argument("-fname", "--fname_prefix", type=str, help="The prefix of your DEM WITHOUT EXTENSION!!! This must be supplied or you will get an error.")
    parser.add_argument("-basin_jns", "--junction_list", type=str, help="If you want, you can pass in a list of junction numbers which I will use to plot the terraces. Useful if you ran the terrace code in parallel. Needs to be a comma separated string!")

    # What sort of analyses you want to do
    parser.add_argument("-LP", "--long_profiler", type=bool, default=False, help="If this is true, I'll make plots of the terrace long profiles (Default = true)")
    parser.add_argument("-PR", "--plot_rasters", type=bool, default=False, help="If this is true, I'll make raster plots of the terrace locations (Default=false)")
    parser.add_argument("-HM", "--heat_map", type=bool, default=False, help="If this is true, I'll make a heat map of the terrace pixel locations")
    parser.add_argument("-dips", "--dips", type=bool,default=False, help="If this is true, I'll calculate the dip and dip direction of each terrace.")
    parser.add_argument("-DT", "--digitised_terraces", type=bool,default=False, help="If this is true I'll filter the terrace points using a shapefile of digitised terraces.")
    parser.add_argument("-shp", "--shapefile_name", type=str, default=None, help="The shapefile of digitised terraces. Must be supplied if you want to filter terraces by shapefile, obvz.")

    parser.add_argument("-ksn", "--colour_by_ksn", type=bool, default=False, help="If this is true I'll colour the main stem channel by ksn. Can only be used if you have the csv file ending in '_MChiSegmented.csv'")

    # These control the format of your figures
    parser.add_argument("-fmt", "--FigFormat", type=str, default='png', help="Set the figure format for the plots. Default is png")
    parser.add_argument("-size", "--size_format", type=str, default='ESURF', help="Set the size format for the figure. Can be 'big' (16 inches wide), 'geomorphology' (6.25 inches wide), or 'ESURF' (4.92 inches wide) (defualt esurf).")

    args = parser.parse_args()

    # get the base directory
    if args.base_directory:
        this_dir = args.base_directory
        # check if you remembered a / at the end of your path_name
        if not this_dir.endswith("/"):
            print("You forgot the '/' at the end of the directory, appending...")
            this_dir = this_dir+"/"
    else:
        this_dir = os.getcwd()

    # check if you supplied the DEM prefix
    if not args.fname_prefix:
        print("WARNING! You haven't supplied your DEM name. Please specify this with the flag '-fname'")
        sys.exit()

    if args.long_profiler:
        if not args.digitised_terraces:
            TerracePlotter.long_profiler_dist(this_dir, args.fname_prefix)
        else:
            #TerracePlotter.long_profiler_dist(this_dir, args.fname_prefix, digitised_terraces=True, shapefile_name = args.shapefile_name)
            TerracePlotter.long_profiler_centrelines(this_dir,args.fname_prefix,args.shapefile_name, args.colour_by_ksn)
            TerracePlotter.MakeTerracePlotChiSpace(this_dir, args.fname_prefix,args.shapefile_name)

    if args.plot_rasters:
        TerracePlotter.MakeRasterPlotTerraceIDs(this_dir, args.fname_prefix, args.FigFormat, args.size_format)
        TerracePlotter.MakeRasterPlotTerraceElev(this_dir, args.fname_prefix, args.FigFormat, args.size_format)
    if args.heat_map:
        TerracePlotter.MakeTerraceHeatMap(this_dir,args.fname_prefix, args.fname_prefix, prec=150, FigFormat=args.FigFormat)
        TerracePlotter.MakeTerraceHeatMapNormalised(this_dir,args.fname_prefix, args.fname_prefix, prec=150, FigFormat=args.FigFormat)
    if args.dips:
        TerracePlotter.write_dip_and_dipdir_to_csv(this_dir,args.fname_prefix, args.digitised_terraces, args.shapefile_name)
        # TerracePlotter.MakeRasterPlotTerraceDips(this_dir,args.fname_prefix,FigFormat=args.FigFormat,size_format=args.size_format)
    if len(args.junction_list) > 0:
        print("You passed in a list of basin junctions. I'll append these to your filename, and make heat plots for each one.")
        jn_list = args.junction_list.split(",")
        for jn in jn_list:
            this_fname = args.fname_prefix+"_"+jn
            print ("This fname is: ", this_fname)
            TerracePlotter.MakeTerraceHeatMap(this_dir,this_fname, args.fname_prefix, prec=150, FigFormat=args.FigFormat)
            TerracePlotter.MakeTerraceHeatMapNormalised(this_dir,this_fname, args.fname_prefix, prec=150, FigFormat=args.FigFormat)

#=============================================================================
# This is just a welcome screen that is displayed if no arguments are provided.
#=============================================================================
def print_welcome():

    print("\n\n=======================================================================")
    print("Hello! Welcome to the terrace long profiler tool.")
    print("You will need to tell me which directory to look in.")
    print("Use the -dir flag to define the working directory.")
    print("If you don't do this I will assume the data is in the same directory as this script.")
    print("For help type:")
    print("   python terrace_profile_plots.py -h\n")
    print("=======================================================================\n\n ")

#=============================================================================
if __name__ == "__main__":
    main(sys.argv[1:])
