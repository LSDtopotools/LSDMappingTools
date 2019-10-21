

# set backend to run on server
import matplotlib
matplotlib.use('Agg')


import geopandas as gpd
import rasterio
from rasterio import features


#from __future__ import print_function
import sys
import os


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
    parser.add_argument("-sfname", "--shapefile_fname", type=str, help="The prefix of your shapefile WITHOUT EXTENSION!!! This must be supplied or you will get an error (unless you're running the parallel plotting).")   
    parser.add_argument("-tfname", "--template_fname", type=str, help="The prefix of your DEM WITH EXTENSION!!! This must be supplied or you will get an error (unless you're running the parallel plotting).") 
    parser.add_argument("-fd", "--field", type=str, help="The name of the field")     
 
    # This is for specifying raster values
    parser.add_argument("-lt", "--use_lookup_table", type=bool, default='false', help="If true, uses a lookup table to set the values of the raster")  
    parser.add_argument("-ltfname", "--lookup_table_fname", type=str, help="The prefix of your lookup table WITH EXTENSION!!!")
    parser.add_argument("-lfd", "--lookup_field", type=str, help="The name of the lookup field")  

    args = parser.parse_args()

    print(argv)
    print(args)

    # get the base directory
    if args.base_directory:
        this_dir = args.base_directory
        # check if you remembered a / at the end of your path_name
        if not this_dir.endswith("/"):
            print("You forgot the '/' at the end of the directory, appending...")
            this_dir = this_dir+"/"
    else:
        this_dir = os.getcwd()
        
    field = args.field    
    

    shp_fn = this_dir+args.shapefile_fname+".shp"
    temp_shp_fn = this_dir+args.shapefile_fname+"_temp_shapefile.shp"
    rst_fn = this_dir+args.template_fname
    out_fn = this_dir+args.fname_prefix+".tif"
    
    
    print("The files being used are:")
    print(shp_fn)
    print(rst_fn)
    print(out_fn)

    # read the shapefle
    tgdf = gpd.read_file(shp_fn)   
    
    if args.lookup_table_fname:
        print("I am going to proceed using a lookup table")
        
        # Load the dataframe
        lookup_fn = this_dir+args.lookup_table_fname
        
 
        # creating the rasterization column
        tgdf['rasterization_equivalent'] = '' 
        
        lookup_df = gpd.read_file(lookup_fn) 
        
        keys = df.values[field]
        values = float(df.values[args.lookup_field])
        
        dictionary = dict(zip(lookup, values))
        
        print("The dictionary is:")
        print(dictionary)
        
        # Loop through the unique values
        for key in keys:
            print("The key is: "+key)
            tgdf['rasterization_equivalent'][tgdf[field] == key] = dictionary[key]
            raster_val.append(dictionary[key])
            shp_val.append(key)     

    else:
        
        print("I am going to proceed using unique values")
        print("This will relsult in a new shapefile")
        # creating the rasterization column
        tgdf['rasterization_equivalent'] = '' 

        #Hosts the equivalences
        eqdic = {}

        # Fill the null values
        tgdf[field].fillna(-5, inplace=True)

        uniquevals = tgdf[field].unique()
        print("The unique values are")
        print(uniquevals)

        # Try to get unique values from the fields
        incre = 0
        print('Creating the correspondances string - values')
        raster_val = []
        shp_val = []

        # Loop through the unique values
        for val in uniquevals:
            print(val)
            eqdic[incre] = val
            tgdf['rasterization_equivalent'][tgdf[field] == val] = incre
            raster_val.append(incre)
            shp_val.append(val)
            incre += 1
        eqdic[-9999] = 'NoData'


    # Get the template raster
    rst = rasterio.open(rst_fn)
    
    # copy and update the metadata from the input raster for the output
    meta = rst.meta.copy()
    meta.update(compress='lzw')
    
    with rasterio.open(out_fn, 'w+', **meta) as out:
        out_arr = out.read(1)

        # this is where we create a generator of geom, value pairs to use in rasterizing
        shapes = ((geom,value) for geom, value in zip(tgdf.geometry, tgdf.rasterization_equivalent))

        burned = features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
        out.write_band(1, burned)

#=============================================================================
if __name__ == "__main__":
    main(sys.argv[1:])