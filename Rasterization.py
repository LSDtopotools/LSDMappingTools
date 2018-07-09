#=============================================================================
# Script to plot the litho information acquired from 
#
# Authors:
#	 Boris Gailleton, Simon Marius Mudd
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
import subprocess
from LSDPlottingTools import LSDMap_VectorTools as VT

#=============================================================================
# This is just a welcome screen that is displayed if no arguments are provided.
#=============================================================================
def print_welcome():

	print("\n\n=======================================================================")
	print("Hello! I'm going to raterise a shapefile for you, the aim is to prepare it for LSDTT.")
	print("The LSDTopoTools MuddChi2014 repository can now ingest lithologic data from shapefiles.")
	print("This script creates a csv file that hold the correspondance between geologic ID and Code.")
	print("You will need to tell me the path and the prefix of your shapefile.")
	print("Use the -wd flag to define the working directory.")
	print("If you don't do this I will assume the data is in the same directory as this script.")
	print("A typical basic command would be:")
	print("python Rasterization.py -dir /home/PhilCollins/DataStore/GIS/US/ -fn beverly_litho -fd geology")
	print("A more complete one with more options:")
	print("python Rasterization.py -dir /home/PhilCollins/DataStore/GIS/Africa/ -fn Tanzania_litho -fd LITH -UTM 33 -SOUTH True -rn GeoRast -res 12")
	print("For help type:")
	print("  python Rasterization.py -h\n")
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
	parser.add_argument("-dir", "--directory", type=str, help="The base directory with the shapefile. If this isn't defined I'll assume it's the same as the current directory.")
	parser.add_argument("-fn", "--fname_prefix", type=str, help="The name of your shapefile without extention (no .shp at the end)")
	parser.add_argument("-fd", "--field_name", type=str, default = "", help="The name of the field that will give the value to the raster")
	parser.add_argument("-rast", "--rasterization", type = bool, default = True, help= "Rasterize your file into a tif raster and produce a lithokey csv file containing the keys for the numerisation of the data.")
	parser.add_argument("-rn","--rename", type=str, default = "", help = "rename your raster after rasterization, optional.")
	parser.add_argument("-UTM", "--UTM_CODE", type = str, default = "", help = "Give the UTM UTM_CODE of your zone to automate the conversion of the raster into an LSDTopoTools compatible raster.")
	parser.add_argument("-S", "--SOUTH", type = bool, default = False, help = "Turn to True if your UTM zone is in the southern hemisphere. A quick way to check that is to go on field and to observe how flow the water in a ")
	parser.add_argument("-res", "--resolution", type = int, default = 30, help = "Precise this argument if you want to change the raster resolution, default 30")
	args = parser.parse_args()
#=============================================================================
	rast_name = ""
	rast_pre = ""
	if((not args.directory) or (not args.fname_prefix)):
		print("You need to feed me with a shapefile, try python Rasterization -dir /path/to/your/file/ -fn shapefile.shp")
		sys.exit()

	if(args.rasterization):
		if(args.field_name == ""):
			print("FATAL ERROR: You need to gives me the Field name you wanna rasterize. This field will give its value or a code value to the raster. The field name is in the attributary table of your shapefile.")
			sys.exit()
		#launching the rasterization
		# fd_name_corr = args.field_name.encode('utf-8')
		VT.rasterize_shapefile(args.directory + args.fname_prefix + ".shp", res = args.resolution, field = args.field_name)
		#getting file names and prefix
		rast_name = args.fname_prefix + "_new2.tif"
		rast_pre = args.fname_prefix + "_new2"
		#Renaming in case you want to
		if(args.rename != ""):
			os.rename(args.directory + rast_name,args.directory+args.rename+".tif")
			# updating the file name and prefix
			rast_name = args.rename+".tif"
			rast_pre = args.rename
			
		if(args.UTM_CODE != ""):
			#Conversion to UTM
			if(args.SOUTH):
				south_string = " +south"
			else:
				south_string = ""
			# preparing the gdalwarp command
			gdal_cmd = "gdalwarp -t_srs '+proj=utm +zone=%s%s +datum=WGS84' -of ENVI -tr %s %s %s %s" % (args.UTM_CODE,south_string,args.resolution,args.resolution,args.directory+rast_name, args.directory+rast_pre+'_rast.bil') 
			print("You choose to convert your file, I am generating a gdalwarp command from your information hopefully you won't have errors.")
			print(gdal_cmd)
			subprocess.call(gdal_cmd, shell=True)
			print("Your raster is now ready to use with LSDTopotools, congrats")




if __name__ == "__main__":
	main(sys.argv[1:])