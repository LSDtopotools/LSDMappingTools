#=============================================================================
# These functions create figures for visualising the m/n selection using the chi
# mapping tool.
#
# It creates separate plots for each basin in the DEM.
#
# Authors:
#	 Simon M. Mudd
#	 Fiona J. Clubb
#	  Boris Gailleton
#=============================================================================
#=============================================================================
# IMPORT MODULES
#=============================================================================
# set backend to run on server
import matplotlib
matplotlib.use('Agg')

import numpy as np
import LSDPlottingTools as LSDP
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.ticker as ticker
import pandas as pd
from matplotlib import colors
import math
import os
import subprocess
#from shapely.geometry import Polygon
from LSDMapFigure import PlottingHelpers as Helper
from LSDMapFigure.PlottingRaster import MapFigure
from LSDMapFigure.PlottingRaster import BaseRaster
from LSDPlottingTools import LSDMap_SAPlotting as SA
from LSDPlottingTools import joyplot
import fnmatch


def litho_pre_check(directory, lk = "", fname = ""):
	"""
	This function check if you have the right files to work with the litho.
	
	Author: Boris

	"""
	dict_file = {}
	# first checking if there is a lithokey csv file to decrypt the lithologic raster
	if lk == "":
		print("I am checking if you have a lithokey file containing the correspondance between lithology and their keys.")
		result = []
		for root, dirs, files in os.walk(directory):
			for name in files:
				if fnmatch.fnmatch(name, "*lithokey.csv"):
					result.append(name)


		if(len(result) == 1):
			print("I found your lithokey file, this last host the conversion information, just check if this is the right one: ")
			print(result[0])
			
			print("now renaming this file with the right prefix_suffix:")
			system_call = "mv "+directory+result[0]+" "+directory + fname + "_lithokey.csv"
			print(directory + fname + "_lithokey.csv")
			subprocess.call(system_call, shell=True)
			dict_file["lithokey"] = [directory,fname + "_lithokey.csv"]


		else:
			if len(result) >1:
				print("I found several lithokey file. Please specify it with -lk flag cause I cannot choose by myself, give a butcher's by yourself")
				quit()
			else:
				if len(result) == 0:
					print("I cannot find the lithokey file to decrypt the lithologic informations, I am out")
					quit()
		result = []
	else:
		print("You gave me the name of the lithokey file, let me check it")
		dict_file["lithokey"] = [directory,lk]


	# Now reading the litho info
	dict_litho = {}
	df_litho = pd.read_csv(dict_file["lithokey"][0]+dict_file["lithokey"][1])
	for i in range(df_litho.shape[0]):
		dict_litho[df_litho["rocktype"][i]] = df_litho["ID"][i]

	dict_file["lithodict"] = dict_litho
	


	return dict_file

def MakeRasterLithoBasinMap(DataDirectory, fname_prefix, lname_prefix, lithodict, size_format='ESURF', FigFormat='png'):
	"""
	This function makes a shaded relief plot of the DEM with lithologic map on the top and basin outline

	Args:
	DataDirectory (str): the data directory with the m/n csv files
	fname_prefix (str): The prefix for the m/n csv files
	size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
	FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.

	Returns:
	Shaded relief plot with the basins coloured by basin ID

	Author: BG, FJC
	"""
		# check if a directory exists for the chi plots. If not then make it.
	raster_directory = DataDirectory+'raster_plots/'
	if not os.path.isdir(raster_directory):
		os.makedirs(raster_directory)

	# Set up fonts for plots
	label_size = 10
	rcParams['font.family'] = 'sans-serif'
	rcParams['font.sans-serif'] = ['arial']
	rcParams['font.size'] = label_size

	# set figure sizes based on format
	if size_format == "geomorphology":
		fig_width_inches = 6.25
	elif size_format == "big":
		fig_width_inches = 16
	else:
		fig_width_inches = 4.92126

	# get the basin IDs to make a discrete colourmap for each ID
	BasinInfoDF = Helper.ReadBasinInfoCSV(DataDirectory, fname_prefix)

	basin_keys = list(BasinInfoDF['basin_key'])
	basin_keys = [int(x) for x in basin_keys]

	basin_junctions = list(BasinInfoDF['outlet_junction'])
	basin_junctions = [float(x) for x in basin_junctions]

	print ('Basin keys are: ')
	print basin_keys

	# get a discrete colormap
	cmap = plt.cm.jet

	# going to make the basin plots - need to have bil extensions.
	print("I'm going to make the basin plots. Your topographic data must be in ENVI bil format or I'll break!!")

	# get the rasters
	raster_ext = '.bil'
	BackgroundRasterName = fname_prefix+raster_ext
	HillshadeName = fname_prefix+'_hs'+raster_ext
	BasinsName = fname_prefix+'_AllBasins'+raster_ext
	LithoMap = lname_prefix+raster_ext

	# create the map figure
	MF = MapFigure(HillshadeName, DataDirectory,coord_type="UTM_km", colourbar_location='None')

	# add the geology drape
	# MF.add_basin_plot(BasinsName,fname_prefix,DataDirectory,
	#				  use_keys_not_junctions = True, show_colourbar = True,
	#				  discrete_cmap=True, n_colours=len(basin_keys), colorbarlabel = "Basin ID",
	#				  colourmap = cmap, adjust_text = False)
	MF.add_drape_image(LithoMap,DataDirectory,colourmap = plt.cm.jet,
						alpha=0.4,
						show_colourbar = False,
						colorbarlabel = "Colourbar", discrete_cmap=True, n_colours=len(lithodict),
						norm = "None",
						colour_min_max = [],
						modify_raster_values=False,
						old_values=[], new_values=[], cbar_type=int,
						NFF_opti = False, custom_min_max = [])
	# add the basin outlines
	Basins = LSDP.GetBasinOutlines(DataDirectory, BasinsName)
	MF.plot_polygon_outlines(Basins, linewidth=0.8)

	# add the channel network
	ChannelDF = Helper.ReadChiDataMapCSV(DataDirectory,fname_prefix)
	ChannelPoints = LSDP.LSDMap_PointData(ChannelDF, data_type = "pandas", PANDEX = True)
	MF.add_point_data(ChannelPoints,show_colourbar="False", scale_points=True, column_for_scaling='drainage_area',alpha=0.5,zorder=100)

	# add the basin labelling
	label_dict = dict(zip(basin_junctions,basin_keys))
	Points = LSDP.GetPointWithinBasins(DataDirectory, BasinsName)
	MF.add_text_annotation_from_shapely_points(Points, text_colour='k', label_dict=label_dict,zorder=200)

	# Save the figure
	ImageName = raster_directory+fname_prefix+'_basin_keys_litho.'+FigFormat
	MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 300)