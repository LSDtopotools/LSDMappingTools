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
from __future__ import division
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
import random


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

	generate_legend_in_csv(dict_file["lithokey"][0]+dict_file["lithokey"][1])
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


def generate_legend_in_csv(name, force = False):
	"""
	This function will generate a random color legend in your lithokey file by assigning a random HTML color hex.
	The idea is to easily vizualise geology by color first, then to be able to directly and manually change the color in the csv file.
	This will generate legend only if this hasn't been done yet to avoid overwritting you existing colors.
	param:
		name (str): name of the csv lithokey file. It should be prefix_lithokey.csv
		force (bool): 
	Do someone read the documentation? if yes just send me an email.

	Author: BG
	Date: 21/09/2017
	"""

	dft = pd.read_csv(name)
	if ("legend" not in dft.columns) or (force == True) :
		print("I am generating random html hex color code for each geology/lithology/whatever you had in your map")
		print("you can manually change it for the next plotting in this file:")
		print(name)
		print("To get a visual version of the color html hex code google it")


		dft["legend"]= pd.Series(np.zeros(dft.shape[0]), index=dft.index)
		for i in range(dft.shape[0]):
			colorTemp = "#"
			for j in range(6):
				colorTemp = colorTemp + (random.choice('0123456789ABCDEF'))
			dft["legend"][i] = colorTemp
		dft.to_csv(name, index = False)
	else:
		print("You already have a legend in your lithokey file! Well done.")

	#now creating and returning the colormap 
	cm = colors.LinearSegmentedColormap.from_list("LithoColorMap", dft.legend, N=dft.shape[0]-1)
	return cm


def movern_two_litho(fname_prefix, DataDirectory, lname_prefix = "", only_MLE = False,SA_channels = False, normalization = True, litho = [], show_legend = True, size_format = "ESURF", basin_list=[], start_movern=0.2, d_movern=0.1, n_movern=7):
	"""
		Function to plot the movern repartition beetween two lithologies per basins, using the percentage of each lithologies in the basin.
		something like this:
		

				|	o	  o o  |
				| o			  |
		   m/n  |	o			o|
				|				|
				|_o______________|
			  lith1			lith2
			  100%			  100%  

		@params:
			fname_prefix (str): the prefix of all your files.
			DataDirectory (str): path to your file
			lname_prefix (str): prefix of your lithology files, leave default "" if you do not know what you are doing.  it should be automatic.
			normalization (bool): Do you want to normalize to 100 %
			litho (list of int/str): list of lithologic code or name [1,54] or ["granit", "anorthosite"]
			show_legend (bool): ...
		@return:
			nothing but create a plot
		@Authors: BG
		@Date: 25/09/2017
		@State: Work in progress, Everything still need to be done.

	"""
	print("Printing a recapitulative chart containing the movern value in function of the basin lithology")

	
	# check if a directory exists for the summary plots. If not then make it.
	summary_directory = DataDirectory+'summary_plots/'
	if not os.path.isdir(summary_directory):
		os.makedirs(summary_directory)
	if(len(litho) != 2):
		print("Sorry mate, but I need two lithogies to plot this. Not 1, not 4, not %s, just 2" %(len(litho)))
		quit()


	# make a figure
	if size_format == "geomorphology":
		fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
		#l_pad = -40
	elif size_format == "big":
		fig = plt.figure(1, facecolor='white',figsize=(16,9))
		#l_pad = -50
	else:
		fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
		#l_pad = -35


	if show_legend:
		gs = plt.GridSpec(100,100,bottom=0.15,left=0.05,right=0.75,top=0.9)
	else:
		gs = plt.GridSpec(100,100,bottom=0.1,left=0.05,right=0.95,top=0.95)

	ax = fig.add_subplot(gs[5:100,10:95])

	# read in the summary csv
	df = Helper.ReadMOverNSummaryCSV(summary_directory,fname_prefix)

	if basin_list != []:
		basin_keys = basin_list
	else:
		# get the basin keys
		basin_keys = df['basin_key'].tolist()
		print basin_keys

	df = df[df['basin_key'].isin(basin_keys)]

	#read the lithologic summary
	try:
		dfl = pd.read_csv(DataDirectory+fname_prefix+"_SBASLITH.csv", sep=",")
		dfl = dfl[dfl["method"] == "percentage"]
	except IOError:
		print("You don't have the requested lithology file that summarize the basin informations") # TODO: Explain this error
		quit()
	# Add the lithologic informations
	df["litho_percent"] = pd.Series(np.zeros(df.shape[0]), index=df.index)
	# getting the litho infos
	for i in df["basin_key"].unique():
		first_percentage = dfl[str(litho[0])][dfl["basin_id"]==i].values
		second_percentage = dfl[str(litho[1])][dfl["basin_id"]==i].values
		tot = first_percentage+second_percentage
		
		# Normalization to 100 %
		if( tot != 100 and normalization):
			first_percentage = first_percentage*100/tot
			second_percentage = second_percentage*100/tot
		# setting the value 
		df["litho_percent"][df["basin_key"]==i] = second_percentage




	#print df
	#ax.scatter(df["litho_percent"],df["Chi_MLE_full"],c = df["basin_key"])
	if (only_MLE):
		# plot the points data
		median_movern = df['Chi_MLE_points'].as_matrix()
		points_max_err = df['Chi_MLE_points_max'].as_matrix()
		points_max_err = points_max_err-median_movern
		points_min_err = df['Chi_MLE_points_min'].as_matrix()
		points_min_err = median_movern-points_min_err
		errors = np.array(zip(points_min_err, points_max_err)).T

		points_chi_keys = df['basin_key'].as_matrix()-0.1
		ax.errorbar(df["litho_percent"], df['Chi_MLE_points'], s=15, marker='o', xerr=None, yerr=errors, ecolor="#fdbb84", fmt='none', elinewidth=1,label='_nolegend_')
		ax.scatter(df["litho_percent"], df['Chi_MLE_points'], s=15, c=df["basin_key"], marker='o', edgecolors='k', lw=0.5,facecolors='#fdbb84', label='Chi Monte Carlo',zorder=200)
	else:
		# plot the full chi data
		full_chi_keys = df['basin_key'].as_matrix()-0.2
		ax.scatter(df["litho_percent"], df['Chi_MLE_full'],marker='o', edgecolors='k', lw=0.5, facecolors='#e34a33', s=15, zorder=200, label='Chi all data')

		# plot the points data
		median_movern = df['Chi_MLE_points'].as_matrix()
		points_max_err = df['Chi_MLE_points_max'].as_matrix()
		points_max_err = points_max_err-median_movern
		points_min_err = df['Chi_MLE_points_min'].as_matrix()
		points_min_err = median_movern-points_min_err
		errors = np.array(zip(points_min_err, points_max_err)).T

		points_chi_keys = df['basin_key'].as_matrix()-0.1
		ax.errorbar(df["litho_percent"], df['Chi_MLE_points'], s=15, marker='x', xerr=None, yerr=errors, ecolor='#fdbb84', fmt='none', elinewidth=1,label='_nolegend_')
		ax.scatter(df["litho_percent"], df['Chi_MLE_points'], s=15, c='#fdbb84', marker='x', edgecolors='k', lw=0.5,facecolors='#fdbb84', label='Chi Monte Carlo',zorder=200)

		# plot the SA data
		SA_keys = df['basin_key'].as_matrix()
		SA_sterr = df['SA_raw_sterr'].as_matrix()
		ax.errorbar(df["litho_percent"], df['SA_raw'], yerr=SA_sterr, c='#2b8cbe', elinewidth=1, fmt='none',label='_nolegend_')
		ax.scatter(df["litho_percent"], df['SA_raw'], s=15, c='#2b8cbe', edgecolors='k',lw=0.5, label='S-A all data', zorder=100)

		if SA_channels:
			# plot the SA data by tribs
			median_movern = df['SA_tribs'].as_matrix()
			points_max_err = df['SA_tribs_max'].as_matrix()
			points_max_err = points_max_err-median_movern
			points_min_err = df['SA_tribs_min'].as_matrix()
			points_min_err = median_movern-points_min_err
			errors = np.array(zip(points_min_err, points_max_err)).T

			SA_tribs_keys = df['basin_key'].as_matrix()+0.1
			ax.errorbar(df["litho_percent"], df['SA_tribs'], s=15, marker='v', facecolors='white', xerr=None, yerr=errors, edgecolors='r', fmt='none', elinewidth=1, linestyle = ":", ecolor='r',label='_nolegend_')
			ax.scatter(df["litho_percent"], df['SA_tribs'], s=15, marker='v', facecolors='white', edgecolors='r', label='S-A by channel',zorder=100)

		# plot the segmented SA data
		median_movern = df['SA_segments'].as_matrix()
		points_max_err = df['SA_segments_max'].as_matrix()
		points_max_err = points_max_err-median_movern
		points_min_err = df['SA_segments_min'].as_matrix()
		points_min_err = median_movern-points_min_err
		errors = np.array(zip(points_min_err, points_max_err)).T

		SA_segment_keys = df['basin_key'].as_matrix()+0.2
		ax.errorbar(df["litho_percent"], df['SA_segments'], s=15, marker='^', facecolors='#a6bddb', xerr=None, yerr=errors, edgecolors='#a6bddb', fmt='none', elinewidth=1, linestyle = ":", ecolor='#a6bddb',label='_nolegend_')
		ax.scatter(df["litho_percent"], df['SA_segments'], s=15, marker='^', facecolors='#a6bddb', edgecolors='k', lw=0.5, label='Segmented S-A', zorder=100)



	# set the axis labels
	ax.set_xlabel('%  Lithologies')
	ax.set_ylabel('Best fit $m/n$')

	if show_legend:
		print "ADDING THE LEGEND"
		# sort both labels and handles by labels
		handles, labels = ax.get_legend_handles_labels()
		labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
		# add the legend
		ax.legend(handles, labels,fontsize=8, bbox_to_anchor=(1.0,0.7),bbox_transform=plt.gcf().transFigure)

	# This gets all the ticks, and pads them away from the axis so that the corners don't overlap
	#ax.tick_params(axis='both', width=1, pad = 2)
	#for tick in ax.xaxis.get_major_ticks():
#		tick.set_pad(2)

	# change tick spacing
	#ax.xaxis.set_major_locator(ticker.MultipleLocator(base=1))

	if(only_MLE):
		plt.savefig(summary_directory+fname_prefix+"_movern_two_litho_MLE.png", dpi = 300)

	else:
		plt.savefig(summary_directory+fname_prefix+"_movern_two_litho.png", dpi = 300)
	
	plt.clf()

	
