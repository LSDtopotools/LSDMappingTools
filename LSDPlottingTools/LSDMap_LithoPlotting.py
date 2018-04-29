#=============================================================================
# These functions create figures for visualising the m/n selection using the chi
# mapping tool.
#
# It creates separate plots for each basin in the DEM.
#
# Authors:
#	 Simon M. Mudd
#	 Fiona J. Clubb
#	 Boris Gailleton
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
from . import LSDMap_MOverNPlotting as MN
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
				print(result)
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


def getLithoColorMap(fname_prefix, DataDirectory, values = "None"):
	"""
		This function get the custom colormap from your lithokey file.
		This has been generated automatically if everything is ok.
		@param:
			fname_prefix (str): the common prefix of all your files
			DataDirectory (str): Path to your file with / at the end
		@return colormap
		@Author: BG
		@Date: 26/09/2017
	"""
	df = pd.read_csv(DataDirectory+fname_prefix+"_lithokey.csv")
	if("legend" not in df.columns.values):
		print("please check your files first with -c True, I cannae find your legend column")
		print("An alternative solution is to add a legend column in your lithokey csv file")
		print("with an html color code for each lithologies, google html color code generator to obtain one")
		quit()
	if(values == "None"):
		colorList = []
		for i in range(df["rocktype"].max()):
			if(i in df["rocktype"].values):
				colorList.append(df["legend"][df["rocktype"] == i].values[0])
			else:
				colorList.append("#FFFFFF")



		#Creating the colormaps
		cm = colors.LinearSegmentedColormap.from_list("LithoColorMap", colorList, N=len(colorList))
	else:
		if values.shape[0]>0:
			colorList = []
			values_np = np.array(values)
			for i in range(values.min(),values.max()+1):
				if(i in df["rocktype"].values):
					colorList.append(df["legend"][df["rocktype"] == i].values[0])
				else:
					colorList.append("#01FA1E")


			#Creating the colormaps
			if len(colorList) == 1:
				cm = colorList[0]
			else:

				cm = colors.LinearSegmentedColormap.from_list("LithoColorMap", colorList, N=len(colorList))
		else:
			print("Something went wrong when I tried to get your colors. try to check your lithokey file.")
			quit()

	return cm

def get_color_litho(fname_prefix, DataDirectory, lithocode):
	"""
		return the color code of a single litho
		@param:
			DataDirectory (str): the data directory with the m/n csv files
			fname_prefix (str): The prefix for the m/n csv files
			litho: code of the lithology

		@Returns:
			an html color code probably for plotting purpose.
		@Authors:
			BG
		@Date: yes

	"""

	df = pd.read_csv(DataDirectory+fname_prefix+"_lithokey.csv")
	if("legend" not in df.columns.values):
		print("please check your files first with -c True, I cannae find your legend column")
		print("An alternative solution is to add a legend column in your lithokey csv file")
		print("with an html color code for each lithologies, google html color code generator to obtain one")
		quit()

	cocode = df["legend"][df["rocktype"] == lithocode].values[0]
	return cocode


def MakeRasterLithoBasinMap(DataDirectory, fname_prefix, lname_prefix, lithodict, size_format='ESURF', FigFormat='png', basins = True, m_chi = False, mancol = [], log_scale_river = False, minmax_m_chi = []):
	"""
	This function makes a shaded relief plot of the DEM with lithologic map on the top and basin outline

	Args:
	DataDirectory (str): the data directory with the m/n csv files
	fname_prefix (str): The prefix for the m/n csv files
	size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
	FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
	Basins (bool): Do you want the basin on top
	minmax_m_chi (list): define a minimum/maximum for plotting m_chi on the top of litho (at the moment this plot is generated from knickpoint dataset)

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

	raster_ext = '.bil'

	# get the basin IDs to make a discrete colourmap for each ID
	if(basins):
		BasinInfoDF = Helper.ReadBasinInfoCSV(DataDirectory, fname_prefix)

		basin_keys = list(BasinInfoDF['basin_key'])
		basin_keys = [int(x) for x in basin_keys]

		basin_junctions = list(BasinInfoDF['outlet_junction'])
		basin_junctions = [float(x) for x in basin_junctions]

		print ('Basin keys are: ')
		print (basin_keys)
		BasinsName = fname_prefix+'_AllBasins.bil'



	# going to make the basin plots - need to have bil extensions.
	print("I'm going to make the basin plots. Your topographic data must be in ENVI bil format or I'll break!!")

	# get the rasters
	raster_ext = '.bil'
	BackgroundRasterName = fname_prefix+raster_ext
	HillshadeName = fname_prefix+'_hs'+raster_ext

	LithoMap = lname_prefix+raster_ext

	# create the map figure
	MF = MapFigure(HillshadeName, DataDirectory,coord_type="UTM_km", colourbar_location='None')
	MF.add_drape_image(HillshadeName,DataDirectory,NFF_opti = True, custom_min_max = [90,240], alpha = 1)
	# add the geology drape
	# MF.add_basin_plot(BasinsName,fname_prefix,DataDirectory,
	#				  use_keys_not_junctions = True, show_colourbar = True,
	#				  discrete_cmap=True, n_colours=len(basin_keys), colorbarlabel = "Basin ID",
	#				  colourmap = cmap, adjust_text = False)
	# getting the right color now

	color_map_litho  = getLithoColorMap(fname_prefix, DataDirectory)
	df_litho_size = pd.read_csv(DataDirectory+fname_prefix+"_lithokey.csv")

	MF.add_drape_image(LithoMap,DataDirectory,colourmap = color_map_litho,
						alpha=0.6,
						show_colourbar = False,
						colorbarlabel = "Colourbar", discrete_cmap=False,
						norm = "None",
						colour_min_max = [0,df_litho_size["rocktype"].max()-1],
						modify_raster_values=False,
						old_values=[], new_values=[], cbar_type=int,
						NFF_opti = True, custom_min_max = [])

	if(basins):
	# add the basin outlines

		Basins = LSDP.GetBasinOutlines(DataDirectory, BasinsName)
		MF.plot_polygon_outlines(Basins, linewidth=0.8)

		# knickpoints!
		if(m_chi):
			ChannelDF = pd.read_csv(DataDirectory+fname_prefix+"_ksnkp_mchi.csv")
			ChannelPoints = LSDP.LSDMap_PointData(ChannelDF, data_type = "pandas", PANDEX = True)
			MF.add_point_data(ChannelPoints,column_for_plotting = 'm_chi',show_colourbar = True, scale_points=True, column_for_scaling='drainage_area',alpha=0.5,zorder=100,this_colourmap = "RdBu_r" ,colour_manual_scale = mancol, scaled_data_in_log = log_scale_river,max_point_size = minmax_m_chi[1], min_point_size = minmax_m_chi[0])
		else:
		# add the channel network
			ChannelDF = Helper.ReadChiDataMapCSV(DataDirectory,fname_prefix)
			ChannelPoints = LSDP.LSDMap_PointData(ChannelDF, data_type = "pandas", PANDEX = True)
			MF.add_point_data(ChannelPoints,show_colourbar="False", scale_points=True, column_for_scaling='drainage_area',alpha=0.5,zorder=100)

	if(basins):
		# add the basin labelling
		label_dict = dict(zip(basin_junctions,basin_keys))
		Points = LSDP.GetPointWithinBasins(DataDirectory, BasinsName)
		MF.add_text_annotation_from_shapely_points(Points, text_colour='k', label_dict=label_dict,zorder=200)

	if(basins):
		# Save the figure
		ImageName = raster_directory+fname_prefix+'_basin_keys_litho.'+FigFormat
	else:
		ImageName = raster_directory+fname_prefix+'_litho.'+FigFormat

	MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 500)


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


def movern_two_litho(fname_prefix, DataDirectory, lname_prefix ='' , color_by_basin = False, only_MLE = True, SA_channels = False, normalization = True, litho = [], show_legend = True, size_format = "ESURF", basin_list=[], start_movern=0.2, d_movern=0.1, n_movern=7, coloured_by_lith = True):
	"""
		Function to plot the movern repartition beetween two lithologies per basins, using the percentage of each lithologies in the basin.
		something like this:


				|	o	  o o    |
				| o			     |
		   m/n  |	o			o|
				|				 |
				|_o______________|
			  lith1			  lith2
			  100%			   100%

		@params:
			fname_prefix (str): the prefix of all your files.
			DataDirectory (str): path to your file
			lname_prefix (str): prefix of your lithology files, leave default "" if you do not know what you are doing.  it should be automatic.
			normalization (bool): Do you want to normalize to 100 %
			litho (list of int/str): list of lithologic code or name [1,54] or ["granit", "anorthosite"]
			show_legend (bool): ...
			coloured by litho (bool): Points will be coloured by lithologic content.
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
	else:
		if isinstance(litho[0],str) and isinstance(litho[1], str):
			#converting the litho string to litho code by reading the lithokey file
			correspondance_pandas_DataFrame = pd.read_csv(DataDirectory+fname_prefix+"_lithokey.csv", sep =",")
			litho[0] = correspondance_pandas_DataFrame["rocktype"][correspondance_pandas_DataFrame["ID"]==litho[0]].values[0] # ID - rocktype - legend
			litho[1] = correspondance_pandas_DataFrame["rocktype"][correspondance_pandas_DataFrame["ID"]==litho[1]].values[0] # ID - rocktype - legend



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


	if show_legend and not only_MLE:
		gs = plt.GridSpec(100,100,bottom=0.15,left=0.05,right=0.75,top=0.9)
	else:
		if (only_MLE):
			gs = plt.GridSpec(100,100,bottom=0.15,left=0.05,right=0.9,top=0.9)
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
		print (basin_keys)

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
		# plot the data with only the MLE, neater plot, I am condisering only keeping this one tbh
		# plot the points data
		median_movern = df['Chi_MLE_points'].as_matrix()
		points_max_err = df['Chi_MLE_points_max'].as_matrix()
		points_max_err = points_max_err-median_movern
		points_min_err = df['Chi_MLE_points_min'].as_matrix()
		points_min_err = median_movern-points_min_err
		errors = np.array(zip(points_min_err, points_max_err)).T


		points_chi_keys = df['basin_key'].as_matrix()-0.1

		# generating Random color for basins
		# Assigning the random color values
		if color_by_basin:
			colormlape = []
			for i in range(df.shape[0]):
				colorTemp = "#"
				for j in range(6):
					colorTemp = colorTemp + (random.choice('0123456789ABCDEF'))
				colormlape.append(colorTemp)

			Cmpalnama = colors.LinearSegmentedColormap.from_list("Basimap", colormlape, N=df.shape[0]-1)
			tpp = ax.scatter(df["litho_percent"], df['Chi_MLE_points'], s=15, c=df["basin_key"], cmap = Cmpalnama, marker='o', edgecolors='k', lw=0.5,facecolors='#fdbb84', label='Chi Monte Carlo',zorder=200)
			ax.errorbar(df["litho_percent"], df['Chi_MLE_points'], s=15, marker='o', xerr=None, yerr=errors, ecolor="#fdbb84", fmt='none', elinewidth=1,label='_nolegend_')
		else:
			if(coloured_by_lith):
				df_lith1 = df[df["litho_percent"]<=50]
				median_movern = df_lith1['Chi_MLE_points'].as_matrix()
				points_max_err = df_lith1['Chi_MLE_points_max'].as_matrix()
				points_max_err = points_max_err-median_movern
				points_min_err = df_lith1['Chi_MLE_points_min'].as_matrix()
				points_min_err = median_movern-points_min_err
				errors = np.array(zip(points_min_err, points_max_err)).T
				tpp = ax.scatter(df_lith1["litho_percent"], df_lith1['Chi_MLE_points'], s=15, c=get_color_litho(fname_prefix,DataDirectory,litho[0]),  marker='o', edgecolors='k', lw=0.5,facecolors=get_color_litho(fname_prefix,DataDirectory,litho[0]), label='Chi Monte Carlo',zorder=200)
				ax.errorbar(df_lith1["litho_percent"], df_lith1['Chi_MLE_points'], s=15, marker='o', xerr=None, yerr=errors, ecolor=get_color_litho(fname_prefix,DataDirectory,litho[0]), fmt='none', elinewidth=1,label='_nolegend_')

				df_lith2 = df[df["litho_percent"]>50]
				median_movern = df_lith2['Chi_MLE_points'].as_matrix()
				points_max_err = df_lith2['Chi_MLE_points_max'].as_matrix()
				points_max_err = points_max_err-median_movern
				points_min_err = df_lith2['Chi_MLE_points_min'].as_matrix()
				points_min_err = median_movern-points_min_err
				errors = np.array(zip(points_min_err, points_max_err)).T
				tpp = ax.scatter(df_lith2["litho_percent"], df_lith2['Chi_MLE_points'], s=15, c=get_color_litho(fname_prefix,DataDirectory,litho[1]),  marker='o', edgecolors='k', lw=0.5,facecolors=get_color_litho(fname_prefix,DataDirectory,litho[1]), label='Chi Monte Carlo',zorder=200)
				ax.errorbar(df_lith2["litho_percent"], df_lith2['Chi_MLE_points'], s=15, marker='o', xerr=None, yerr=errors, ecolor=get_color_litho(fname_prefix,DataDirectory,litho[1]), fmt='none', elinewidth=1,label='_nolegend_')

				# for q in range(df.shape[0]):

				# 	if df["litho_percent"].values[q]<=50:
				# 		coloritemp = get_color_litho(fname_prefix, DataDirectory, litho[0])
				# 	else:
				# 		coloritemp = get_color_litho(fname_prefix, DataDirectory, litho[1])
				# 	#print coloritemp

				# 	tpp = ax.scatter(df["litho_percent"].values[q], df['Chi_MLE_points'].values[q], s=15, c=coloritemp,  marker='o', edgecolors='k', lw=0.5,facecolors=coloritemp, label='Chi Monte Carlo',zorder=200)
				# 	ax.errorbar(df["litho_percent"].values[q], df['Chi_MLE_points'].values[q], s=15, marker='o', xerr=None, yerr=errors[0][q], ecolor=coloritemp, fmt='none', elinewidth=1,label='_nolegend_')
			else:
				tpp = ax.scatter(df["litho_percent"], df['Chi_MLE_points'], s=15, c="#fdbb84",  marker='o', edgecolors='k', lw=0.5,facecolors='#fdbb84', label='Chi Monte Carlo',zorder=200)
				ax.errorbar(df["litho_percent"], df['Chi_MLE_points'], s=15, marker='o', xerr=None, yerr=errors, ecolor="#fdbb84", fmt='none', elinewidth=1,label='_nolegend_')

		if(color_by_basin and show_legend):
			gs2 = plt.GridSpec(100,100,bottom=0,left=0,right=1,top=1)
			cax = fig.add_subplot(gs2[10:90,90:93])

			plt.colorbar(tpp, cax = cax,ticks =[0,df.basin_key.max()], orientation = "vertical",label = '')
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
	ax.set_ylabel('Best fit '+r'$\theta$')

	#Now deal with the Lithology names

	correspondance_pandas_DataFrame = pd.read_csv(DataDirectory+fname_prefix+"_lithokey.csv", sep =",")
	litho_wanne = correspondance_pandas_DataFrame["ID"][correspondance_pandas_DataFrame["rocktype"]==litho[0]].values[0] # ID - rocktype - legend
	litho_tou = correspondance_pandas_DataFrame["ID"][correspondance_pandas_DataFrame["rocktype"]==litho[1]].values[0] # ID - rocktype - legend

	ax.set_xticks([0,25,50,75,100])
	ax.set_xticklabels([litho_wanne,'',"50",'',litho_tou])

	if show_legend and not only_MLE:
		print ("ADDING THE LEGEND")
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

def MakeChiPlotsByLith(DataDirectory, fname_prefix, basin_list=[0], start_movern=0.2, d_movern=0.1, n_movern=7,
					size_format='ESURF', FigFormat='png', animate=False, keep_pngs=False, plot_colorbar = False):
	"""
	This function makes chi-elevation plots for each basin and each value of m/n
	where the channels are coloured by the K value (for model runs with spatially varying K).

	Args:
		DataDirectory (str): the data directory with the m/n csv files
		fname_prefix (str): The prefix for the m/n csv files
		basin_list: a list of the basins to make the plots for. If an empty list is passed then
		all the basins will be analysed. Default = basin 0.
		start_movern (float): the starting m/n value. Default is 0.2
		d_movern (float): the increment between the m/n values. Default is 0.1
		n_movern (float): the number of m/n values analysed. Default is 7.
		size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
		FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
		animate (bool): If this is true then it creates a movie of the chi-elevation plots coloured by MLE.
		keep_pngs (bool): If this is false and the animation flag is true, then the pngs are deleted and just the video is kept.

	Returns:
		Plot of each m/n value for each basin.

	Author: FJC
	"""
	from LSDPlottingTools import colours

	# check if a directory exists for the chi plots. If not then make it.
	K_directory = DataDirectory+'chi_plots_Lith/'
	if not os.path.isdir(K_directory):
		os.makedirs(K_directory)

	# Set up fonts for plots
	label_size = 10
	rcParams['font.family'] = 'sans-serif'
	rcParams['font.sans-serif'] = ['arial']
	rcParams['font.size'] = label_size

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


	#colorbar axis
	gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.95,top=1.0)
	if plot_colorbar:
		ax2 = fig.add_subplot(gs[10:95,82:85])
		#gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.95,top=1.0)
		ax = fig.add_subplot(gs[10:95,5:80])
	else:
		ax = fig.add_subplot(gs[5:95,5:95])
	# read in the csv files
	ProfileDF = Helper.ReadChiProfileCSV(DataDirectory, fname_prefix)
	BasinStatsDF = Helper.ReadBasinStatsCSV(DataDirectory, fname_prefix)

	# get the number of basins
	basin_keys = list(BasinStatsDF['basin_key'])
	basin_keys = [int(x) for x in basin_keys]

	# get the list of basins
	if basin_list == []:
		print("You didn't give me a list of basins, so I'll just run the analysis on all of them!")
		basin_list = basin_keys

	# loop through each m over n value
	end_movern = start_movern+d_movern*(n_movern-1)
	m_over_n_values = np.linspace(start_movern,end_movern,n_movern)

	# best fit moverns
	best_fit_moverns = MN.SimpleMaxMLECheck(BasinStatsDF)

	for m_over_n in m_over_n_values:
		# read in the full stats file
		print("This m/n is: "+str(m_over_n))
		FullStatsDF = Helper.ReadFullStatsCSV(DataDirectory,fname_prefix,m_over_n)

		# loop through all the basins in the basin list
		for basin_key in basin_list:
			print("This basin key is %s") %str(basin_key)

			# mask the data frames for this basin
			ProfileDF_basin = ProfileDF[ProfileDF['basin_key'] == basin_key]
			FullStatsDF_basin = FullStatsDF[FullStatsDF['basin_key'] == basin_key]

			#print FullStatsDF_basin
			print ("Getting the reference_source_key")

			print(FullStatsDF_basin.iloc[0]['reference_source_key'])

			# get the data frame for the main stem
			ProfileDF_MS = ProfileDF_basin[ProfileDF_basin['source_key'] == FullStatsDF_basin.iloc[0]['reference_source_key']]

			# get the data frame for the tributaries
			ProfileDF_basin = ProfileDF_basin[ProfileDF_basin['source_key'] != FullStatsDF_basin.iloc[0]['reference_source_key']]
			# merge with the full data to get the MLE for the tributaries
			ProfileDF_tribs = ProfileDF_basin.merge(FullStatsDF_basin, left_on = "source_key", right_on = "test_source_key")

			# get the chi and elevation data for the main stem
			movern_key = 'm_over_n = %s' %(str(m_over_n))
			MainStemX = list(ProfileDF_MS[movern_key])
			MainStemElevation = list(ProfileDF_MS['elevation'])
			MainStemK = list(ProfileDF_MS[fname_prefix+"_geol"])


			# get the chi, elevation, and MLE for the tributaries
			TributariesX = list(ProfileDF_tribs[movern_key])
			TributariesElevation = list(ProfileDF_tribs['elevation'])
			TributariesK = list(ProfileDF_tribs[fname_prefix+"_geol"])


			# get the colourmap to colour channels by the MLE value
			#NUM_COLORS = len(MLE)
			K_array = np.asarray(TributariesK)
			#min_K = np.min(K_array)
			#max_K = np.max(K_array)
			#seal_the_seal = pd.concat(TributariesK,MainStemK)

			this_cmap = getLithoColorMap(fname_prefix, DataDirectory, values = ProfileDF_basin[fname_prefix+"_geol"].unique())
			n_colours = 10

			#cNorm  = colors.Normalize(vmin=min_K, vmax=max_K)
			#plt.cm.ScalarMappable(norm=cNorm, cmap=this_cmap)

			# now plot the data with a colourmap
			if(isinstance(this_cmap,str)):
				sc = ax.scatter(TributariesX,TributariesElevation,c=this_cmap, s=2.5, edgecolors='none')
				sc = ax.scatter(MainStemX, MainStemElevation,c=this_cmap, s=2.5, edgecolors='none')
			else:
				sc = ax.scatter(TributariesX,TributariesElevation,c=TributariesK,cmap=this_cmap, s=2.5, edgecolors='none')
				sc = ax.scatter(MainStemX, MainStemElevation,c=MainStemK,cmap=this_cmap, s=2.5, edgecolors='none')

			# some formatting of the figure
			ax.spines['top'].set_linewidth(1)
			ax.spines['left'].set_linewidth(1)
			ax.spines['right'].set_linewidth(1)
			ax.spines['bottom'].set_linewidth(1)

			# make the labels
			ax.set_xlabel("$\chi$ (m)")
			ax.set_ylabel("Elevation (m)")

			# the best fit m/n
			best_fit_movern = best_fit_moverns[basin_key]
			print ("BEST FIT M/N IS: "+ str(best_fit_movern))
			print ("THIS M/N IS: "+str(m_over_n))

			# label with the basin and m/n
			title_string = "Basin "+str(basin_key)+", "+r"$\theta$ = "+str(m_over_n)
			if best_fit_movern == m_over_n:
				ax.text(0.05, 0.95, title_string,
						verticalalignment='top', horizontalalignment='left',
						transform=ax.transAxes,
						color='red', fontsize=10)
			else:
				ax.text(0.05, 0.95, title_string,
						verticalalignment='top', horizontalalignment='left',
						transform=ax.transAxes,
						color='black', fontsize=10)

			# add the colorbar
			if(plot_colorbar):
				colorbarlabel = "$Lith$"
				cbar = plt.colorbar(sc,cmap=this_cmap,spacing='uniform', orientation='vertical',cax=ax2)
				cbar.set_label(colorbarlabel, fontsize=10)
				ax2.set_ylabel(colorbarlabel, fontname='Arial', fontsize=10)

				# #change labels to scientific notation
				# colours.fix_colourbar_ticks(cbar,n_colours, cbar_type=float, min_value = min_K, max_value = max_K, cbar_label_rotation=0, cbar_orientation='vertical')
				# # we need to get linear values between min and max K
				# these_labels = np.linspace(min_K,max_K,n_colours)
				# # now round these and convert to scientific notation
				# these_labels = [str('{:.2e}'.format(float(x))) for x in these_labels]
				# new_labels = []
				# for label in these_labels:
				#	 a,b = label.split("e")
				#	 b = b.replace("0", "")
				#	 new_labels.append(a+' x 10$^{%s}$' % b)

				# ax2.set_yticklabels(new_labels, fontsize=8)

			#save the plot
			newFilename = K_directory+"Chi_profiles_by_Lith_"+str(basin_key)+"_"+str(m_over_n)+"."+str(FigFormat)

			# This gets all the ticks, and pads them away from the axis so that the corners don't overlap
			ax.tick_params(axis='both', width=1, pad = 2)
			for tick in ax.xaxis.get_major_ticks():
				tick.set_pad(2)
			# WTF is there an indent issue there

			plt.savefig(newFilename,format=FigFormat,dpi=300)
			ax.cla()
			if(plot_colorbar):
				ax2.cla()

	if animate:
		# animate the pngs using ffmpeg
		system_call = "ffmpeg -framerate 3 -pattern_type glob -i '"+K_directory+"Chi_profiles_by_Lith*.png' -y -vcodec libx264 -s 1230x566 -pix_fmt yuv420p "+K_directory+"Chi_profiles_by_Lith.mp4"
		print (system_call)
		subprocess.call(system_call, shell=True)
		# delete the pngs if you want
		if not keep_pngs:
			system_call = "rm "+K_directory+"Chi_profiles_by_Lith*.png"
			subprocess.call(system_call, shell=True)
	plt.close(fig)
