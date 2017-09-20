#=============================================================================
# These functions create figures for visualising the m/n selection using the chi
# mapping tool.
#
# It creates separate plots for each basin in the DEM.
#
# Authors:
#     Simon M. Mudd
#     Fiona J. Clubb
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


def litho_pre_check(directory, lk = ""):
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
			dict_file["lithokey"] = [directory,result[0]]
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
	print dict_file 

	return dict_file