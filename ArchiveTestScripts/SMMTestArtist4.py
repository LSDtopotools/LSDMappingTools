# -*- coding: utf-8 -*-
"""
Created on Mon Jun 05 12:28:29 2017

@author: smudd
"""


import matplotlib

# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')


import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rcParams
from LSDMapFigure.PlottingRaster import MapFigure


#label_size = 100
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
#rcParams['font.size'] = label_size
#rcParams['lines.linewidth']  = 1.5   


DataDirectory = "/home/smudd/SMMDataStore/analysis_for_papers/movern_testing/"
Base_file = "Irian_Jaya_PP"

#Directory = "/home/s1563094/Datastore/DATA/UK/LiDAR_DTM_1m/HIN/"
#Base_file = "HIN_"

#BackgroundRasterName = Base_file+".bil"
BackgroundRasterName = Base_file+".bil"
DrapeRasterName = Base_file+"_hs.bil"
BasinRasterName = Base_file+"_AllBasins.bil"
DischargeRasterName= Base_file+"_Q.bil"



plt.clf() 
cbar_loc = "right"
MF = MapFigure(BackgroundRasterName, DataDirectory,coord_type="UTM_km",colourbar_location = cbar_loc)
#MF.add_drape_image(BackgroundRasterName,DataDirectory,alpha = 1)
MF.add_drape_image(DischargeRasterName,DataDirectory,colourmap = "cubehelix", alpha = 0.6)
ImageName = DataDirectory+"TestNewArtist.png" 
fig_size_inches = 12
ax_style = "Normal"
MF.save_fig(fig_width_inches = fig_size_inches, FigFileName = ImageName, axis_style = ax_style, Fig_dpi = 250)
