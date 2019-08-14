# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 16:57:22 2017

@author: Yes
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time as clock

# Force matplotlib to not use any Xwindows backend.


from matplotlib import rcParams
import matplotlib.cm as cm
from LSDMapFigure.PlottingRaster import MapFigure
from LSDMapFigure.PlottingRaster import BaseRaster

from LSDPlottingTools import colours as lsdcolours
from LSDPlottingTools import init_plotting_DV
import LSDPlottingTools as LSDP
import sys
#import pandas as bamboo_bears

#sys.path.append("PATH/TO/LSDPlottingTools/")
#
#init_plotting_DV()

label_size = 12
#rcParams['font.family'] = 'sans-serif'
#rcParams['font.sans-serif'] = ['Liberation Sans']
#rcParams['font.size'] = label_size
#rcParams['lines.linewidth']  = 1.5


Directory = "/home/s1675537/PhD/DataStoreBoris/Emma/"
wDirectory = "/home/s1675537/PhD/DataStoreBoris/Emma/"
Base_file = "Betics_UTM30clip_PP"

#df = bamboo_bears.read_csv(rDirectory2 + fname2, sep=",")
csv_file = "/home/s1675537/PhD/DataStoreBoris/Emma/new.csv"
BackgroundRasterName = Base_file + ".bil"
DrapeRasterName = "Betics_UTM30clip_hs.bil"

#BlackRasterD = "/home/s1675537/PhD/DataStoreBoris/GIS/Data/Carpathian/Best_DEM/"
#BlackRaster = "black_35.bil"

alpha_black = 1
#CurveRasterName = Base_file+"_curvature.bil"

#BR = BaseRaster(BackgroundRasterName, Directory)
#BR.set_raster_type("Terrain")
#print(BR._colourmap)
#BR.show_raster()

#BR.set_colourmap("RdYlGn")
#BR.show_raster()
thisPointData = LSDP.LSDMap_PointData(csv_file)
river_network = LSDP.LSDMap_PointData("/home/s1675537/PhD/DataStoreBoris/Emma/rv.csv")

#names = ['cubehelix','CMRmap','RdBu']
#names = ['spring_r', 'autumn_r','YlOrRd','YlOrRd_r']
names = ['autumn_r']
for nami in names:
    plt.clf()
    MF = MapFigure(BackgroundRasterName, Directory,coord_type="UTM_km", alpha = 1)
    MF.add_drape_image(DrapeRasterName,Directory,alpha = 0.5)
    #MF.add_drape_image(BlackRaster, BlackRasterD, alpha = alpha_black)
    MF.add_point_data(river_network,scale_points = True, max_point_size = 1, min_point_size = 0.1, column_for_scaling ="drainage area",  scaled_data_in_log = True)
    MF.add_point_data( thisPointData,column_for_plotting = "elevation",
                       this_colourmap = nami, colorbarlabel = "knickpoint sign", scale_points = True, max_point_size = 100, min_point_size = 0, column_for_scaling ="knickpoints",  scaled_data_in_log = True, minimum_log_scale_cut_off = 2 )


    #MF.add_drape_image(CurveRasterName,Directory,colourmap = "cubehelix",alpha = 0.4, show_colourbar = True, colorbarlabel= "Colourbar")
    #MF.show_plot()
    ImageName = wDirectory+str(int(clock.time()))+nami+".png"
    fig_size_inches = 12
    ax_style = "Normal"
    MF.save_fig(fig_width_inches = fig_size_inches, FigFileName = ImageName, axis_style = ax_style, Fig_dpi = 500)



# Customise the DrapePlot
#dp.make_drape_colourbar(cbar_label=colourbar_label)
#dp.set_fig_axis_labels()

#dp.show_plot()
