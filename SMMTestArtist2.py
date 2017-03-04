# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 16:57:22 2017

@author: smudd
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 12:55:23 2017

@author: smudd
"""

import matplotlib

# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')


import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rcParams
import matplotlib.cm as cm
from LSDMapFigure.PlottingRaster import MapFigure
from LSDMapFigure.PlottingRaster import BaseRaster

from LSDPlottingTools import colours as lsdcolours
from LSDPlottingTools import init_plotting_DV
from LSDPlottingTools import LSDMap_PointTools
import sys

#sys.path.append("PATH/TO/LSDPlottingTools/")
#
#init_plotting_DV()

#label_size = 100
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
#rcParams['font.size'] = label_size
#rcParams['lines.linewidth']  = 1.5   


#DataDirectory = "/home/smudd/SMMDataStore/analysis_for_papers/Meghalaya/chi_analysis/"
#Directory = "C:\\Vagrantboxes\\LSDTopoTools\\Topographic_projects\\Meghalaya\\Divides\\"
Directory = "C:\\Vagrantboxes\\LSDTopoTools\\Topographic_projects\\Meghalaya\\test_segments\\"
#DataDirectory = "T:\\analysis_for_papers\\Meghalaya/chi_analysis\\"
Base_file = "Mega_divide_segments"

#Directory = "/home/s1563094/Datastore/DATA/UK/LiDAR_DTM_1m/HIN/"
#Base_file = "HIN_"

BackgroundRasterName = Base_file+".bil"
DrapeRasterName = Base_file+"_hs.bil"
ChiRasterName = Base_file+"_chi_coord.bil"


#BR = BaseRaster(BackgroundRasterName, Directory)
#BR.set_raster_type("Terrain")
#print(BR._colourmap)
#BR.show_raster()

#BR.set_colourmap("RdYlGn")
#BR.show_raster()

#PD_file = Base_file+"_chi_coord_basins.csv"  
PD_file = Base_file+"_MChiSegmented.csv"
PointData = LSDMap_PointTools.LSDMap_PointData(Directory+PD_file)

plt.clf() 
cbar_loc = "Bottom"
MF = MapFigure(BackgroundRasterName, Directory,coord_type="UTM_km",colourbar_location = cbar_loc)
MF.add_drape_image(DrapeRasterName,Directory,alpha = 0.4)
#MF.add_drape_image(ChiRasterName,Directory,colourmap = "cubehelix",alpha = 0.4)
MF.add_point_data(PointData,column_for_plotting = "chi",colorbarlabel = "I am point data",
                  scale_points = True,column_for_scaling = "drainage area",
                       scaled_data_in_log = False)
#MF.show_plot()
ImageName = Directory+"TestNewArtist.png" 
fig_size_inches = 6
ax_style = "Normal"
MF.save_fig(fig_width_inches = fig_size_inches, FigFileName = ImageName, axis_style = ax_style, Fig_dpi = 250)



# Customise the DrapePlot
#dp.make_drape_colourbar(cbar_label=colourbar_label)
#dp.set_fig_axis_labels()

#dp.show_plot()