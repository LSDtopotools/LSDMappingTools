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
matplotlib.use('Agg')

import matplotlib.cm as cm
from LSDMapFigure.PlottingRaster import MapFigure
from LSDMapFigure.PlottingRaster import BaseRaster
import matplotlib.pyplot as plt
from LSDPlottingTools import colours as lsdcolours
from LSDPlottingTools import init_plotting_DV
import sys

#sys.path.append("PATH/TO/LSDPlottingTools/")
#
#init_plotting_DV()

Directory = "/home/s0923330/Datastore/5m_dems/scotland/DEM_comparison/river_tay/"
#Directory = "T:\\analysis_for_papers\\Meghalaya\\divide_migration\\"
Base_file = "Tay_Nextmap_HS"

BackgroundRasterName = Base_file+".bil"
#DrapeRasterName = Base_file+"_HS.bil"
#ChiRasterName = Base_file+"_SO.bil"

#BR = BaseRaster(BackgroundRasterName, Directory)
#BR.set_raster_type("Terrain")
#print(BR._colourmap)
#BR.show_raster()

#BR.set_colourmap("RdYlGn")
#BR.show_raster()

#label_size = 100
#rcParams['font.family'] = 'sans-serif'
#rcParams['font.sans-serif'] = ['arial']
#rcParams['font.size'] = label_size
#rcParams['lines.linewidth']  = 1.5

plt.clf()
MF = MapFigure(BackgroundRasterName, Directory,coord_type="UTM_km",colourbar_location='None')
#MF.add_drape_image(DrapeRasterName,Directory,alpha = 0.4)
#MF.add_drape_image(ChiRasterName,Directory,colourmap = "cubehelix",alpha = 0.4, show_colourbar = True)
#MF.show_plot()
MF.save_fig(fig_width_inches = 12, dpi=300, FileName=Base_file+'.png', format='png')

# Customise the DrapePlot
#dp.make_drape_colourbar(cbar_label=colourbar_label)
#dp.set_fig_axis_labels()

#dp.show_plot()
