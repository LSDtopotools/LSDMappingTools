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


import matplotlib.pyplot as plt
import matplotlib

# Force matplotlib to not use any Xwindows backend.
#matplotlib.use('Agg')

from matplotlib import rcParams
import matplotlib.cm as cm
from LSDMapFigure.PlottingRaster import MapFigure
from LSDMapFigure.PlottingRaster import BaseRaster

from LSDPlottingTools import colours as lsdcolours
from LSDPlottingTools import init_plotting_DV
import sys

#sys.path.append("PATH/TO/LSDPlottingTools/")
#
#init_plotting_DV()

#label_size = 100
#rcParams['font.family'] = 'sans-serif'
#rcParams['font.sans-serif'] = ['arial']
#rcParams['font.size'] = label_size
#rcParams['lines.linewidth']  = 1.5   

#Directory = "C:\\VagrantBoxes\\LSDTopoTools\\Topographic_projects\\Meghalaya\\divides\\"
#Directory = "T:\\analysis_for_papers\\Meghalaya\\divide_migration\\"
Directory = "/home/smudd/SMMDataStore/analysis_for_papers/Meghalaya/divide_migration/"
Base_file = "Mega_divide"


BackgroundRasterName = Base_file+".bil"
DrapeRasterName = Base_file+"_hs.bil"
ChiRasterName = Base_file+"_chi_coord.bil"

#BR = BaseRaster(BackgroundRasterName, Directory)
#BR.set_raster_type("Terrain")
#print(BR._colourmap)
#BR.show_raster()

#BR.set_colourmap("RdYlGn")
#BR.show_raster()

   


plt.clf() 
MF = MapFigure(BackgroundRasterName, Directory,coord_type="UTM_km")
MF.add_drape_image(DrapeRasterName,Directory,alpha = 0.4)
MF.add_drape_image(ChiRasterName,Directory,colourmap = "cubehelix",alpha = 0.4, show_colourbar = True)
#MF.show_plot()
ImageName = Directory+"TestNewArtist.png" 
fig_size_inches = 6
ax_style = "Madhouse"
MF.save_fig(fig_width_inches = fig_size_inches, FigFileName = ImageName, axis_style = ax_style)



# Customise the DrapePlot
#dp.make_drape_colourbar(cbar_label=colourbar_label)
#dp.set_fig_axis_labels()

#dp.show_plot()