# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 12:55:23 2017

@author: smudd
"""


import matplotlib.cm as cm
from LSDMapArtist.drapeplot_experimental import DrapeAxes

from LSDPlottingTools import colours as lsdcolours
from LSDPlottingTools import init_plotting_DV
import sys

#sys.path.append("PATH/TO/LSDPlottingTools/")
#
#init_plotting_DV()

Directory = "T:\\analysis_for_papers\\Meghalaya\\divide_migration\\"
Base_file = "Mega_divide"
BackgroundRasterName = Base_file+".bil"
DrapeRasterName = Base_file+"_hs.bil"
#DrapeRasterName = "BoscastleErodeDiff_GRID_UNI_TLIMM.bil"





# Standard colourmap
Colourmap = "RdYlGn"

#Non-linear colourmap
##ColourLevels = lsdcolours.nonlinear_colourmap.create_levels(-3.0, 3.0, -0.2, 0.2, -0.5, 0.5)
##Colourmap = lsdcolours.nonlinear_colourmap("seismic", ColourLevels)

# Transformed colourmap
#c = lsdcolours.TransformedColourmap(lambda x: x/2+0.5, cm.jet)

drape_min_threshold = None
drape_max_threshold = None
colourbar_label = "Yo (m)"

#raster = BaseRaster(RasterName, DataDirectory)
dp = DrapeAxes(DrapeRasterName, BackgroundRasterName, Directory,
                      Colourmap, background_type="Hillshade", 
                      show_background_colourbar=False,
                      show_drape = True,
                      colourbar_label=colourbar_label)



# Customise the DrapePlot
#dp.make_drape_colourbar(cbar_label=colourbar_label)
#dp.set_fig_axis_labels()

dp.show_plot()