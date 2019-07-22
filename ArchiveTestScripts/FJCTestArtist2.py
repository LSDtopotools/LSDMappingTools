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
from LSDPlottingTools import LSDMap_VectorTools as LSDMap_VT

# Generate random colormap
def rand_cmap(nlabels, type='bright', first_color_black=True, last_color_black=False, verbose=True):
    """
    Creates a random colormap to be used together with matplotlib. Useful for segmentation tasks
    :param nlabels: Number of labels (size of colormap)
    :param type: 'bright' for strong colors, 'soft' for pastel colors
    :param first_color_black: Option to use first color as black, True or False
    :param last_color_black: Option to use last color as black, True or False
    :param verbose: Prints the number of labels and shows the colormap. True or False
    :return: colormap for matplotlib
    """
    from matplotlib.colors import LinearSegmentedColormap
    import colorsys
    import numpy as np

    if type not in ('bright', 'soft'):
        print ('Please choose "bright" or "soft" for type')
        return

    if verbose:
        print('Number of labels: ' + str(nlabels))

    # Generate color map for bright colors, based on hsv
    if type == 'bright':
        randHSVcolors = [(np.random.uniform(low=0.0, high=1),
                          np.random.uniform(low=0.2, high=1),
                          np.random.uniform(low=0.9, high=1)) for i in xrange(nlabels)]

        # Convert HSV list to RGB
        randRGBcolors = []
        for HSVcolor in randHSVcolors:
            randRGBcolors.append(colorsys.hsv_to_rgb(HSVcolor[0], HSVcolor[1], HSVcolor[2]))

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]

        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Generate soft pastel colors, by limiting the RGB spectrum
    if type == 'soft':
        low = 0.6
        high = 0.95
        randRGBcolors = [(np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high)) for i in xrange(nlabels)]

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]
        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Display colorbar
    if verbose:
        from matplotlib import colors, colorbar
        from matplotlib import pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(15, 0.5))

        bounds = np.linspace(0, nlabels, nlabels + 1)
        norm = colors.BoundaryNorm(bounds, nlabels)

        cb = colorbar.ColorbarBase(ax, cmap=random_colormap, norm=norm, spacing='proportional', ticks=None,
                                   boundaries=bounds, format='%1i', orientation=u'horizontal')

    return random_colormap

#sys.path.append("PATH/TO/LSDPlottingTools/")
#
#init_plotting_DV()

Directory = "/media/fionaclubb/terrace_lidar/Mendocino_TJ/mendocino/coastal_california/"
#Directory = "T:\\analysis_for_papers\\Meghalaya\\divide_migration\\"
Base_file = "California_clip"

#BackgroundRasterName = Base_file+".bil"
DrapeRasterName = Base_file+"_hs.bil"
BasinsName = Base_file+"_basins.bil"

#BR = BaseRaster(BackgroundRasterName, Directory)
#BR.set_raster_type("Terrain")
#print(BR._colourmap)
#BR.show_raster()

#BR.set_colourmap("RdYlGn")
#BR.show_raster()

#label_size = 100
#rcParams['font.family'] = 'sans-serif'
#rcParams['font.sans-serif'] = ['Liberation Sans']
#rcParams['font.size'] = label_size
#rcParams['lines.linewidth']  = 1.5

plt.clf()
MF = MapFigure(DrapeRasterName, Directory,coord_type="UTM_km",colourbar_location='None')
# add the basins drape
cmap = cm.Set1
MF.add_drape_image(BasinsName, Directory, colourmap = cmap, alpha = 0.5, colorbarlabel='Basin ID', show_colourbar = False)
# add the basin outlines
Basins = LSDMap_VT.GetBasinOutlines(Directory, BasinsName)
MF.plot_polygon_outlines(Basins, linewidth=0.8)
#MF.add_drape_image(ChiRasterName,Directory,colourmap = "cubehelix",alpha = 0.4, show_colourbar = True)
#MF.show_plot()
MF.save_fig(fig_width_inches = 12, FigFileName=Directory+Base_file+'.png', FigFormat='png', Fig_dpi=300, transparent=True)

# Customise the DrapePlot
#dp.make_drape_colourbar(cbar_label=colourbar_label)
#dp.set_fig_axis_labels()

#dp.show_plot()
