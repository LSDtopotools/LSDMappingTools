## LSDMap_BasemapTools.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## These functions are tools to create a basemap
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## SMM
## 24/01/2018
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#from __future__ import absolute_import, division, print_function, unicode_literals
from __future__ import absolute_import, division, print_function


import matplotlib
matplotlib.use('Agg')

import LSDPlottingTools.LSDMap_GDALIO as LSDMGDAL
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.cm 
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize
import os



def GenerateExtentShapefile(DataDirectory, RasterFile):
    """
    This just wraps a LSDMap_GDALIO script
    
    Args:
        DataDirectory (str): the data directory with the basin raster
        RasterFile (str): the name of the raster

    Returns:
        Shapefile of the raster footprint. Has "_footprint" in filename.

    Author: SMM
    
    Date: 24/01/2018
    
    """
    LSDMGDAL.CreateShapefileOfRasterFootprint(DataDirectory, RasterFile)
    
    
def GenerateBasemapImage(DataDirectory, RasterFile, FigWidthInches = 4, FigHeightInches = 3, bm_width = 2000000, bm_height = 2000000, projection = 'lcc',resolution = 'h', lat_0 = 0, lon_0 = 0, lat_1 = 45,lat_2 = 55, satellite_height = 10000000, FigFormat = "png", fig_dpi = 500):
    """
    This makes the basemap image. 
    
    
    Author: SMM
    
    Date: 24/01/2018
    """
    
    # Make sure data directory is in correct format
    if not DataDirectory.endswith(os.sep):
        print("You forgot the separator at the end of the directory, appending...")
        DataDirectory = DataDirectory+os.sep
    
    # Set up the figure. This is needed to both size the figure and get the axis handle for plotting polygons
    fig, ax = plt.subplots(figsize=(FigWidthInches, FigHeightInches))
    
    # get some filenames
    RasterSplit = RasterFile.split(".")
    Raster_prefix = RasterSplit[0]
    Shape_name = DataDirectory+Raster_prefix+"_footprint"
    SName = "Shape"
    
    FigFileName = DataDirectory+Raster_prefix+"_basemap."+FigFormat
    
    
    
    # Now for the basemap 
    # setup Lambert Conformal basemap.
    m = Basemap(width=bm_width,height=bm_width,projection=projection,
                resolution=resolution,lat_1=lat_1 ,lat_2=lat_2,lat_0=lat_0,lon_0=lon_0, satellite_height = satellite_height)    
 
    # create the shapefile
    LSDMGDAL.CreateShapefileOfRasterFootprint(DataDirectory, RasterFile)
    m.readshapefile(Shape_name, "footprint")

    # draw coastlines.
    m.drawcoastlines()
    # draw a boundary around the map, fill the background.
    # this background will end up being the ocean color, since
    # the continents will be drawn on top.
    m.drawmapboundary(fill_color='whitesmoke')
    # fill continents, set lake color same as ocean color.
    m.fillcontinents(color='silver',lake_color='whitesmoke')
    # draw parallels and meridians.
    # label parallels on right and top
    # meridians on bottom and left
    parallels = np.arange(0.,90,5.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels,labels=[False,True,True,False])
    meridians = np.arange(10.,351.,5.)
    m.drawmeridians(meridians,labels=[True,False,False,True])
    m.drawcountries()

    # Make a patch from the shapefile
    # All this stuff from:
    # http://www.datadependence.com/2016/06/creating-map-visualisations-in-python/
    df_poly = pd.DataFrame({
            'shapes': [Polygon(np.array(shape), True) for shape in m.footprint]})    
    
    #df_poly = df_poly.merge(new_areas, on='area', how='left')
    #cmap = plt.get_cmap('Oranges')   
    pc = PatchCollection(df_poly.shapes, zorder=2, alpha = 0.5)
    pc.set_facecolor("dimgrey")
    ax.add_collection(pc)  

    plt.savefig(FigFileName,format=FigFormat,dpi=fig_dpi)    
    
    

