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

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
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



def BasemapExtentSizer(FigWidthInches, FigHeightInches):
    """
    The returns the aspect ratio that you need for the main figure given the dimensions of the figure.
    Basically this just calculates the padding needed for the axis labels
    
    Args:
        FigWidthInches (float): How wide you want the basemap
        FigHeightInches (float): How high you want your basemap
        
    Returns: The aspect ratio (a float) of the basemap extent.
    Author: SMM
    
    Date 03/02/2018
    """
    
    x_width = FigWidthInches - 0.3  # The 0.3 inches is the text)
    y_width = FigHeightInches - 0.2
    
    aspect_ratio = x_width/y_width
    return aspect_ratio

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
    
    
def GenerateBasemapImage(DataDirectory, RasterFile, FigWidthInches = 4, FigHeightInches = 3, bm_width = 2000000, bm_height = 2000000, projection = 'lcc',resolution = 'l', lat_0 = 0, lon_0 = 0, lat_1 = 45,lat_2 = 55, satellite_height = 10000000, FigFormat = "png", fig_dpi = 500, out_fname_prefix = ""):
    """
    This makes the basemap image. 
    
    Args:
        DataDirectory (str): The directory of the raster file
        RasterFile (str): the name of the raster file (include extension)
        FigWidthInches (float): How wide you want the basemap
        FigHeightInches (float): How high you want your basemap
        bm_width (float): The width in metres of your basemap
        bm_height (float): The height in metres covered by your basemap
        projection (str): The projection of your basemap. See basemap docs for details
        resolution (str): Resolution. See basmap documentation. Default is "l" (for low) since higher resolution ("high" or "full") must be installed separately with basemap (and is very large)
        lat_0 (flt): Latitude of centre of your map
        lon_0 (flt): Longitude of centre of your map
        lat_1 (flt): See basemap documentation 
        lon_1 (flt): See basemap documentation
        satellite_height (flt): The satellite height in metres for geostationary projections
        FigFormat (str): Figure format, can be `png`, `svg`, `pdf`, etc.
        fig_dpi (int): Dots per inch of your figure 
        out_fname_prefix (str): The prefix of the image file. If blank uses the fname_prefix
    
    
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

    # Get the name of the image
    if len(out_fname_prefix) == 0:
        FigFileName = DataDirectory+Base_file+"_basemap."+fig_format
    else:
        FigFileName = DataDirectory+out_fname_prefix+"_basemap."+fig_format    

    
    # Now for the basemap 
    # setup Lambert Conformal basemap.
    m = Basemap(width=bm_width,height=bm_width,projection=projection,
                resolution=resolution,lat_1=lat_1 ,lat_2=lat_2,lat_0=lat_0,lon_0=lon_0, satellite_height = satellite_height, area_thresh = 100000)    
 
    # create the shapefile
    LSDMGDAL.CreateShapefileOfRasterFootprint(DataDirectory, RasterFile)
    m.readshapefile(Shape_name, "footprint")

    # draw coastlines.
    m.drawcoastlines(linewidth = 0.5)
    # draw a boundary around the map, fill the background.
    # this background will end up being the ocean color, since
    # the continents will be drawn on top.
    m.drawmapboundary(fill_color='snow')
    # fill continents, set lake color same as ocean color.
    m.fillcontinents(color='lightgray',lake_color='snow')
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
    pc.set_facecolor("crimson")
    ax.add_collection(pc)  

    plt.savefig(FigFileName,format=FigFormat,dpi=fig_dpi)    
    
    
def GenerateBasemapImageAutomated(DataDirectory, RasterFile, FigWidthInches = 4, FigHeightInches = 3, FigFormat = "png", fig_dpi = 500, regional_extent_multiplier = 5, label_spacing_multiplier = 0.5, out_fname_prefix = ""):
    """
    This makes the basemap image. Uses data from the raster to size the figure and locate the centrepoint 
    
    Args:
        DataDirectory (str): The directory of the raster file
        RasterFile (str): the name of the raster file (include extension)
        FigWidthInches (flt): How wide you want the basemap
        FigHeightInches (float): How high you want your basemap
        FigFormat (str): Figure format, can be `png`, `svg`, `pdf`, etc.
        fig_dpi (int): Dots per inch of your figure 
        regional_extent_multiplier (float): How much bigger you want the extent vs the size of the raster 
        label_spacing_multiplier (float): If the meridians and parallels are too close, increase this number. Default of 0.5
        out_fname_prefix (str): The prefix of the image file. If blank uses the fname_prefix
    
    
    Author: SMM
    
    Date: 01/02/2018
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
    
    # Get the name of the image
    if len(out_fname_prefix) == 0:
        FigFileName = DataDirectory+Base_file+"_basemap."+FigFormat
    else:
        FigFileName = DataDirectory+out_fname_prefix+"_basemap."+FigFormat   
    
    # Now we get the extents from the raster
    centre_lat, centre_long, extent_lat, extent_long, xproj_extent, yproj_extent = LSDMGDAL.GetCentreAndExtentOfRaster(DataDirectory, RasterFile)
    

    
    
    # Calculate the aspect ratio
    aspect_ratio = BasemapExtentSizer(FigWidthInches,  FigHeightInches)
    
    # Figure out the longest dimension
    long_dimension = xproj_extent
    if yproj_extent > long_dimension:
        long_dimension = yproj_extent
    
    # Get the full extent by mulitplying the longest extent by the multiplier
    full_dimension = long_dimension*regional_extent_multiplier
         
        
    # now get the two dimensions for the extent of the figure
    if aspect_ratio > 1:
        # This is when the figure is wider than tall
        x_ext = full_dimension*aspect_ratio
        y_ext = full_dimension
    else:
        x_ext = full_dimension
        y_ext = full_dimension*aspect_ratio
       
    # Now for the basemap 
    # setup Lambert Conformal basemap.
    m = Basemap(width=x_ext,height=y_ext,projection="lcc",
                resolution="l",lat_1=45 ,lat_2=55,lat_0=centre_lat,lon_0=centre_long, satellite_height = 100000, area_thresh = 100000)    
 
    # create the shapefile
    LSDMGDAL.CreateShapefileOfRasterFootprint(DataDirectory, RasterFile)
    m.readshapefile(Shape_name, "footprint")

    # draw coastlines.
    m.drawcoastlines(linewidth = 0.5)
    # draw a boundary around the map, fill the background.
    # this background will end up being the ocean color, since
    # the continents will be drawn on top.
    m.drawmapboundary(fill_color='snow')
    # fill continents, set lake color same as ocean color.
    m.fillcontinents(color='lightgray',lake_color='snow')
    m.drawcountries()
    
    #==========================================
    # draw parallels and meridians.    
    # Calculate the spacing of the meridians and parallels
    max_latlong_ext = extent_lat
    if extent_long < max_latlong_ext:
        max_latlong_ext = extent_long
    max_latlong_ext = max_latlong_ext*regional_extent_multiplier
    print("The maximum extent is:"+str(max_latlong_ext))
    
    latlong_label_spacing = int(label_spacing_multiplier*max_latlong_ext+0.5)
    print("And the label spacing is: "+str(latlong_label_spacing))
    start_lat = int(centre_lat - max_latlong_ext*2 - 0.5)
    end_lat = int(centre_lat + max_latlong_ext*2+0.5)
    start_long = int(centre_long - max_latlong_ext*2 -0.5)
    end_long = int(centre_long + max_latlong_ext*2+0.5)
          
    # label parallels on right and top
    # meridians on bottom and left   
    parallels = np.arange(start_lat,end_lat,latlong_label_spacing)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels,labels=[False,True,True,False])
    meridians = np.arange(start_long,end_long,latlong_label_spacing)
    m.drawmeridians(meridians,labels=[True,False,False,True])
    #==========================================
    

    # Make a patch from the shapefile
    # All this stuff from:
    # http://www.datadependence.com/2016/06/creating-map-visualisations-in-python/
    df_poly = pd.DataFrame({
            'shapes': [Polygon(np.array(shape), True) for shape in m.footprint]})    
    
    #df_poly = df_poly.merge(new_areas, on='area', how='left')
    #cmap = plt.get_cmap('Oranges')   
    pc = PatchCollection(df_poly.shapes, zorder=2, alpha = 0.5)
    pc.set_facecolor("crimson")
    ax.add_collection(pc)  

    plt.savefig(FigFileName,format=FigFormat,dpi=fig_dpi)    
    
