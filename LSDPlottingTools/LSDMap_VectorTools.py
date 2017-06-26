## LSDMap_BasicManipulation.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## These functions are tools to deal with vector data using shapely
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## FJC
## 26/06/17
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

import numpy as np
from . import LSDMap_GDALIO as LSDMap_IO

#--------------------------------------------------------------------------------#
# BASIN FUNCTIONS
# These functions do various operations on basin polygons
#
#--------------------------------------------------------------------------------#
def GetBasinOutlines(DataDirectory, fname_prefix):
    """
    This function takes in the raster of basins and gets a dict of basin polygons,
    where the key is the basin key and the value is a shapely polygon of the basin

    Args:
        DataDirectory (str): the data directory with the basin raster
        fname_prefix (str): the prefix for the DEM

    Returns:
        list of shapely polygons with the basins

    Author: FJC
    """

    # read in the basins raster
    basin_name = fname_prefix+"_AllBasins.bil"
    OutputShapefile = fname_prefix+'_basins.shp'

    # polygonise the raster
    BasinDict = LSDMap_IO.PolygoniseRaster(DataDirectory, basin_name, OutputShapefile)
    return BasinDict

def GetBasinCentroids(DataDirectory, fname_prefix):
    """
    This function takes in the raster of basins and returns a dict where the
    key is the basin key and the value is the shapely point of the centroid

    Args:
        DataDirectory (str): the data directory with the basin raster
        fname_prefix (str): the prefix for the DEM

    Returns:
        dict of centroid points

    Author: FJC
    """
    from shapely.geometry import Point

    # get the basin polygons
    BasinDict = GetBasinOutlines(DataDirectory, fname_prefix)

    # get the centroids
    CentroidDict = {}
    for basin_key, basin in BasinDict.iteritems():
        CentroidDict[basin_key] = Point(basin.centroid)

    return CentroidDict

def GetPointWithinBasins(DataDirectory,fname_prefix):
    """
    This function takes in the raster of basin and returns a dict where the
    key is the basin key and the value is a shapely point that is representative
    of the basin (guaranteed to be within the polygon)

    Args:
        DataDirectory (str): the data directory with the basin raster
        fname_prefix (str): the prefix for the DEM

    Returns:
        dict of representative points

    Author: FJC
    """
    from shapely.geometry import Point

    # get the basin polygons
    BasinDict = GetBasinOutlines(DataDirectory, fname_prefix)

    # get the centroids
    PointDict = {}
    for basin_key, basin in BasinDict.iteritems():
        PointDict[basin_key] = Point(basin.representative_point())

    return PointDict
