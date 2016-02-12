## LSDMap_BasicManipulation.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## These functions are tools to deal with rasters
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## SMM
## 26/07/2014
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

import osgeo.gdal as gdal
import numpy as np
from osgeo import osr
from os.path import exists
from osgeo.gdalconst import GA_ReadOnly
from numpy import uint8
import LSDMap_GDALIO as LSDMap_IO

#==============================================================================
# THis function takes a raster an writes a new raster where everything below a
# threshold is set to nodata
#==============================================================================
def SetNoDataBelowThreshold(raster_filename,new_raster_filename, threshold = 0, driver_name = "ENVI", NoDataValue = -9999):
    
    # read the data
    rasterArray = LSDMap_IO.ReadRasterArrayBlocks(raster_filename)
    print "Read the data"
    
    # set any point on the raster below the threshold as nodata
    rasterArray[rasterArray <= threshold] = NoDataValue
    print "Reset raster values"
    
    # write the data to a new file
    LSDMap_IO.array2raster(raster_filename,new_raster_filename,rasterArray,driver_name, NoDataValue)
    print "Wrote raster"
#==============================================================================

#==============================================================================
# This function sets all nodata values to a constant value
#==============================================================================
def SetToConstantValue(raster_filename,new_raster_filename, constant_value, driver_name = "ENVI"):

    # get the nodata value
    NoDataValue =  LSDMap_IO.getNoDataValue(raster_filename)
               
    # read the data
    rasterArray = LSDMap_IO.ReadRasterArrayBlocks(raster_filename)
    print "Read the data"   
    
    # set any nodata to a constant value
    rasterArray[rasterArray != NoDataValue] = constant_value
    print "Changed to a constant value"
    
    # write the data to a new file
    LSDMap_IO.array2raster(raster_filename,new_raster_filename,rasterArray,driver_name, NoDataValue)
    print "Wrote raster"
    
    
    
    
    