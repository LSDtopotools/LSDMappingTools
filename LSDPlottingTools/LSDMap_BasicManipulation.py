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
from glob import glob
import LSDOSystemTools as LSDOst
import LSDMap_GDALIO as LSDMap_IO
import LSDMap_BasicPlotting as LSDMBP
import LSDMap_PointData as LSDMPD


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

#==============================================================================
# This function calcualtes a hillshade and writes to file
#==============================================================================    
def GetHillshade(raster_filename,new_raster_filename, azimuth = 315, angle_altitude = 45, driver_name = "ENVI", NoDataValue = -9999):

    # get the hillshade
    hillshade_raster = LSDMBP.Hillshade(raster_filename, azimuth, angle_altitude)

    # write to file
    LSDMap_IO.array2raster(raster_filename,new_raster_filename,hillshade_raster,driver_name, NoDataValue)         
    
#==============================================================================
# This function takes all the csv files in a directory and converts to 
# GeoJSON files
#==============================================================================     
def ConvertAllCSVToGeoJSON(path):
    
    # make sure names are in correct format
    NewPath = LSDOst.AppendSepToDirectoryPath(path)
    
    print "The formatted path is: " + NewPath
    
        
    for FileName in glob(NewPath+"*.csv"): 
        print "filename is: " + FileName
        
        thisPointData = LSDMPD.LSDMap_PointData(FileName)
        thisPointData.TranslateToReducedGeoJSON(FileName)
        
        
#==============================================================================
# This function takes all the csv files in a directory and converts to 
# GeoJSON files
#==============================================================================     
def ConvertAllCSVToShapefile(path):
    
    # make sure names are in correct format
    NewPath = LSDOst.AppendSepToDirectoryPath(path)
    
    print "The formatted path is: " + NewPath
    
        
    for FileName in glob(NewPath+"*.csv"): 
        print "filename is: " + FileName
        
        thisPointData = LSDMPD.LSDMap_PointData(FileName)
        thisPointData.TranslateToReducedShapefile(FileName)   

#==============================================================================
# This does a basic mass balance. 
# Assumes all units are metres
#==============================================================================         
def RasterMeanValue(path, file1):
    
    # make sure names are in correct format
    NewPath = LSDOst.AppendSepToDirectoryPath(path)
    
    raster_file1 = NewPath+file1
    
    NPixels = LSDMap_IO.GetNPixelsInRaster(raster_file1)

    Raster1 = LSDMap_IO.ReadRasterArrayBlocks(raster_file1,raster_band=1)
    
    mean_value = np.sum(Raster1)/float(NPixels)     
  
    return mean_value   



        
#==============================================================================
# This does a basic mass balance. 
# Assumes all units are metres
#==============================================================================         
def BasicMassBalance(path, file1, file2):
    
    # make sure names are in correct format
    NewPath = LSDOst.AppendSepToDirectoryPath(path)
    
    raster_file1 = NewPath+file1
    raster_file2 = NewPath+file2
    
    PixelArea = LSDMap_IO.GetPixelArea(raster_file1)
    print "PixelArea is: " + str(PixelArea) 
    
    print "The formatted path is: " + NewPath
    Raster1 = LSDMap_IO.ReadRasterArrayBlocks(raster_file1,raster_band=1)
    Raster2 = LSDMap_IO.ReadRasterArrayBlocks(raster_file2,raster_band=1)
    
    NewRaster = np.subtract(Raster2,Raster1)
    
    mass_balance = np.sum(NewRaster)*PixelArea

    print "linear dif " + str(np.sum(NewRaster))    
        
    return mass_balance       