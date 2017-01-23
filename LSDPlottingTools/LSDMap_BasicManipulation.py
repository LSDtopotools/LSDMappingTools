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
# This takes a grouping list of basin keys and transforms it into a list
# of junction names
#==============================================================================
def BasinKeyToJunction(grouped_data_list,basin_info_csv):
    
    thisPointData = LSDMPD.LSDMap_PointData(basin_info_csv)
    
    thisJunctionData = thisPointData.QueryData("outlet_junction")
    
    junction_grouped_list = []

    if not grouped_data_list:
        return grouped_data_list
    else:
        for group in grouped_data_list:
            this_list = []
            for element in group:
                this_list.append(thisJunctionData[element])
            junction_grouped_list.append(this_list)

    print junction_grouped_list 
    return junction_grouped_list
    
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## This function orders basins in a sequence, so that plots can be made
## as a function of distance, for example
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        
def BasinOrderer(Basin_csv_name, FileName, criteria_string,reverse=False, 
                 exclude_criteria_string = "None", exlude_criteria_greater = False, 
                 exclude_criteria_value = 0 ):        

    #========================================================
    # now we need to label the basins
    # Now we get the chi points
    EPSG_string = LSDMap_IO.GetUTMEPSG(FileName)
    print "EPSG string is: " + EPSG_string
    
    thisPointData = LSDMPD.LSDMap_PointData(Basin_csv_name) 
    
    # convert to easting and northing
    [easting,northing] = thisPointData.GetUTMEastingNorthingFromQuery(EPSG_string,"outlet_latitude","outlet_longitude")
 
    these_data = thisPointData.QueryData("outlet_junction")
    #print M_chi
    these_data = [int(x) for x in these_data]

    wanted_data = thisPointData.QueryData(criteria_string)
    wd = np.asarray(wanted_data)

    # Now we exclude some of the basins
    #if exclude_criteria_string != "None":
    #    exclude_data = thisPointData.QueryData(exclude_criteria_string)
        
    indices = np.argsort(wd)
    
    print("data is: ")
    print wd

    if reverse:
        indices = indices[::-1]
           
    sorted_basins = np.array(these_data)[indices]


    print sorted_basins       


#==============================================================================
# This function takes groups of data and then resets values in a
# raster to mimic these values
#============================================================================== 
def RedefineIntRaster(rasterArray,grouped_data_list,spread):
    
    counter = 0
    if not grouped_data_list:
        return rasterArray
    else:
        for group in grouped_data_list:
            for element in group:
                rasterArray[rasterArray == element] = counter
                counter= counter+1
                            
            counter = counter+spread
    return rasterArray
            
            
    
    

    # You need to manipulate the groupings 



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
# This does a very basic swath analysis in one direction
# if axis is 0, this is along x axis, if axis is 1, is along y axis
# otherwise will throw error
#==============================================================================         
def SimpleSwath(path, file1, axis):
   
    # make sure names are in correct format
    NewPath = LSDOst.AppendSepToDirectoryPath(path)
    
    raster_file1 = NewPath+file1

    # get some information about the raster 
    NDV, xsize, ysize, GeoT, Projection, DataType = LSDMap_IO.GetGeoInfo(raster_file1)
    
    print "NDV is: "
    print NDV
    
    if NDV == None:
        NDV = -9999
        print "No NDV defined"
    
    Raster1 = LSDMap_IO.ReadRasterArrayBlocks(raster_file1,raster_band=1)
    
    
    #nan_raster = Raster1[Raster1==NDV]=np.nan 
    #print nan_raster
    
    #now mask the nodata
    masked_Raster1  = np.ma.masked_values(Raster1, NDV)
    
    means = np.mean(masked_Raster1, axis)
    medians = np.median(masked_Raster1, axis)
    std_deviations = np.std(masked_Raster1, axis)
    twentyfifth_percentile = np.percentile(masked_Raster1, 25, axis)
    seventyfifth_percentile = np.percentile(masked_Raster1, 75, axis)
 
    # This stuff only works with numpy 1.8 or later, wich we don't have
    #means = np.nanmean(nan_raster, axis)
    #medians = np.nanmedian(nan_raster, axis)
    #std_deviations = np.nanstd(nan_raster, axis)
    #twentyfifth_percentile = np.nanpercentile(nan_raster, 25, axis)
    #seventyfifth_percentile = np.nanpercentile(nan_raster, 75, axis)
   
    #print means
    #print medians
    #print std_deviations
    #print twentyfifth_percentile
    #print seventyfifth_percentile    
    
    
    return means,medians,std_deviations,twentyfifth_percentile,seventyfifth_percentile   

        
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