## LSDMap_BasicManipulation.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## These functions are tools to deal with rasters
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## SMM
## 26/07/2014
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from . import LSDMap_OSystemTools as LSDOst
from . import LSDMap_GDALIO as LSDMap_IO
from pyproj import Proj, transform


def GetUTMEastingNorthing(EPSG_string,latitude,longitude):
    """This returns the easting and northing for a given latitude and longitide

    Args:
        ESPG_string (str): The ESPG code. 326XX is for UTM north and 327XX is for UTM south
        latitude (float): The latitude in WGS84
        longitude (float): The longitude in WGS84

    Returns:
        easting,northing The easting and northing in the UTM zone of your selection

    Author:
        Simon M Mudd
    """

    #print "Yo, getting this stuff: "+EPSG_string
    # The lat long are in epsg 4326 which is WGS84
    inProj = Proj(init='epsg:4326')
    outProj = Proj(init=EPSG_string)
    ea,no = transform(inProj,outProj,longitude,latitude)

    return ea,no


def ConvertNorthingForImshow(RasterName,Northing):
    """This returns a northing that is inverted using the minimum and maximum values from the raster for use in imshow (because imshow inverts the raster)

    Args:
        RasterName (str): The raster's name with full path and extension
        Northing (float): The northing coordinate in metres (from UTM WGS84)

    Returns:
        ConvertedNorthing The northing inverted from top to bottom

    Author: SMM
    """

    extent_raster = LSDMap_IO.GetRasterExtent(RasterName)


    Northing_converted = extent_raster[3]-Northing
    Northing_converted = Northing_converted+extent_raster[2]

    return Northing_converted

#==============================================================================
# THis function takes a raster an writes a new raster where everything below a
# threshold is set to nodata
#==============================================================================
def SetNoDataBelowThreshold(raster_filename,new_raster_filename, threshold = 0, driver_name = "ENVI", NoDataValue = -9999):
    """This takes a raster an then converts all data below a threshold to nodata, it then prints the resulting raster.

    Args:
        raster_filename (str): The raster's name with full path and extension
        new_raster_filename (str): The name of the raster to be printed
        threshold (float): Data below this in the original raster will be converted to nodata.
        driver_name (str): The raster format (see gdal documentation for options. LSDTopoTools used "ENVI" format.)
        NoDataValue (float): The nodata value. Usually set to -9999.

    Returns:
        None, but prints a new raster to file

    Author: SMM
    """



    # read the data
    rasterArray = LSDMap_IO.ReadRasterArrayBlocks(raster_filename)
    print("Read the data")

    # set any point on the raster below the threshold as nodata
    rasterArray[rasterArray <= threshold] = NoDataValue
    print("Reset raster values")

    # write the data to a new file
    LSDMap_IO.array2raster(raster_filename,new_raster_filename,rasterArray,driver_name, NoDataValue)
    print("Wrote raster")
#==============================================================================

#==============================================================================
# This function sets all nodata values to a constant value
#==============================================================================
def SetToConstantValue(raster_filename,new_raster_filename, constant_value, driver_name = "ENVI"):
    """This takes a raster an then converts all non-nodata to a constant value.

    This is useful if you want to make masks, for example to have blocks of single erosion rates for cosmogenic calculations.

    Args:
        raster_filename (str): The raster's name with full path and extension
        new_raster_filename (str): The name of the raster to be printed
        constant_value (float): All non-nodata will be converted to this value in a new raster.
        driver_name (str): The raster format (see gdal documentation for options. LSDTopoTools used "ENVI" format.)
        NoDataValue (float): The nodata value. Usually set to -9999.

    Returns:
        None, but prints a new raster to file

    Author: SMM
    """

    # get the nodata value
    NoDataValue =  LSDMap_IO.getNoDataValue(raster_filename)

    # read the data
    rasterArray = LSDMap_IO.ReadRasterArrayBlocks(raster_filename)
    print("Read the data")

    # set any nodata to a constant value
    rasterArray[rasterArray != NoDataValue] = constant_value
    print("Changed to a constant value")

    # write the data to a new file
    LSDMap_IO.array2raster(raster_filename,new_raster_filename,rasterArray,driver_name, NoDataValue)
    print("Wrote raster")

#==============================================================================
# This function calcualtes a hillshade and writes to file
#==============================================================================
def GetHillshade(raster_filename,new_raster_filename, azimuth = 315, angle_altitude = 45, driver_name = "ENVI", NoDataValue = -9999):
    """This calls the hillshade function from the basic manipulation package, but then prints the resulting raster to file.

   Args:
        raster_filename (str): The raster's name with full path and extension
        new_raster_filename (str): The name of the raster to be printed
        azimuth (float): Azimuth angle (compass direction) of the sun (in degrees).
        angle_altitude (float):Altitude angle of the sun.
        driver_name (str): The raster format (see gdal documentation for options. LSDTopoTools used "ENVI" format.)
        NoDataValue (float): The nodata value. Usually set to -9999.

    Returns:
        None, but prints a new raster to file.

    Author: SMM
    """
    # avoid circular import
    from . import LSDMap_BasicPlotting as LSDMBP
    # get the hillshade
    hillshade_raster = LSDMBP.Hillshade(raster_filename, azimuth, angle_altitude)

    # write to file
    LSDMap_IO.array2raster(raster_filename,new_raster_filename,hillshade_raster,driver_name, NoDataValue)



#==============================================================================
# This takes a grouping list of basin keys and transforms it into a list
# of junction names
#==============================================================================
def BasinKeyToJunction(grouped_data_list,thisPointData):
    """This takes a basin_info_csv file (produced by several LSDTopoTools routies) and spits out lists of the junction numbers (it converts basin numbers to junction numbers).


    Args:
        grouped_data_list (int list): A list of list of basin numbers
        thisPointData (str): A point data object with the basins

    Returns:
        Junction_grouped_list: the junction numbers of the basins.

    Author: SMM
    """

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

    print(junction_grouped_list)
    return junction_grouped_list


def BasinOrderToBasinRenameList(basin_order_list):
    """When we take data from the basins they will be numbered accoring to their junction rank, which is controlled by flow routing.

    The result is often numbered basins that have something that appears random to human eyes.
    We have developed a routine to renumber these basins.
    However, the way this works is to find a basin number and rename in the profile plots,
    such that when it finds a basin number it will rename that.
    So if you want to rename the seventh basin 0, you need to give a list where the seventh element is 0.

    This is a pain because how one would normally order basins would be to look at the image of the basin numbers,
    and then write the order in which you want those basins to appear.

    This function converts between these two lists. You give the function the order you want the basins to appear,
    and it gives a renaming list.

    Args:
        basin_order_list (int): the list of basins in which you want them to appear in the numbering scheme

    Return:
        The index into the returned basins

    Author: SMM
    """

    #basin_dict = {}
    max_basin = max(basin_order_list)
    basin_rename_list = [0] * max_basin
    basin_rename_list.append(0)

    print(("length is: "+str(len(basin_rename_list))))

    # Swap the keys
    for idx,basin in enumerate(basin_order_list):
        print(("Index: "+str(idx)+", basin: "+str(basin)))
        basin_rename_list[basin] = idx

    #print basin_rename_list
    return basin_rename_list


##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## This function orders basins in a sequence, so that plots can be made
## as a function of distance, for example
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def BasinOrderer(thisPointData, FileName, criteria_string,reverse=False,
                 exclude_criteria_string = "None", exlude_criteria_greater = False,
                 exclude_criteria_value = 0 ):

    #========================================================
    # now we need to label the basins
    # Now we get the chi points
    EPSG_string = LSDMap_IO.GetUTMEPSG(FileName)
    print("EPSG string is: " + EPSG_string)

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
    print(wd)

    if reverse:
        indices = indices[::-1]

    sorted_basins = np.array(these_data)[indices]


    print(sorted_basins)


#==============================================================================
# This function takes groups of data and then resets values in a
# raster to mimic these values
#==============================================================================
def RedefineIntRaster(rasterArray,grouped_data_list,spread):
    """This function takes values from an integer raster and renames them based on a list.

    It is useful for renaming basin numbers.

    Args:
        rasterArray (np.array): The raster array
        grouped_data_list (int): A list of lists containing groups to be redefined
        spread (int): How big of a difference between groups. For plotting this helps to generate differenc colours.

    Returns:
        np.array: The new array

    Author: SMM
    """

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

#==============================================================================
# This function takes groups of data and then resets values in a
# raster to mimic these values
#==============================================================================
def MaskByCategory(rasterArray,rasterForMasking,data_list):
    """This function takes values from an integer raster and renames them based on a list.

    It is useful for renaming basin numbers.

    Args:
        rasterArray (np.array): The raster array
        grouped_data_list (int): A list of lists containing groups to be redefined
        spread (int): How big of a difference between groups. For plotting this helps to generate differenc colours.

    Returns:
        np.array: The new array

    Author: SMM
    """

    # The -9090 is just a placeholder
    for item in data_list:
        rasterForMasking[rasterForMasking == item] = -9090

    rasterArray[rasterForMasking != -9090] = np.nan


    return rasterArray

#==============================================================================
# This function takes groups of data and then resets values in a
# raster to mimic these values
#==============================================================================
def NanBelowThreshold(rasterArray,threshold):
    """This function takes an array and turns any element below threshold to a nan

    It is useful for renaming basin numbers.

    Args:
        rasterArray (np.array): The raster array
        threshold (int): The threshold value

    Returns:
        np.array: The new array

    Author: SMM
    """

    # The -9090 is just a placeholder
    rasterArray[rasterArray < threshold] = np.nan

    return rasterArray



#==============================================================================
# This does a basic mass balance.
# Assumes all units are metres
#==============================================================================
def RasterMeanValue(path, file1):
    """This takes the average of a raster.

    Args:
        path (str): The path to the raster
        file1 (str): The name of the file

    Returns:
        mean_value: The mean

    Author: SMM
    """

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
    """This function averages all the data along one of the directions

    Args:
        path (str): The path to the files
        file1 (str): The name of the first raster.
        axis (int): Either 0 (rows) or 1 (cols)

    Returns:
        float: A load of information about the swath.

        * means
        * medians
        * std_deviations
        * twentyfifth_percentile
        * seventyfifth_percentile

        at each node across the axis of the swath.

    Author: SMM
    """

    # make sure names are in correct format
    NewPath = LSDOst.AppendSepToDirectoryPath(path)

    raster_file1 = NewPath+file1

    # get some information about the raster
    NDV, xsize, ysize, GeoT, Projection, DataType = LSDMap_IO.GetGeoInfo(raster_file1)

    print("NDV is: ")
    print(NDV)

    if NDV == None:
        NDV = -9999
        print("No NDV defined")

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
    """This function checks the difference in "volume" between two rasters.

    Args:
        path (str): The path to the files
        file1 (str): The name of the first raster.
        file2 (str): The name of the second raster

    Returns:
        float: The differnece in the volume betweeen the two rasters

    Author: SMM
    """

    # make sure names are in correct format
    NewPath = LSDOst.AppendSepToDirectoryPath(path)

    raster_file1 = NewPath+file1
    raster_file2 = NewPath+file2

    PixelArea = LSDMap_IO.GetPixelArea(raster_file1)
    print("PixelArea is: " + str(PixelArea))

    print("The formatted path is: " + NewPath)
    Raster1 = LSDMap_IO.ReadRasterArrayBlocks(raster_file1,raster_band=1)
    Raster2 = LSDMap_IO.ReadRasterArrayBlocks(raster_file2,raster_band=1)

    NewRaster = np.subtract(Raster2,Raster1)

    mass_balance = np.sum(NewRaster)*PixelArea

    print("linear dif " + str(np.sum(NewRaster)))

    return mass_balance
