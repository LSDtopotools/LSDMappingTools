## LSDMap_GDALIO.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## These functions are tools to deal with rasters
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## SMM
## 26/07/2014
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#from __future__ import absolute_import, division, print_function, unicode_literals
from __future__ import absolute_import, division, print_function

import osgeo.gdal as gdal
import osgeo.gdal_array as gdal_array
import numpy as np
from osgeo import osr
from osgeo import ogr
import os
from os.path import exists
from osgeo.gdalconst import GA_ReadOnly

#==============================================================================
def getNoDataValue(rasterfn):
    """This gets the nodata value from the raster

    Args:
        rasterfn (str): The filename (with path and extension) of the raster

    Returns:
        float: nodatavalue; the nodata value

    Author: SMM
    """
    raster = gdal.Open(rasterfn)
    band = raster.GetRasterBand(1)
    return band.GetNoDataValue()
#==============================================================================

#==============================================================================
def setNoDataValue(rasterfn):
    """This sets the nodata value from the raster

    Args:
        rasterfn (str): The filename (with path and extension) of the raster

    Returns:
        None

    Author: SMM
    """

    raster = gdal.Open(rasterfn)
    band = raster.GetRasterBand(1)
    return band.SetNoDataValue()
#==============================================================================

#==============================================================================
def GetUTMMaxMin(FileName):
    """This gets the minimum and maximum UTM values.

    *WARNING* it assumes raster is already projected into UTM, and is in ENVI format! It reads from an ENVI header file.

    Args:
        FileName (str): The filename (with path and extension) of the raster

    Returns:
        float: The cell size in metres
        float: The X minimum (easting) in metres
        float: The X maximum (easting) in metres
        float: The Y minimum (northing) in metres
        float: The Y maximum (northing) in metres

    Author: SMM
    """


    if exists(FileName) is False:
            raise Exception('[Errno 2] No such file or directory: \'' + FileName + '\'')

    NDV, xsize, ysize, GeoT, Projection, DataType = GetGeoInfo(FileName)
    CellSize = GeoT[1]
    XMin = GeoT[0]
    XMax = XMin+CellSize*xsize

    YMax = GeoT[3]
    YMin = YMax-CellSize*ysize

    return CellSize,XMin,XMax,YMin,YMax
#==============================================================================

#==============================================================================
# Gets the pixel area, assumes units are projected
#==============================================================================
def GetPixelArea(FileName):
    """Gets the area in m^2 of the pixels

    Args:
        rasterfn (str): The filename (with path and extension) of the raster

    Returns:
        float: Pixel_area (float): The area of each pixel

    Author: SMM
    """

    if exists(FileName) is False:
            raise Exception('[Errno 2] No such file or directory: \'' + FileName + '\'')

    NDV, xsize, ysize, GeoT, Projection, DataType = GetGeoInfo(FileName)
    CellSize = GeoT[1]

    return CellSize*CellSize
#==============================================================================


#==============================================================================
# this takes rows and columns of minium and maximum values and converts them
# to UTM
def GetUTMMaxMinFromRowsCol(FileName,x_max_col,x_min_col,y_max_row,y_min_row):
    """This gets the minimum and maximum UTM values but you give it the row and column numbers.

    Note:
        This assumes raster is already projected into UTM, and is in ENVI format! It reads from an ENVI header file.

    Args:
        FileName (str): The filename (with path and extension) of the raster
        x_max_col (int): The column to use as the maximum
        x_min_col (int): The column to use as the minimum
        y_max_row (int): The row to use as the maximum
        y_min_row (int): The row to use as the minimum

    Returns:

        float: The X maximum (easting) in metres
        float: The X minimum (easting) in metres
        float: The Y maximum (northing) in metres
        float: The Y minimum (northing) in metres

    Author: SMM
    """


    if exists(FileName) is False:
            raise Exception('[Errno 2] No such file or directory: \'' + FileName + '\'')

    NDV, xsize, ysize, GeoT, Projection, DataType = GetGeoInfo(FileName)
    CellSize = GeoT[1]
    XMin = GeoT[0]

    YMax = GeoT[3]
    YMin = YMax-CellSize*ysize

    xmax_UTM = XMin+x_max_col*CellSize
    xmin_UTM = XMin+x_min_col*CellSize

    # need to be careful with the ymax_UTM since the rows go from the top
    # but the header index is to bottom corner

    print("yll: "+str(YMin)+" and nrows: " +str(ysize) + " dx: "+str(CellSize))

    ymax_from_bottom = ysize-y_min_row
    ymin_from_bottom = ysize-y_max_row
    ymax_UTM = YMin+ymax_from_bottom*CellSize
    ymin_UTM = YMin+ymin_from_bottom*CellSize

    return xmax_UTM,xmin_UTM,ymax_UTM,ymin_UTM
#==============================================================================

#==============================================================================
# This gets the x and y vectors of the data
#==============================================================================
def GetLocationVectors(FileName):
    """This gets a vector of the x and y locations of the coordinates

    Note:
        This assumes raster is already projected into UTM, and is in ENVI format! It reads from an ENVI header file.

    Args:
        FileName (str): The filename (with path and extension) of the raster.

    Return:
        float: A vector of the x locations (eastings)
        float: A vector of the y locations (northings)

    Author: SMM
    """



    NDV, xsize, ysize, GeoT, Projection, DataType = GetGeoInfo(FileName)

    CellSize,XMin,XMax,YMin,YMax = GetUTMMaxMin(FileName)



    x_vec = np.arange(XMin,XMax,CellSize)
    y_vec = np.arange(YMin,YMax,CellSize)

    return x_vec,y_vec
#==============================================================================




#==============================================================================
# This gets the extent of the raster
def GetRasterExtent(FileName):
    """This gets a vector of the minimums and maximums of the coordinates

    Note:
        This assumes raster is already projected into UTM, and is in ENVI format! It reads from an ENVI header file.

    Args:
        FileName (str): The filename (with path and extension) of the raster.

    Return:
        float: A vector that contains

            * extent[0]: XMin
            * extent[1]: XMax
            * extent[2]: YMin
            * extent[3]: YMax

    Author: SMM
    """
    CellSize,XMin,XMax,YMin,YMax = GetUTMMaxMin(FileName)
    extent = [XMin,XMax,YMin,YMax]
    return extent

#==============================================================================
# Function to read the original file's projection:
def GetGeoInfo(FileName):
    """This gets information from the raster file using gdal

    Args:
        FileName (str): The filename (with path and extension) of the raster.

    Return:
        float: A vector that contains:
            * NDV: the nodata values
            * xsize: cellsize in x direction
            * ysize: cellsize in y direction
            * GeoT: the tranform (a string)
            * Projection: the Projection (a string)
            * DataType: The type of data (an int explaing the bits of each data element)

    Author: SMM
    """


    if exists(FileName) is False:
            raise Exception('[Errno 2] No such file or directory: \'' + FileName + '\'')


    SourceDS = gdal.Open(FileName, gdal.GA_ReadOnly)
    if SourceDS == None:
        raise Exception("Unable to read the data file")

    NDV = SourceDS.GetRasterBand(1).GetNoDataValue()
    xsize = SourceDS.RasterXSize
    ysize = SourceDS.RasterYSize
    GeoT = SourceDS.GetGeoTransform()
    Projection = osr.SpatialReference()
    Projection.ImportFromWkt(SourceDS.GetProjectionRef())
    DataType = SourceDS.GetRasterBand(1).DataType
    DataType = gdal.GetDataTypeName(DataType)

    return NDV, xsize, ysize, GeoT, Projection, DataType
#==============================================================================

#==============================================================================
# This gets the UTM zone, if it exists
def GetUTMEPSG(FileName):
    """Uses GDAL to get the EPSG string from the raster.

    Args:
        FileName (str): The filename (with path and extension) of the raster.

    Return:
        str: The EPSG string

    Author: SMM
    """
    if exists(FileName) is False:
            raise Exception('[Errno 2] No such file or directory: \'' + FileName + '\'')

    # see if the file exists and get the dataset
    SourceDS = gdal.Open(FileName, gdal.GA_ReadOnly)
    if SourceDS == None:
        raise Exception("Unable to read the data file")

    EPSG_string = 'NULL'

    # get the projection
    print("Let me get that projection for you")
    print("In this function I will extract the UTM zone")
    prj=SourceDS.GetProjection()
    srs=osr.SpatialReference(wkt=prj)

    if srs.IsProjected:
        print("The dataset is projected.")
        #print("Trying projcs")
        #print(str(srs.GetAttrValue(str('PROJCS'),0)))
        #print(srs.GetAttrValue(str('projcs')))
        proj_str = srs.GetAttrValue(str('projcs'))
        print("The projection string is: "+proj_str)

        print(proj_str)


        if proj_str != None:
            
            

            # extract the UTM information
            if "UTM Zone" in proj_str:
                print("Found the string UTM zone")

                first_split = proj_str.split(',')
                first_half = first_split[0]
                second_half = first_split[1]
                if "Northern" in second_half:
                    N_or_S = "N"
                else:
                    N_or_S = "S"
                second_split = first_half.split(' ')
                zone = second_split[2]
                    
                                   
            elif "UTM zone" in proj_str:
  
                print("This seems to be from the new gdal version")
                first_split = proj_str.split(' ')
                zone_str = first_split[-1]
                print("Zone string is: "+zone_str)
                zone = zone_str[:-1]
                print("The zone is: "+zone)

                if zone_str[-1]=="S":
                    N_or_S = "S"
                else:
                    N_or_S = "N"    

                print("And the hemisphere is: "+N_or_S)


            elif "_Hemisphere" in proj_str:
                if "Southern" in proj_str:
                    N_or_S = "S"
                else:
                    N_or_S = "N"
                proj_split = proj_str.split('_')
                zone = proj_split[2]
                
            
            else:

                proj_split = proj_str.split('_')
                zone = proj_split[-1]

                N_or_S = zone[-1]
                zone = zone[:-1]


            # adding some logic for zones < 10
            if len(zone) < 2:
                zone = '0'+zone
                
            EPSG_string = 'epsg:'
            if N_or_S == 'S':
                EPSG_string = EPSG_string+'327'+zone
            else:
                EPSG_string = EPSG_string+'326'+zone
            print("The EPSG string is: "+EPSG_string)
    else:
        raise Exception("This is not a projected coordinate system!")



    print(EPSG_string)
    return EPSG_string


#==============================================================================
# Function to read the original file's projection:
def GetNPixelsInRaster(FileName):
    """This gets the total number of pixels in the raster

    Args:
        FileName (str): The filename (with path and extension) of the raster.

    Return:
        int: The total number of pixels

    Author: SMM
    """

    NDV, xsize, ysize, GeoT, Projection, DataType = GetGeoInfo(FileName)

    return xsize*ysize

#==============================================================================

#==============================================================================
# Function to read the original file's projection:
def CheckNoData(FileName):
    """This looks through the head file of an ENVI raster and if it doesn't find the nodata line it rewrites the file to include the nodata line.

    Args:
        FileName (str): The filename (with path and extension) of the raster.

    Return:
        int: The total number of pixels (although what it is really doing is updating the header file. The return is just to check if it is working and yes I know this is stupid. )

    Author: SMM
    """


    if exists(FileName) is False:
        raise Exception('[Errno 2] No such file or directory: \'' + FileName + '\'')

    # read the file, and check if there is a no data value
    SourceDS = gdal.Open(FileName, gdal.GA_ReadOnly)
    if SourceDS == None:
        raise Exception("Unable to read the data file")
    NoDataValue = SourceDS.GetRasterBand(1).GetNoDataValue()

    print("In the check nodata routine. Nodata is: ")
    print(NoDataValue)

    if NoDataValue == None:
        print("This raster does not have no data. Updating the header file")
        header_name = FileName[:-4]
        header_name = header_name+".hdr"

        # read the header
        if exists(header_name) is False:
            raise Exception('[Errno 2] No such file or directory: \'' + header_name + '\'')
        else:
            this_file = open(header_name, 'r')
            lines = this_file.readlines()
            lines.append("data ignore value = -9999")
            this_file.close()

            this_file = open(header_name, 'w')
            for item in lines:
                this_file.write("%s" % item)        # no newline since a newline command character comes with the lines

            this_file.close()

    NDV, xsize, ysize, GeoT, Projection, DataType = GetGeoInfo(FileName)

    return xsize*ysize

#==============================================================================


#==============================================================================
def ReadRasterArrayBlocks(raster_file,raster_band=1):
    """This reads a raster file (from GDAL) into an array. The "blocks" bit makes it efficient.
    Args:
        FileName (str): The filename (with path and extension) of the raster.
        raster_band (int): the band of the raster (almost all uses with LSDTopoTools will have a 1 band raster)

    Return:
        np.array: A numpy array with the data from the raster.

    Author: SMM
    """


    if exists(raster_file) is False:
            raise Exception('[Errno 2] No such file or directory: \'' + raster_file + '\'')

    dataset = gdal.Open(raster_file, GA_ReadOnly )
    if dataset == None:
        raise Exception("Unable to read the data file")

    band = dataset.GetRasterBand(raster_band)
    NoDataValue = dataset.GetRasterBand(1).GetNoDataValue()

    block_sizes = band.GetBlockSize()
    x_block_size = block_sizes[0]
    y_block_size = block_sizes[1]

    #If the block y size is 1, as in a GeoTIFF image, the gradient can't be calculated,
    #so more than one block is used. In this case, using8 lines gives a similar
    #result as taking the whole array.
    if y_block_size < 8:
        y_block_size = 8

    xsize = band.XSize
    ysize = band.YSize

    print("xsize: " +str(xsize)+" and y size: " + str(ysize))

    max_value = band.GetMaximum()
    min_value = band.GetMinimum()

    # now initiate the array
    data_array = np.zeros((ysize,xsize))

    #print "data shape is: "
    #print data_array.shape

    if max_value == None or min_value == None:
        stats = band.GetStatistics(0, 1)
        max_value = stats[1]
        min_value = stats[0]

    for i in range(0, ysize, y_block_size):
        if i + y_block_size < ysize:
            rows = y_block_size
        else:
            rows = ysize - i

        for j in range(0, xsize, x_block_size):
            if j + x_block_size < xsize:
                cols = x_block_size
            else:
                cols = xsize - j

            # get the values for this block
            values = band.ReadAsArray(j, i, cols, rows)

            # move these values to the data array
            data_array[i:i+rows,j:j+cols] = values

    print("NoData is:", NoDataValue)
    if NoDataValue is not None:
        nodata_mask = data_array == NoDataValue
        data_array[nodata_mask] = np.nan

    return data_array
#==============================================================================

#==============================================================================
def ReadRasterArrayBlocks_numpy(raster_file,raster_band=1):
    """
    This reads a raster file into an array using numpy. The "blocks" bit makes it efficient.
    Args:
        FileName (str): The filename (with path and extension) of the raster.
        raster_band (int): the band of the raster (almost all uses with LSDTopoTools will have a 1 band raster) Doesn't work yet for more than one band.

    Return:
        np.array: A numpy array with the data from the raster.

    Author: SMM
    """

    print("I will now use numpy.fromfile to load your raster, I still need to be tested, If something goes wrong switch back to the classic method by switching NFF_opti = False. Classic method = unefficient.")
    if exists(raster_file) is False:
            raise Exception('[Errno 2] No such file or directory: \'' + raster_file + '\'')

    with open(raster_file[:-3]+"hdr","r") as hdr_file:
        print("I am opening your raster info")
        no_data_hdr = -9999 # setting the default noData in case there is no value setted in the hdr file
        for line in hdr_file:
            #testing the data type
            if(line[0:9] == "data type"):
                info_dtype = line[-2:]
                info_dtype = int(info_dtype)
            else:
                if(line[0:8] == "map info"):
                    info = line[12:-2]
                    info = info.split(",")
                    x_min = float(info[3])
                    y_max = float(info[4])
                    x_res = float(info[5])
                    y_res = float(info[6])
                    utm_zone = int(info[7])
                    utm_hemisphere = info[8]
                else:
                    if(line[0:7] == "samples"):
                        num_col = line.replace(" ","").split("=")[1]
                        print("there are " + str(num_col) + " columns")
                        num_col = int(num_col)
                    else:
                        if(line[0:5] == "lines"):
                            num_lines = line.replace(" ","").split("=")[1]
                            print("there are " + str(num_lines) + " lines")
                            num_lines = int(num_lines)
                        else:
                            if(line[0:17] == "data ignore value"):
                                no_data_hdr = line.replace(" ","").split("=")[1]
                                print("No data value is: " + str(no_data_hdr) )
                                no_data_hdr = float(no_data_hdr)


        #The type of data representation:
        #1 = Byte: 8-bit unsigned integer
        #2 = Integer: 16-bit signed integer
        #3 = Long: 32-bit signed integer
        #4 = Floating-point: 32-bit single-precision
        #5 = Double-precision: 64-bit double-precision floating-point
        #6 = Complex: Real-imaginary pair of single-precision floating-point
        #9 = Double-precision complex: Real-imaginary pair of double precision floating-point
        #12 = Unsigned integer: 16-bit
        #13 = Unsigned long integer: 32-bit
        #14 = 64-bit long integer (signed)
        #15 = 64-bit unsigned long integer (unsigned)

    if(info_dtype == 1):
        data_type = np.dtype('uint8')
    else:
        if(info_dtype == 2):
            data_type = np.dtype('int16')
        else:
            if(info_dtype == 3):
                data_type = np.dtype('int32')
            else:
                if(info_dtype == 4):
                    data_type = np.dtype('Float32')
                else:
                    if(info_dtype == 5):
                        data_type = np.dtype('Float64')
                    else:
                        if(info_dtype == 12):
                            data_type = np.dtype('uint16')
    print("your data type is "  + str(data_type))
    # Alright now loading and converting the data
    print("I am now ingesting your raster")
    data_array = np.fromfile(raster_file,data_type).reshape(num_lines,num_col)
    if(info_dtype in [1,2,3,12]):
        data_array = data_array.astype(float)
    print("I nailed it")
    x_max = x_min + x_res*num_col
    y_min = y_max - y_res*num_lines
    print("I am returning the raster array and info")
    NoDataValue = no_data_hdr
    if NoDataValue is not None:
        data_array[data_array == NoDataValue] = np.nan


    return data_array
#==============================================================================

#==============================================================================
def array2raster(rasterfn,newRasterfn,array,driver_name = "ENVI", noDataValue = -9999):
    """Takes an array and writes to a GDAL compatible raster. It needs another raster to map the dimensions.

    Args:
        FileName (str): The filename (with path and extension) of a raster that has the same dimensions as the raster to be written.
        newRasterfn (str): The filename (with path and extension) of the new raster.
        array (np.array): The array to be written
        driver_name (str): The type of raster to write. Default is ENVI since that is the LSDTOpoTools format
        noDataValue (float): The no data value

    Return:
        np.array: A numpy array with the data from the raster.

    Author: SMM
    """

    raster = gdal.Open(rasterfn)
    geotransform = raster.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    cols = raster.RasterXSize
    rows = raster.RasterYSize

    driver = gdal.GetDriverByName(driver_name)
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outRaster.GetRasterBand(1).SetNoDataValue( noDataValue )
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromWkt(raster.GetProjectionRef())
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()
#==============================================================================


def RasterDifference(RasterFile1, RasterFile2, raster_band=1, OutFileName="Test.outfile", OutFileType="ENVI"):
    """
    Takes two rasters of same size and subtracts second from first,
    e.g. Raster1 - Raster2 = raster_of_difference
    then writes it out to file
    """

    Raster1 = gdal.Open(RasterFile1)
    Raster2 = gdal.Open(RasterFile2)

    print("RASTER 1: ")
    print(Raster1.GetGeoTransform())
    print(Raster1.RasterCount)
    print(Raster1.GetRasterBand(1).XSize)
    print(Raster1.GetRasterBand(1).YSize)
    print(Raster1.GetRasterBand(1).DataType)

    print("RASTER 2: ")
    print(Raster2.GetGeoTransform())
    print(Raster2.RasterCount)
    print(Raster2.GetRasterBand(1).XSize)
    print(Raster2.GetRasterBand(1).YSize)
    print(Raster2.GetRasterBand(1).DataType)

    raster_array1 = np.array(Raster1.GetRasterBand(raster_band).ReadAsArray())
    raster_array2 = np.array(Raster2.GetRasterBand(raster_band).ReadAsArray())

    assert(raster_array1.shape == raster_array2.shape )
    print("Shapes: ", raster_array1.shape, raster_array2.shape)

    difference_raster_array = raster_array1 - raster_array2

#    import matplotlib.pyplot as plt
#
#    plt.imshow(difference_raster_array)
#

    driver = gdal.GetDriverByName(OutFileType)

    dsOut = driver.Create(OutFileName,
                          Raster1.GetRasterBand(1).XSize,
                          Raster1.GetRasterBand(1).YSize,
                          1,
                          gdal.GDT_Float32)
                          #Raster1.GetRasterBand(raster_band).DataType)


    gdal_array.CopyDatasetInfo(Raster1,dsOut)
    bandOut = dsOut.GetRasterBand(1)
    gdal_array.BandWriteArray(bandOut, difference_raster_array)

#==============================================================================
def PolygoniseRaster(DataDirectory, RasterFile, OutputShapefile='polygons'):
    """
    This function takes in a raster and converts to a polygon shapefile using rasterio
    from https://gis.stackexchange.com/questions/187877/how-to-polygonize-raster-to-shapely-polygons/187883#187883?newreg=8b1f507529724a8488ce4789ba787363

    Args:
        DataDirectory (str): the data directory with the basin raster
        RasterFile (str): the name of the raster
        OutputShapefile (str): the name of the output shapefile WITHOUT EXTENSION. Default = 'polygons'

    Returns:
        Dictionary where key is the raster value and the value is a shapely polygon

    Author: FJC
    """
    # import modules
    import rasterio
    from rasterio.features import shapes
    from shapely.geometry import shape, Polygon, mapping
    import fiona

    # define the mask
    #mask = None
    raster_band = 1

    # get raster no data value
    NDV = getNoDataValue(DataDirectory+RasterFile)

    # load in the raster using rasterio
    with rasterio.open(DataDirectory+RasterFile) as src:
        image = src.read(raster_band, masked=False)

        msk = src.read_masks(1)

        results = (
        {'properties': {'raster_val': v}, 'geometry': s}
        for i, (s, v)
        in enumerate(
            shapes(image, mask=msk, transform=src.transform)))

    # define shapefile attributes
    # crs = src.crs.wkt
    # print (crs)
    print("Let me grab the coordinate reference system.")
    crs = GetUTMEPSG(DataDirectory+RasterFile)
    schema = {'geometry': 'Polygon',
              'properties': { 'ID': 'float'}}



    # This is necessary to filter the basin results
    geoms = list(results)
    #print("Geom size is: "+str(len(geoms)))
    
    filtered_geoms = {}
    area_dict = {}
    for f in geoms:
        this_shape = Polygon(shape(f['geometry']))
        this_val = float(f['properties']['raster_val'])
        #print("ID is: "+str(this_val))
        this_area = this_shape.area
        if this_val in filtered_geoms.keys():
            print("Whoops. Found a repeated ID. Getting rid of the smaller one.")
            if area_dict[this_val] < this_area:
                filtered_geoms[this_val] = f
                area_dict[this_val] = this_area               
                print("Found a repeated ID. Keeping the one with area of "+str(this_area))
            else:
                print("Keeping the initial ID.")
        else:
            filtered_geoms[this_val] = f
            area_dict[this_val] = this_area
    
    new_geoms = []
    for key,item in filtered_geoms.items():
        this_shape = Polygon(shape(item['geometry']))
        this_val = float(item['properties']['raster_val'])
        #print("ID is: "+str(this_val)) 
        this_area = this_shape.area
        #print("Area is: "+str(this_area))
        new_geoms.append(item)
    #print("Geom size is: "+str(len(new_geoms)))
            
    # transform results into shapely geometries and write to shapefile using fiona
    PolygonDict = {}
    with fiona.open(DataDirectory+OutputShapefile, 'w', crs=crs, driver='ESRI Shapefile', schema=schema) as output:
        for f in new_geoms:
            this_shape = Polygon(shape(f['geometry']))
            this_val = float(f['properties']['raster_val'])
            print("ID is: "+str(this_val))
            if this_val != NDV: # remove no data values
                output.write({'geometry': mapping(this_shape), 'properties':{'ID': this_val}})
            PolygonDict[this_val] = this_shape

    return PolygonDict

#==============================================================================
def PolygoniseRasterMerge(DataDirectory, RasterFile, OutputShapefile='polygons'):
    """
    This function takes in a raster and converts to a polygon shapefile using rasterio
    from https://gis.stackexchange.com/questions/187877/how-to-polygonize-raster-to-shapely-polygons/187883#187883?newreg=8b1f507529724a8488ce4789ba787363

    This version recognises where there are multiple polygons with the same key and merges
    them to a MultiPolygon using cascaded_union

    Args:
        DataDirectory (str): the data directory with the basin raster
        RasterFile (str): the name of the raster
        OutputShapefile (str): the name of the output shapefile WITHOUT EXTENSION. Default = 'polygons'

    Returns:
        Dictionary where key is the raster value and the value is a shapely polygon

    Author: FJC
    """
    # import modules
    import rasterio
    from rasterio.features import shapes
    from shapely.geometry import shape, Polygon, mapping
    from shapely.ops import cascaded_union
    import fiona

    # define the mask
    #mask = None
    raster_band = 1

    # get raster no data value
    NDV = getNoDataValue(DataDirectory+RasterFile)

    # load in the raster using rasterio
    with rasterio.open(DataDirectory+RasterFile) as src:
        image = src.read(raster_band, masked=False)

        msk = src.read_masks(1)

        results = (
        {'properties': {'raster_val': v}, 'geometry': s}
        for i, (s, v)
        in enumerate(
            shapes(image, mask=msk, transform=src.transform)))

    # define shapefile attributes
    # crs = src.crs.wkt
    # print (crs)
    crs = GetUTMEPSG(DataDirectory+RasterFile)
    schema = {'geometry': 'Polygon',
              'properties': { 'ID': 'float'}}



    # transform results into shapely geometries and write to shapefile using fiona
    geoms = list(results)
    PolygonDict = {}
    with fiona.open(DataDirectory+OutputShapefile, 'w', crs=crs, driver='ESRI Shapefile', schema=schema) as output:
        for f in geoms:
            this_shape = Polygon(shape(f['geometry']))
            this_val = float(f['properties']['raster_val'])
            if this_val in PolygonDict:
                Polygons = [this_shape, PolygonDict[this_val]]
                this_shape = cascaded_union(Polygons)
            if this_val != NDV: # remove no data values
                output.write({'geometry': mapping(this_shape), 'properties':{'ID': this_val}})

            PolygonDict[this_val] = this_shape

    return PolygonDict


def CreateShapefileOfRasterFootprint(DataDirectory, RasterFile):
    """
    This function takes a raster and creates a shapefile that is the footprint
    of the raster. Used for plotting the raster footprint on regional maps using
    basemap.
    Variously put together from:
    http://osgeo-org.1560.x6.nabble.com/gdal-dev-Creating-a-simple-shapefile-with-ogr-td3749101.html
    
    
    Args:
        DataDirectory (str): the data directory with the basin raster
        RasterFile (str): the name of the raster

    Returns:
        Shapefile of the raster footprint

    Author: SMM
    
    Date: 23/01/2018
    """

    print("Trying to create a shapefile.")
    print("The Data directory is: "+DataDirectory+ " and the raster is: "+ RasterFile)
    driver_name = "ESRI shapefile"
    driver = ogr.GetDriverByName(driver_name)

    # get the filename of the outfile.
    if not DataDirectory.endswith(os.sep):
        print("You forgot the separator at the end of the directory, appending...")
        DataDirectory = DataDirectory+os.sep
        
    # Get the raster prefix
    SplitRasterfile = RasterFile.split(".")
    RasterPrefix = SplitRasterfile[0]
 
    # get the espg of the raster
    FullFilename = DataDirectory+RasterFile
    ESPG_this_raster = GetUTMEPSG(FullFilename)  
    ESPG_this_raster = str(ESPG_this_raster)
    print("The raster has coordinate of: "+ESPG_this_raster)
    ESPG_this_raster_split = ESPG_this_raster.split(":")
    ESPG_this_raster = ESPG_this_raster_split[-1] 
    print ("This ESPG is: "+str(ESPG_this_raster))
    
    # Get extent of raster
    [xmin,xmax,ymin,ymax] = GetRasterExtent(FullFilename)

    # Create ring
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(xmin, ymin)
    ring.AddPoint(xmin, ymax)
    ring.AddPoint(xmax, ymax)
    ring.AddPoint(xmax, ymin)
    ring.AddPoint(xmin, ymin)

    # Create polygon
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)    
    
    # Create a coordinate transformation
    source = osr.SpatialReference()
    source.ImportFromEPSG(int(ESPG_this_raster))
    
    target = osr.SpatialReference()
    target.ImportFromEPSG(4326)
    
    transform = osr.CoordinateTransformation(source, target)
    
    # now transformt the polygon
    poly.Transform(transform)
    
    # see what you got
    #print("The polygon is:")
    #print(poly.ExportToWkt()) 

    # create the data source
    OutFileName = DataDirectory+RasterPrefix+"_footprint.shp" 
    print("The output shapefile is: "+OutFileName)
    datasource = driver.CreateDataSource(OutFileName)    


    # create the layer
    layer = datasource.CreateLayer(OutFileName, target, ogr.wkbPolygon)
    feature = ogr.Feature(layer.GetLayerDefn())   
    feature.SetGeometry(poly)
    layer.CreateFeature(feature)
    
    # Clean up
    feature.Destroy()
    datasource.Destroy()
    
    
def GetCentreAndExtentOfRaster(DataDirectory, RasterFile):
    """
    This function takes a raster and returns the centrepoint and the extent in both degrees and metres. 
    
    Args:
        DataDirectory (str): the data directory with the basin raster
        RasterFile (str): the name of the raster

    Returns:
        The lat-long of the centrepoint, the x-y- extent in both degrees and metres 

    Author: SMM
    
    Date: 01/02/2018
    """

    print("Trying to create a shapefile.")
    print("The Data directory is: "+DataDirectory+ " and the raster is: "+ RasterFile)
    driver_name = "ESRI shapefile"
    driver = ogr.GetDriverByName(driver_name)

    # get the filename of the outfile.
    if not DataDirectory.endswith(os.sep):
        print("You forgot the separator at the end of the directory, appending...")
        DataDirectory = DataDirectory+os.sep
        
    # Get the raster prefix
    SplitRasterfile = RasterFile.split(".")
    RasterPrefix = SplitRasterfile[0]
 
    # get the espg of the raster
    FullFilename = DataDirectory+RasterFile
    ESPG_this_raster = GetUTMEPSG(FullFilename)  
    ESPG_this_raster = str(ESPG_this_raster)
    print("The raster has coordinate of: "+ESPG_this_raster)
    ESPG_this_raster_split = ESPG_this_raster.split(":")
    ESPG_this_raster = ESPG_this_raster_split[-1] 
    print ("This ESPG is: "+str(ESPG_this_raster))
    
    # Get extent of raster
    [xmin,xmax,ymin,ymax] = GetRasterExtent(FullFilename)
    xproj_extent = xmax-xmin
    yproj_extent = ymax-ymin
        
    # Create ring
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(xmin, ymin)
    ring.AddPoint(xmin, ymax)
    ring.AddPoint(xmax, ymax)
    ring.AddPoint(xmax, ymin)
    
    # Create a coordinate transformation
    source = osr.SpatialReference()
    source.ImportFromEPSG(int(ESPG_this_raster))
    
    target = osr.SpatialReference()
    target.ImportFromEPSG(4326)
    
    transform = osr.CoordinateTransformation(source, target)
    
    # now transform the ring so you can get coordinates in lat-long
    ring.Transform(transform)
    
    # now get the xmin,ymin, and xmax, ymax coords in lat-long
    pt1 = ring.GetPoint(0)
    min_long = pt1[0]
    min_lat = pt1[1]
    
    pt2 = ring.GetPoint(2)
    max_long = pt2[0]
    max_lat = pt2[1]
    
    extent_long = max_long-min_long
    extent_lat = max_lat-min_lat
    
    centre_long = min_long+extent_long*0.5
    centre_lat = min_lat+extent_lat*0.5
    
    return centre_lat, centre_long, extent_lat, extent_long, xproj_extent, yproj_extent
    
    
    
    # Leaving this here to show how to loop through points
    # Now try with the geometry tools
    #print("Trying ring")
    #for i in range(0, ring.GetPointCount()):
    #    pt = ring.GetPoint(i)
    #    print("The point is: ")
    #    print(str(pt[0])+" "+str(pt[1]))
    