#-*- coding: utf-8 -*-
"""
Script that take the xyz output of Lithochild (Ask Emma for more info) and convert it to BIL raster.

Created Today
@author: Boris
"""
import pandas
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import gdal, ogr, os, osr


########## Ignore this ########
def array2raster(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,array, nodata, EPSG):
    """This function take a regular array to create a raster with it"""
    print("I am dealing with nodata values")
    array[np.isnan(array)] = nodata # Dealing with Nodata values
    print("I am writing the raster")
    cols = array.shape[1]
    rows = array.shape[0]

    originX = rasterOrigin[0]
    originY = rasterOrigin[1]

    driver = gdal.GetDriverByName('ENVI')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Float64)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    #outband.SetNoDataValue(nodata)
    outband.WriteArray(array)
    #outband.SetNoDataValue(nodata)

    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(EPSG)

    outRaster.SetProjection(outRasterSRS.ExportToWkt())

    outband.FlushCache()
##############################################################



if __name__ == "__main__": # Ignore this line

    ########## Parameters ########
    # Parameters I need

    print(" This script should take the output of lithochild and convert it to a bil ENVI raster ready for use with LSDTopoTools/MappingTool")

    Directory = "/home/s1675537/PhD/DataStoreBoris/Emma/" # needed if the file is not in the same directory as the script
    csv_file = "beticstest5_ts11_xyz.txt" # name of the test file containing the data
    separator = ", " # string that separate the data /!\ in case, tabulation is \t /!\
    raster_name = "Final.bil"
    EPSG = 32630 # EPSG code of your projection
    write_Directory = "/home/s1675537/PhD/DataStoreBoris/Emma/" # the path were you want the file to be placed


    Xres = 10 # resolution on the X axis
    Yres = 10 # resolution for the Y axis (I haven't tested it yet on a non-square matrix)
    Xmin = 560488 # Xmin of your raster
    Xmax = 605518 # Xmax of your raster
    Ymin = 4086007 # Ymin of your raster
    Ymax = 4120537 # Ymax of your raster
    nodata = -1 # No data value

    ########## Python code, no more parameters to set ########
    rasterOrigin = (Xmin, Ymax) # The origin of the raster is the upper left corner
    raster_output = write_Directory+raster_name
    hdr_name = raster_name[:-4]+".hdr" # Needed to deal with no data

    # Importing the files

    print("I am importing the data from csv with Pandas, ignore the warning if there is any")
    df = pandas.read_csv(Directory + csv_file, sep=separator) # Pandas is impressively fast
    print("I am done, I am  now converting it into workable array")
    X = np.linspace(0,Xmax-Xmin,(Xmax-Xmin)/Xres) # creating the extent of the mesh
    Y = np.linspace(0,Ymax-Ymin,(Ymax-Ymin)/Yres)
    gridx,gridy = np.meshgrid(X,Y) # creating the mesh
    data = np.array(df.values) # Numpyisation of the data
    print("I have everything, I am now converting this irregular grid to a regular one %s/%s"%(Xres,Yres))
    grid_z0 = griddata(data[:,0:2], data[:,2],(gridx,gridy) , method='linear') # interpolation, Linear look the more "natural"
    print("the minimum value is %s and the maximum is %s"%(np.nanmin(grid_z0),np.nanmax(grid_z0)))

    print("I am now creating and saving the raster in EPSG %s" %(EPSG))
    grid_z0 = grid_z0[::-1] # inverting the array, required to make the raster
    array2raster(raster_output, rasterOrigin, Xres,Yres, grid_z0, nodata,EPSG) # Creating the raster
    print("OK, now I have to add nodata value at the hand of the hdr file because Gdal does not do it for some reason")
    #data ignore value = 0
    with open(write_Directory+hdr_name, "a") as myfile:
        myfile.write("data ignore value = " +str(nodata))

    print("I am done. If your raster is empty (nan or -1 values) it might means that your X/Y coord./res. are not right. If the problem persists, well, try other things or ask Simon")
