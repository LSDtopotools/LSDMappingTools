# -*- coding: utf-8 -*-
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
    print("I amDealing with nodata values")
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
    outband.WriteArray(array)

    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(EPSG)

    outRaster.SetProjection(outRasterSRS.ExportToWkt())

    outband.FlushCache()
##############################################################



if __name__ == "__main__": # Ignore this line

    ########## Parameters ########
    # Parameters I need

    print(" This script should take the output of lithochild and convert it to a bil ENVI raster ready for use with LSDTopoTools/MappingTool")

    Directory = "" # needed if the file is not in the same directory as the script
    csv_file = "beticstest5_ts11_xyz.txt" # name of the test file containing the data
    separator = "\t" # string that separate the data /!\ tabulation is \t /!\
    raster_name = "output.bil"
    EPSG = 32630 # EPSG code of your projection


    Xres = 10 # resolution on the X axis
    Yres = 10 # resolution for the Y axis (I haven't tested it yet on a non-square matrix)
    Xmin = 560488 # Xmin of your raster
    Xmax = 605518 # Xmax of your raster
    Ymin = 4086007 # Ymin of your raster
    Ymax = 4120537 # Ymax of your raster
    nodata = -1 # No data value

    ########## Python code, no more parameters to set ########
    rasterOrigin = (Xmin, Ymax) # The origin of the raster is the upper left corner

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
    array2raster(raster_name, rasterOrigin, Xres,Yres, grid_z0, nodata,EPSG) # Creating the raster

    print("I am done")
