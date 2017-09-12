# -*- coding: utf-8 -*-
"""
Created on Tue Feb 07 13:34:02 2017

@author: smudd
"""

# This script modified from
# http://geoinformaticstutorial.blogspot.co.uk/2012/11/convert-shapefile-to-raster-with-gdal.html

# Importing needed modules
import os
from os.path import exists
from osgeo import ogr, osr
import LSDPlottingTools as LSDPT
import gdal as gdal
import numpy as np
from osgeo.gdalconst import GA_ReadOnly

def readFile(filename):
    print("Hey buddy, Reading the file: "+filename)

    filehandle = gdal.Open(filename, GA_ReadOnly )
    if filehandle == None:
        raise Exception("Unable to read the data file")

    band1 = filehandle.GetRasterBand(1)
    geotransform = filehandle.GetGeoTransform()
    geoproj = filehandle.GetProjection()
    Z = band1.ReadAsArray()
    xsize = filehandle.RasterXSize
    ysize = filehandle.RasterYSize
    return xsize,ysize,geotransform,geoproj,Z

def writeFile(filename,geotransform,geoprojection,data):
    (x,y) = data.shape
    format = "GTiff"
    noDataValue = -9999
    driver = gdal.GetDriverByName(format)
    # you can change the dataformat but be sure to be able to store negative values including -9999
    dst_datatype = gdal.GDT_Float32

    #print(data)

    dst_ds = driver.Create(filename,y,x,1,dst_datatype)
    dst_ds.GetRasterBand(1).WriteArray(data)
    dst_ds.GetRasterBand(1).SetNoDataValue( noDataValue )
    dst_ds.SetGeoTransform(geotransform)
    dst_ds.SetProjection(geoprojection)
    return 1



def Rasterize_BGS_geologic_maps(shapefile_name):

    # The shapefile to be rasterized:
    print('Rasterize ' + shapefile_name)
    #get path and filename seperately
    shapefilefilepath = LSDPT.GetPath(shapefile_name)
    shapefilename = LSDPT.GetFileNameNoPath(shapefile_name)
    shapefileshortname = LSDPT.GetFilePrefix(shapefile_name)

    print("Shapefile name is: "+shapefilename)

    # now get the the fields from the shapefile
    daShapefile = shapefile_name

    dataSource = ogr.Open(daShapefile)
    daLayer = dataSource.GetLayer(0)

    # lets see what the layers are
    print("Let me tell you what the names of the fields are!")
    layerDefinition = daLayer.GetLayerDefn()
    for i in range(layerDefinition.GetFieldCount()):
        print(layerDefinition.GetFieldDefn(i).GetName())


    # The raster file to be created and receive the rasterized shapefile
    outrastername = shapefileshortname + '.tif'
    outraster = shapefilefilepath+os.sep+ outrastername
    outcsv = shapefilefilepath+os.sep+shapefileshortname+'_lithokey.csv'
    print("Full name of out raster is: "+outraster)

    # Rasterize!!
    system_call = 'gdal_rasterize -a BGSREF -l ' + shapefileshortname +' -tr 90 -90 -a_nodata -9999 ' +  shapefile_name + ' ' + outraster
    print("System call is: ")
    print(system_call)
    os.system(system_call)

    # now convert the raster to UTM, as well as delete the stupid TIF
    # The raster file to be created and receive the rasterized shapefile
    outrastername_bil = shapefileshortname + '.bil'
    outraster_bil = shapefilefilepath+os.sep+ outrastername_bil
    print("Full name of out raster is: "+outraster_bil)

    # This assumes UTM zone 30, because why would we do any work in East Anglia?
    system_call2 = 'gdalwarp -t_srs EPSG:32630 -of ENVI -dstnodata -9999 ' +  outraster + ' ' + outraster_bil
    os.system(system_call2)

    # Now get rid of the tif
    system_call3 = 'rm '+ outraster
    os.system(system_call3)


    # Make a key for the bedrock
    geol_dict = dict()
    for feature in daLayer:
        ID = feature.GetField("BGSREF")
        GEOL = feature.GetField("RCS_D")

        if ID not in geol_dict:
            print("I found a new rock type, ID: "+ str(ID)+ " and rock type: " + str(GEOL))
            geol_dict[ID] = GEOL

    print("The rocks are: ")
    print(geol_dict)

    with open(outcsv, 'wb') as f:
        f.write('ID,rocktype\n')
        for key in geol_dict:
            f.write(str(key)+','+ str(geol_dict[key])+'\n')

    print("All done")


def Rasterize_geologic_maps_pythonic(shapefile_name, raster_resolution = 400, geol_field = "xx"):

    # The shapefile to be rasterized:
    print('Rasterize ' + shapefile_name)
    #get path and filename seperately
    shapefilefilepath = LSDPT.GetPath(shapefile_name)
    shapefilename = LSDPT.GetFileNameNoPath(shapefile_name)
    shapefileshortname = LSDPT.GetFilePrefix(shapefile_name)

    print("Shapefile name is: "+shapefilename)

    # now get the the fields from the shapefile
    daShapefile = shapefile_name

    dataSource = ogr.Open(daShapefile)
    daLayer = dataSource.GetLayer(0)

    # lets see what the layers are
    print("Let me tell you what the names of the fields are!")
    layerDefinition = daLayer.GetLayerDefn()
    for i in range(layerDefinition.GetFieldCount()):
        print(layerDefinition.GetFieldDefn(i).GetName())

    # The raster file to be created and receive the rasterized shapefile
    outrastername = shapefileshortname + '.tif'
    outraster = shapefilefilepath+os.sep+ outrastername
    outcsv = shapefilefilepath+os.sep+shapefileshortname+'_lithokey.csv'
    print("Full name of out raster is: "+outraster)

    # Create the destination data source
    inGridSize=float(raster_resolution)
    xMin, xMax, yMin, yMax = daLayer.GetExtent()

    xRes = int((xMax - xMin) / inGridSize)
    yRes = int((yMax - yMin) / inGridSize)
    rasterDS =  gdal.GetDriverByName('GTiff').Create(outraster, xRes, yRes, 1,  gdal.GDT_Byte)

    # Define spatial reference
    NoDataVal = -9999
    rasterDS.SetProjection(daLayer.GetSpatialRef().ExportToWkt())
    rasterDS.SetGeoTransform((xMin, inGridSize, 0, yMax, 0, -inGridSize))
    rBand = rasterDS.GetRasterBand(1)
    rBand.SetNoDataValue(NoDataVal)
    rBand.Fill(NoDataVal)

    # Rasterize
    gdal.RasterizeLayer(rasterDS, [1], daLayer, options = ["ATTRIBUTE=GEOL_CODE"])

    # Make a key for the bedrock
    geol_dict = dict()
    for feature in daLayer:
        ID = feature.GetField(geol_field)
        GEOL = feature.GetField("GEOL_CODE")

        if ID not in geol_dict:
            print("I found a new rock type, ID: "+ str(ID)+ " and rock type: " + str(GEOL))
            geol_dict[ID] = GEOL

    print("The rocks are: ")
    print(geol_dict)

    with open(outcsv, 'wb') as f:
        f.write('ID,rocktype\n')
        for key in geol_dict:
            f.write(str(key)+','+ str(geol_dict[key])+'\n')

    print("Done rasterizing!")
    return outraster

def Correct_Raterized_GLIM_map(tifname):
    # And now for a hack that converts to
    print("The raster name is: "+tifname)

    [xsize,ysize,geotransform,geoproj,Z] = readFile(tifname)

    print("Before data check")
    print(Z)

    print("Data type is: "+ str(Z.dtype))
    X = Z.astype(int)
    # Set large negative values to -9999
    X[X<=0] = -9999
    #Z[np.isnan(Z)]= -9999

    print("After_data_check")
    print(X)

    #get path and filename seperately
    filepath = LSDPT.GetPath(tifname)
    #filename = LSDPT.GetFileNameNoPath(tifname)
    fileshortname = LSDPT.GetFilePrefix(tifname)

    outraster2 = filepath+fileshortname + '2.tif'
    writeFile(outraster2,geotransform,geoproj,X)

def geologic_maps_modify_shapefile(shapefile_name, geol_field = "xx"):

    # The shapefile to be rasterized:
    print('Rasterize ' + shapefile_name)
    #get path and filename seperately
    shapefilefilepath = LSDPT.GetPath(shapefile_name)
    #shapefilename = LSDPT.GetFileNameNoPath(shapefile_name)
    shapefileshortname = LSDPT.GetFilePrefix(shapefile_name)

    # get the new shapefile name
    new_shapefile_name = shapefilefilepath+os.sep+shapefileshortname+"_new.shp"

    # copy the shapefile into the new shapefile--we don't wwant to mess up the original data
    print("The New Shapefile name is: "+new_shapefile_name)
    Copy_Shapefile(shapefile_name,new_shapefile_name)

    # New shapefile is opened for writing.
    dataSource = ogr.Open(new_shapefile_name,1)
    daLayer = dataSource.GetLayer(0)

    # add a new field
    new_field = ogr.FieldDefn("GEOL_CODE", ogr.OFTInteger)
    daLayer.CreateField(new_field)

    # lets see what the layers are
    print("Let me tell you what the names of the fields are after I added one!")
    layerDefinition = daLayer.GetLayerDefn()
    for i in range(layerDefinition.GetFieldCount()):
        print(layerDefinition.GetFieldDefn(i).GetName())


    # Make a key for the bedrock
    geol_dict = dict()
    geol_iterator = 0
    for feature in daLayer:
        GEOL = feature.GetField(geol_field)

        if GEOL not in geol_dict:
            geol_iterator = geol_iterator+1
            print("I found a new rock type, GEOL: "+ str(GEOL)+ " and rock type: " + str(geol_iterator))
            geol_dict[GEOL] = geol_iterator

        # now get the geol code
        this_geol_code = geol_dict[GEOL]
        # set the feature
        feature.SetField("GEOL_CODE", this_geol_code)

        # need to update the layer
        daLayer.SetFeature(feature)

    print("The rocks are: ")
    print(geol_dict)

    print("All done")


    return new_shapefile_name, geol_dict


def Copy_Shapefile(shapefile_name,new_shapefile_name):
    """
    Sweet Jesus why is this so difficult?
    """

    if exists(shapefile_name) is False:
        raise Exception('[Errno 2] No such file or directory: \'' + shapefile_name + '\'')

    # get the short name of the new shapefile
    shapefileshortname = LSDPT.GetFilePrefix(new_shapefile_name)
    print("The shortname is: "+shapefileshortname)

    # read in the data
    src = ogr.Open(shapefile_name)
    daLayer = src.GetLayer(0)

    # lets see what the layers are
    print("Let me tell you what the names of the fields are!")
    layerDefinition = daLayer.GetLayerDefn()
    for i in range(layerDefinition.GetFieldCount()):
        print(layerDefinition.GetFieldDefn(i).GetName())

    geom_type = layerDefinition.GetGeomType()

    # get rid of previous copies
    if exists(new_shapefile_name):
        os.remove(new_shapefile_name)

    # get the driver and create a new data source
    driver = ogr.GetDriverByName('ESRI Shapefile')
    #src.Destroy()


    # Now write to the the outfile
    out_ds = driver.CreateDataSource(new_shapefile_name)
    # create the output layer
    #out_lyr = out_ds.CreateLayer("yo",srs = daLayer.GetSpatialRef(),geom_type=ogr.wkbPolygon)
    out_lyr = out_ds.CreateLayer("yo",srs = daLayer.GetSpatialRef(),geom_type=geom_type)

    # Add input Layer Fields to the output Layer if it is the one we want
    for i in range(0, layerDefinition.GetFieldCount()):
        fieldDefn = layerDefinition.GetFieldDefn(i)
        #fieldName = fieldDefn.GetName()
        out_lyr.CreateField(fieldDefn)

    # Get the output Layer's Feature Definition
    outLayerDefn = out_lyr.GetLayerDefn()

    # Add features to the ouput Layer
    for inFeature in daLayer:
        # Create output Feature
        outFeature = ogr.Feature(outLayerDefn)

        # add in geometries
        geom = inFeature.GetGeometryRef()
        outFeature.SetGeometry(geom.Clone())
        # Add new feature to output Layer

        # Add field values from input Layer
        for i in range(0, outLayerDefn.GetFieldCount()):
            fieldDefn = outLayerDefn.GetFieldDefn(i)
            #fieldName = fieldDefn.GetName()
            outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(),
                inFeature.GetField(i))

        out_lyr.CreateFeature(outFeature)
        #out_ds.Destroy()



if __name__ == "__main__":

    #shapefile_name = '/home/smudd/SMMDataStore/analysis_for_papers/Geology_raster/bgs-50k_1726879/sc034/sc034_eyemouth_bedrock.shp'
    #shapefile_name = 'T:\\analysis_for_papers\\Geology_raster\\bgs-50k_1726879\\sc034\\sc034_eyemouth_bedrock.shp'

    #shapefile_name = '/home/smudd/SMMDataStore/analysis_for_papers/Iberia_geology/SouthernSpain_geology.shp'
    #tifname = '/home/smudd/SMMDataStore/analysis_for_papers/Iberia_geology/SouthernSpain_geology_new2.tif'

    shapefile_name = '/home/s1675537/PhD/LSDTopoData/movern_analysis/santa_cruz/GS.shp'


    #shapefile_name = 'C:\\VagrantBoxes\\LSDTopoTools\\Topographic_projects\\Iberia\\TipOfSpain.shp'
    #new_shapefile_name = 'C:\\VagrantBoxes\\LSDTopoTools\\Topographic_projects\\Iberia\\New_TipOfSpain.shp'
    #tifname = 'C:\\VagrantBoxes\\LSDTopoTools\\Topographic_projects\\Iberia\\TipOfSpain_new.tif'

    new_shapefile_name, geol_dict = geologic_maps_modify_shapefile(shapefile_name, geol_field = "LITH62")
    tifname = Rasterize_geologic_maps_pythonic(new_shapefile_name,raster_resolution = 30, geol_field = "LITH62")
    Correct_Raterized_GLIM_map(tifname)
