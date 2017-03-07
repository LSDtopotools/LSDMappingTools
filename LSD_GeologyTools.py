# -*- coding: utf-8 -*-
"""
Created on Tue Feb 07 13:34:02 2017

@author: smudd
"""

# This script modified from 
# http://geoinformaticstutorial.blogspot.co.uk/2012/11/convert-shapefile-to-raster-with-gdal.html

# Importing needed modules 
import os 
from osgeo import ogr
import LSDPlottingTools as LSDPT
import gdal as gdal


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

def Rasterize_GLIM_geologic_maps(shapefile_name):

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
    system_call = 'gdal_rasterize -a OBJECTID -l ' + shapefileshortname +' -tr 400 -400 -a_nodata -9999 ' +  shapefile_name + ' ' + outraster
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
    #system_call3 = 'rm '+ outraster
    #os.system(system_call3)
    
 
    # Make a key for the bedrock
    geol_dict = dict()
    for feature in daLayer:
        ID = feature.GetField("xx")
        GEOL = feature.GetField("xx")
        
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
    
def Rasterize_GLIM_geologic_maps_pythonic(shapefile_name):

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
    inGridSize=float(400)
    xMin, xMax, yMin, yMax = daLayer.GetExtent()
    xRes = int((xMax - xMin) / inGridSize)
    yRes = int((yMax - yMin) / inGridSize)
    rasterDS = gdal.GetDriverByName('GTiff').Create(outraster, xRes, yRes, 1,    gdal.GDT_Byte)
    
    # Define spatial reference
    NoDataVal = -9999
    rasterDS.SetProjection(daLayer.GetSpatialRef().ExportToWkt())
    rasterDS.SetGeoTransform((xMin, inGridSize, 0, yMax, 0, -inGridSize))
    rBand = rasterDS.GetRasterBand(1)
    rBand.SetNoDataValue(NoDataVal)
    rBand.Fill(NoDataVal)

    # Rasterize
    err = gdal.RasterizeLayer(rasterDS, [1], daLayer, burn_values=[200], options = ["ALL_TOUCHED=TRUE", "ATTRIBUTE=OBJECTID"])
      
    # Make a key for the bedrock
    geol_dict = dict()
    for feature in daLayer:
        ID = feature.GetField("xx")
        GEOL = feature.GetField("xx")
        
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

def GLIM_geologic_maps_modify_shapefile(shapefile_name):

    # The shapefile to be rasterized:     
    print('Rasterize ' + shapefile_name) 
    #get path and filename seperately 
    shapefilefilepath = LSDPT.GetPath(shapefile_name)
    shapefilename = LSDPT.GetFileNameNoPath(shapefile_name)
    shapefileshortname = LSDPT.GetFilePrefix(shapefile_name)

    print("Shapefile name is: "+shapefilename)

    # now get the the fields from the shapefile
    daShapefile = shapefile_name

    dataSource = ogr.Open(daShapefile,1)
    daLayer = dataSource.GetLayer(0)
    
    # lets see what the layers are
    print("Let me tell you what the names of the fields are!")
    layerDefinition = daLayer.GetLayerDefn()
    for i in range(layerDefinition.GetFieldCount()):
        print(layerDefinition.GetFieldDefn(i).GetName())        
        
    # add a new field
    new_field = ogr.FieldDefn("GEOL_CODE", ogr.OFTInteger)
    daLayer.CreateField(new_field)

    # lets see what the layers are
    print("Let me tell you what the names of the fields are after I added one!")
    layerDefinition = daLayer.GetLayerDefn()
    for i in range(layerDefinition.GetFieldCount()):
        print(layerDefinition.GetFieldDefn(i).GetName())      
    
      
    # Make a key for the bedrock
    #geol_dict = dict()
    #for feature in daLayer:
    #    ID = feature.GetField("xx")
    #    
    #    if ID == "va":
    #        feature.SetField("GEOL_CODE", 0)
    #    elif ID == "vb":
    #        feature.SetField("GEOL_CODE", 1)
    #    elif ID == "ss":
    #        feature.SetField("GEOL_CODE", 3)
    #    else:
    #        feature.SetField("GEOL_CODE", 4)
      
    print("All done")   

    
if __name__ == "__main__":
    
    #shapefile_name = '/home/smudd/SMMDataStore/analysis_for_papers/Geology_raster/bgs-50k_1726879/sc034/sc034_eyemouth_bedrock.shp'
    shapefile_name = '/home/smudd/SMMDataStore/analysis_for_papers/Iberia_geology/SouthernSpain_geology.shp'
    #shapefile_name = 'T:\\analysis_for_papers\\Geology_raster\\bgs-50k_1726879\\sc034\\sc034_eyemouth_bedrock.shp' 
    
    #Rasterize_BGS_geologic_maps(shapefile_name)
    #Rasterize_GLIM_geologic_maps_pythonic(shapefile_name)
    GLIM_geologic_maps_modify_shapefile(shapefile_name)