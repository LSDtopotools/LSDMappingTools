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


def Rasterize_BGS_geologic_maps(shapefile_name):

    # The shapefile to be rasterized:     
    print('Rasterize ' + shapefile_name) 
    #get path and filename seperately 
    shapefilefilepath = LSDPT.GetPath(shapefile_name)
    shapefilename = LSDPT.GetFileNameNoPath(shapefile_name)
    shapefileshortname = LSDPT.GetFilePrefix(shapefile_name)

    print("Shapefile name is: "+shapefilename)

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
    outraster_bil = shapefilefilepath+"os.sep"+ outrastername_bil
    print("Full name of out raster is: "+outraster_bil)
    
    # This assumes UTM zone 30, because why would we do any work in East Anglia?
    system_call2 = 'gdalwarp -t_srs EPSG:32630 -of ENVI -dstnodata -9999 ' +  outraster + ' ' + outraster_bil
    os.system(system_call2)
    
    # Now get rid of the tif
    system_call3 = 'rm '+ outraster
    os.system(system_call3)
    
    # now get the the fields from the shapefile
    daShapefile = shapefile_name

    dataSource = ogr.Open(daShapefile)
    daLayer = dataSource.GetLayer(0)
 
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

if __name__ == "__main__":
    
    shapefile_name = '/home/smudd/SMMDataStore/analysis_for_papers/Geology_raster/bgs-50k_1726879/sc034/sc034_eyemouth_bedrock.shp'
    
    #shapefile_name = 'T:\\analysis_for_papers\\Geology_raster\\bgs-50k_1726879\\sc034\\sc034_eyemouth_bedrock.shp' 
    
    Rasterize_BGS_geologic_maps(shapefile_name)