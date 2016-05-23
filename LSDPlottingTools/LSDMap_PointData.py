## LSDMap_Points.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## These functions are tools to deal with point files
## These files come in csv and can be read so that they can be output as
## Shapefiles or GeoJSON files
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
import LSDOSystemTools as LSDOst
import os

class LSDMap_PointData(object):
    
    # The constructor: it needs a filename to read    
    def __init__(self,FileName):
        
        # This gets the filename without the .csv
        file_prefix = LSDOst.GetFilePrefix(FileName)
        
        self.FilePrefix = file_prefix        
        print "The object file prefix is: " + self.FilePrefix
        
        #See if the parameter files exist
        if os.access(FileName,os.F_OK):
            this_file = open(FileName, 'r')
            lines = this_file.readlines()

            # get rid of the control characters
            this_line = LSDOst.RemoveEscapeCharacters(lines[0])
            
            # Now get a list with the names of the parameters
            self.VariableList = []                        
            TestList = this_line.split(',')
            
            for name in TestList:
                self.VariableList.append(name.lower())
                
            
            print "Variable list is: "
            print self.VariableList
            
            # get rid of the names
            del lines[0]             
            
            # now you need to make a dict that contains a list for each varaible name
            DataDict = {}  
            for name in self.VariableList:
                DataDict[name] = []
                
            # now get the data into the dict
            for line in lines:
                this_line = LSDOst.RemoveEscapeCharacters(line)
                split_line = this_line.split(',')
                
                for index,name in enumerate(self.VariableList):
                    DataDict[name].append(float(split_line[index]))                    
          
            self.PointData = DataDict
        else:
            print "Uh oh I could not open that file"
            self.VariableList = []
            self.PointData = {}
            
        # now make sure the data has latitude and longitude entries
        if "latitude" not in self.VariableList:
            print "Something has gone wrong, latitude is not in the variable list"
            print "Here is the variable list: "
            print self.VariableList
        if "longitude" not in self.VariableList:
            print "Something has gone wrong, longitude is not in the variable list"
            print "Here is the variable list: "
            print self.VariableList            

        # Add the latitude and longitude to their own data members and get rid 
        # of those from the VariableList        
        self.Latitude = self.PointData["latitude"]
        self.Longitude = self.PointData["longitude"]
        
 

##==============================================================================
##==============================================================================
## DATA ACCESS
##==============================================================================
##==============================================================================       
    # Get data elements
    def GetParameterNames(self,PrintToScreen = False):
        
        if PrintToScreen:        
            print self.VariableList
            
        return self.VariableList 
        
    # Get data elements
    def GetLatitude(self,PrintToScreen = False):
        
        if PrintToScreen:        
            print self.Latitude
            
        return self.Latitude     
        
    # Get data elements
    def GetLongitude(self,PrintToScreen = False):
        
        if PrintToScreen:        
            print self.Longitude
            
        return self.Longitude            
        
##==============================================================================
##==============================================================================
## Format conversion
##==============================================================================
##==============================================================================         
    # This translates the CRNData object to an Esri shapefile
    def TranslateToReducedShapefile(self,FileName):

        import osgeo.ogr as ogr

        #  set up the shapefile driver
        driver = ogr.GetDriverByName("ESRI Shapefile")

        # Get the path to the file
        this_path = LSDOst.GetPath(FileName)
        DataName = self.FilePrefix
        
        FileOut = this_path+DataName+".shp"
        
        print "The filename will be: " + FileOut

        # delete the existing file
        if os.path.exists(FileOut):
            driver.DeleteDataSource(FileOut)
            print "That file exists, I am deleting it in order to start again."
        else:
            print "I am making a new shapefile for you"

        # create the data source
        data_source = driver.CreateDataSource(FileOut)
        
        # create the spatial reference, in this case WGS84 (which is ESPG 4326)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)

        print "Creating the layer"

        # create the layer
        layer = data_source.CreateLayer(DataName, srs, ogr.wkbPoint)

        print "Adding the field names"
        
        # Add the field names
        for name in self.VariableList:
            layer.CreateField(ogr.FieldDefn(name, ogr.OFTReal))
        
        # Process the text file and add the attributes and features to the shapefile
        for index,lat in enumerate(self.Latitude):
            
            # create the feature
            feature = ogr.Feature(layer.GetLayerDefn())
            
            for name in self.VariableList:
                feature.SetField(name, self.PointData[name][index])    

            # create the WKT for the feature using Python string formatting
            wkt = "POINT(%f %f)" %  (float(self.Longitude[index]), float(self.Latitude[index]))

            # Create the point from the Well Known Txt
            point = ogr.CreateGeometryFromWkt(wkt)

            # Set the feature geometry using the point
            feature.SetGeometry(point)
            # Create the feature in the layer (shapefile)
            layer.CreateFeature(feature)
            # Destroy the feature to free resources
            feature.Destroy()

        # Destroy the data source to free resources
        data_source.Destroy()        
        
        
    # This translates the CRNData object to an GeoJSON
    def TranslateToReducedGeoJSON(self,FileName):
        # Parse a delimited text file of volcano data and create a shapefile

        import osgeo.ogr as ogr
        import osgeo.osr as osr


        #  set up the shapefile driver
        driver = ogr.GetDriverByName("GeoJSON")

        # Get the path to the file
        this_path = LSDOst.GetPath(FileName)
        DataName = self.FilePrefix
        
        FileOut = this_path+DataName+".geojson"
        
        print "The filename will be: " + FileOut

        # delete the existing file
        if os.path.exists(FileOut):
            driver.DeleteDataSource(FileOut)

        # create the data source
        data_source = driver.CreateDataSource(FileOut)
        
        # create the spatial reference,  in this case WGS84 (which is ESPG 4326)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)

        print "Creating the layer"
        
        # create the layer
        layer = data_source.CreateLayer("PointData", srs, ogr.wkbPoint)

        print "Adding the field names"
        
        # Add the field names
        for name in self.VariableList:
            layer.CreateField(ogr.FieldDefn(name, ogr.OFTReal))
        
        # Process the text file and add the attributes and features to the shapefile
        for index,lat in enumerate(self.Latitude):
            
            # create the feature
            feature = ogr.Feature(layer.GetLayerDefn())
            
            for name in self.VariableList:
                feature.SetField(name, self.PointData[name][index])   

            # create the WKT for the feature using Python string formatting
            wkt = "POINT(%f %f)" %  (float(self.Longitude[index]), float(self.Latitude[index]))

            # Create the point from the Well Known Txt
            point = ogr.CreateGeometryFromWkt(wkt)

            # Set the feature geometry using the point
            feature.SetGeometry(point)
            # Create the feature in the layer (shapefile)
            layer.CreateFeature(feature)
            # Destroy the feature to free resources
            feature.Destroy()

        # Destroy the data source to free resources
        data_source.Destroy()        