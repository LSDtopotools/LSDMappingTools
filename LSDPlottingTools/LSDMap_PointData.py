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
from pyproj import Proj, transform

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
                this_name = LSDOst.RemoveEscapeCharacters(name)
                self.VariableList.append(this_name.lower())
                
            
            print "Variable list is: "
            print self.VariableList
            
            # get rid of the names
            del lines[0]             
            
            # now you need to make a dict that contains a list for each varaible name
            DataDict = {} 
            TypeList = []
            for name in self.VariableList:
                DataDict[name] = []
                
            # now get the data into the dict
            #firstline = True
            for line in lines:
                this_line = LSDOst.RemoveEscapeCharacters(line)
                split_line = this_line.split(',')
                
                for index,name in enumerate(self.VariableList):
                    this_var = LSDOst.RemoveEscapeCharacters(split_line[index])
                    #this_variable = LSDOst.ParseStringToType(this_var)
                    DataDict[name].append(this_var) 
            
            # now go back and get the correct type             
            DataDictTyped = {}    
            for name in self.VariableList:
                this_list = DataDict[name]
                typed_list = LSDOst.ParseListToType(this_list)
                DataDictTyped[name] = typed_list
                
                TypeList.append(type(typed_list[0])) 
                             
            self.PointData = DataDictTyped
            self.DataTypes = TypeList
        else:
            print "Uh oh I could not open that file"
            self.VariableList = []
            self.DataTypes = []
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

    # Get data types
    def GetParameterTypes(self,PrintToScreen = False):
        
        if PrintToScreen:        
            print self.DataTypes
            
        return self.DataTypes 

        
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

    def QueryData(self,data_name,PrintToScreen = False):

        if data_name not in self.VariableList:
            print "The data " + data_name + " is not one of the data elements in this point data"
        else:
        
            if PrintToScreen:
                print "The " + data_name + "data is: "
                print self.PointData[data_name]  
                
            return self.PointData[data_name]   
 
    def GetUTMEastingNorthing(self,EPSG_string):
        
        print "Yo, getting this stuff: "+EPSG_string
        # The lat long are in epsg 4326 which is WGS84
        inProj = Proj(init='epsg:4326')
        outProj = Proj(init=EPSG_string)
        this_Lat = self.Latitude[0]
        this_Lon = self.Longitude[0]
        
        print "Lat-long: "
        print this_Lat
        print this_Lon
        
        easting =[]
        northing = []
        
        for idx, Lon in enumerate(self.Longitude):
            Lat = self.Latitude[idx]
            ea,no = transform(inProj,outProj,Lon,Lat)
            easting.append(ea)
            northing.append(no)
            
        return easting,northing

    def GetUTMEastingNorthingFromQuery(self,EPSG_string,Latitude_string,Longitude_string):
        
        print "Yo, getting this stuff: "+EPSG_string
        # The lat long are in epsg 4326 which is WGS84
        inProj = Proj(init='epsg:4326')
        outProj = Proj(init=EPSG_string)
        
        
        this_Lat = self.QueryData(Latitude_string)
        this_Lon = self.QueryData(Longitude_string)

        
        easting =[]
        northing = []
        
        for idx, Lon in enumerate(this_Lon):
            Lat = this_Lat[idx]
            ea,no = transform(inProj,outProj,Lon,Lat)
            easting.append(ea)
            northing.append(no)
            
        return easting,northing

        
    
##==============================================================================
##==============================================================================
## Data manipulation
##==============================================================================
##==============================================================================    
    def ThinData(self,data_name,Threshold_value):
        
        print "I am thinning the data for you!"
        
        # Get the data for thinning
        if data_name not in self.VariableList:
            print "The data " + data_name + " is not one of the data elements in this point data"
        else:       
            this_data = self.PointData[data_name]
        
        this_data = [float(x) for x in this_data]
        
        # Start a new data dict
        NewDataDict = {}
        NewLat = []
        NewLon = []
        for name in self.VariableList:
            NewDataDict[name] = []

        
        # Get all the data to be delelted
        delete_indices = []
        for index, data in enumerate(this_data):
            if data<Threshold_value:
                delete_indices.append(index)
            else:
                NewLat.append(self.Latitude[index])
                NewLon.append(self.Longitude[index])
                for name in self.VariableList:
                    this_element = self.PointData[name][index]
                    NewDataDict[name].append(this_element)    
                
        # Now reset the data dict
        self.PointData = NewDataDict
        self.Latitude = NewLat
        self.Longitude = NewLon
    
    
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

        # Add the field names
        for index,name in enumerate(self.VariableList):
            print "The variable name is " + name + " and the type is: " + str(self.DataTypes[index])           
            
            
            if self.DataTypes[index] is int:           
                layer.CreateField(ogr.FieldDefn(name, ogr.OFTInteger))
            elif self.DataTypes[index] is float:
                layer.CreateField(ogr.FieldDefn(name, ogr.OFTReal))
            elif self.DataTypes[index] is str:
                print "Making a sting layer for layer " + name
                layer.CreateField(ogr.FieldDefn(name, ogr.OFTString))
            else:
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

        #print "Adding the field names"
        
        # Add the field names
        for index,name in enumerate(self.VariableList):
            print "The variable name is " + name + " and the type is: " + str(self.DataTypes[index])           
            
            
            if self.DataTypes[index] is int:           
                layer.CreateField(ogr.FieldDefn(name, ogr.OFTInteger))
            elif self.DataTypes[index] is float:
                layer.CreateField(ogr.FieldDefn(name, ogr.OFTReal))
            elif self.DataTypes[index] is str:
                print "Making a sting layer for layer " + name
                layer.CreateField(ogr.FieldDefn(name, ogr.OFTString))
            else:
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