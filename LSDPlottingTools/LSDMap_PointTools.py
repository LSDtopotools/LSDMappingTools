## LSDMap_Points.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## These functions are tools to deal with point files
## These files come in csv and can be read so that they can be output as
## Shapefiles or GeoJSON files
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## SMM
## 26/07/2014
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from __future__ import absolute_import, division, print_function, unicode_literals

from osgeo import osr
from . import LSDMap_OSystemTools as LSDOst
import os
import glob
import pandas
import numpy as np
from pyproj import Proj, transform


#==============================================================================
# This function takes all the csv files in a directory and converts to
# GeoJSON files
#==============================================================================
def ConvertAllCSVToGeoJSON(path):
    """This looks in a directory and converts all .csv files to GeoJSON.

    This is handy if, for example, you want to display data on the web using leaflet or D3.js

    Note:
        This assumes your csv files have latitude and longitude columns. If the LSDMap_PointData object will not be able to read them.

    Args:
        path (str): The path in which you want to convert the csv files

    Returns:
        None, but you will get a load of GeoJSON files.

    Author: SMM
    """

    # make sure names are in correct format
    NewPath = LSDOst.AppendSepToDirectoryPath(path)

    print("The formatted path is: " + NewPath)


    for FileName in glob(NewPath+"*.csv"):
        print("filename is: " + FileName)

        thisPointData = LSDMap_PointData(FileName)
        thisPointData.TranslateToReducedGeoJSON(FileName)


#==============================================================================
# This function takes all the csv files in a directory and converts to
# Shapefiles files
#==============================================================================
def ConvertAllCSVToShapefile(path):
    """This looks in a directory and converts all .csv files to shapefiles

    This is handy if, for example, you want to display data using ArcMap of QGIS

    Note:
        This assumes your csv files have latitude and longitude columns. If the LSDMap_PointData object will not be able to read them.

    Args:
        path (str): The path in which you want to convert the csv files

    Returns:
        None, but you will get a load of GeoJSON files.

    Author: SMM
    """

    # make sure names are in correct format
    NewPath = LSDOst.AppendSepToDirectoryPath(path)

    print("The formatted path is: " + NewPath)


    for FileName in glob(NewPath+"*.csv"):
        print("filename is: " + FileName)

        thisPointData = LSDMap_PointData(FileName)
        thisPointData.TranslateToReducedShapefile(FileName)


class LSDMap_PointData(object):

    # The constructor: it needs a filename to read
    def __init__(self,FileName):
        """This is the LSDMap_pointdata object. It loads csv files that have latitude and longitude data (in WGS84) and keeps other data records.

        The object can convert to UTM, and it also can print data to other file formats like GeoJSON and shapefiles.

        Args:
            Filename (str): The name of the csv file (with path and extension) that contains the point data. It should have columns labelled "latitude" and "longitude".

        Author: SMM
        """

        # This gets the filename without the .csv
        file_prefix = LSDOst.GetFilePrefix(FileName)

        self.FilePrefix = file_prefix
        print("The object file prefix is: " + self.FilePrefix)

        #See if the parameter files exist
        native_way = True
        if os.access(FileName,os.F_OK):
            if(native_way):
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


                print("Variable list is: ")
                print(self.VariableList)

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
                print("I am loading data using pandas, I haven't been widely tested yet, you can switch to the old way if you have troubles by looking for native_way in LSDMap_PointData")
                #Loading the file
                print("Loading")
                data = pandas.read_csv(FileName, sep=",")
                print("Loaded")
                #Extracting the headers
                self.VariableList = list(data.columns.values)
                #Correcting the names
                for names in self.VariableList:
                    names =  LSDOst.RemoveEscapeCharacters(names).lower()

                print("Your Variable list is : ")
                print(self.VariableList)
                #ingesting data values
                dataPoints = np.array(data.values)
                DataDict = {}
                #Feeding the dictionnary[header][table of values]
                print("I am ingesting the data")
                for i in range(len(self.VariableList)):
                    DataDict[self.VariableList[i]] = dataPoints[:,i]
                #dealing with the types
                DataDictTyped = {}
                TypeList = []
                for name in self.VariableList:
                    this_list = DataDict[name]
                    typed_list = LSDOst.ParseListToType(this_list)
                    DataDictTyped[name] = typed_list

                    TypeList.append(type(typed_list[0]))

                self.PointData = DataDictTyped
                self.DataTypes = TypeList
                print("The points data are successfully loaded")

        else:
            print("Uh oh I could not open that file")
            self.VariableList = []
            self.DataTypes = []
            self.PointData = {}

        # now make sure the data has latitude and longitude entries
        if "latitude" not in self.VariableList:
            print("Something has gone wrong, latitude is not in the variable list")
            print("Here is the variable list: ")
            print(self.VariableList)
        if "longitude" not in self.VariableList:
            print("Something has gone wrong, longitude is not in the variable list")
            print("Here is the variable list: ")
            print(self.VariableList)

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
        """Gets the list of parameter names.

        Args:
            PrintToScreen (bool): If true, prints to screen

        Return:
            str: A list of the variable names

        Author: SMM
        """

        if PrintToScreen:
            print(self.VariableList)

        return self.VariableList

    # Get data types
    def GetParameterTypes(self,PrintToScreen = False):
        """Gets the tyes of each names.

        Args:
            PrintToScreen (bool): If true, prints to screen

        Return:
            str: A list of the variable types

        Author: SMM
        """

        if PrintToScreen:
            print(self.DataTypes)

        return self.DataTypes


    # Get data elements
    def GetLatitude(self,PrintToScreen = False):
        """Gets the latitude list.

        Args:
            PrintToScreen (bool): If true, prints to screen

        Return:
            float: A list of the latitudes

        Author: SMM
        """


        if PrintToScreen:
            print(self.Latitude)

        return self.Latitude

    # Get data elements
    def GetLongitude(self,PrintToScreen = False):
        """Gets the longitude list.

        Args:
            PrintToScreen (bool): If true, prints to screen

        Return:
            float: A list of the longitudes

        Author: SMM
        """


        if PrintToScreen:
            print(self.Longitude)

        return self.Longitude

    def QueryData(self,data_name,PrintToScreen = False):

        """Returns the list of the data that has the column header data_name

        Args:
            PrintToScreen (bool): If true, prints to screen.
            data_name (str): The header of the column you want

        Return:
            float: A list of the data

        Author: SMM
        """
        if data_name not in self.VariableList:
            print("The data " + data_name + " is not one of the data elements in this point data")
            empty_list = []
            return empty_list
        else:

            if PrintToScreen:
                print("The " + data_name + "data is: ")
                print(self.PointData[data_name])

            return self.PointData[data_name]

    def GetUTMEastingNorthing(self,EPSG_string):
        """Returns two lists: the latitude and longitude converted to northing and easting.

        Args:
            PrintToScreen (bool): If true, prints to screen.
            EPSG_string (str): The EPSG code of the UTM coordinates you want (326XX) with zone XX is for north, 327XX is for south.

        Return:
            float: Two lists containing easting and northing

        Author: SMM
        """
        print("Yo, getting this stuff: "+EPSG_string)
        # The lat long are in epsg 4326 which is WGS84
        inProj = Proj(init='epsg:4326')
        outProj = Proj(init=EPSG_string)
        this_Lat = self.Latitude[0]
        this_Lon = self.Longitude[0]

        print("Lat-long: ")
        print(this_Lat)
        print(this_Lon)

        easting =[]
        northing = []

        for idx, Lon in enumerate(self.Longitude):
            Lat = self.Latitude[idx]
            ea,no = transform(inProj,outProj,Lon,Lat)
            easting.append(ea)
            northing.append(no)

        return easting,northing

    def GetUTMEastingNorthingFromQuery(self,EPSG_string,Latitude_string,Longitude_string):
        """Returns two lists: the latitude and longitude converted to northing and easting. But you can define the columns if there are more than one latitude and longitude columns.

        Note:
            This is used mainly if there are multple lat-long coordinates in the csv file. For example when you have basin centroids and basin outlets in the same file.
        Args:
            PrintToScreen (bool): If true, prints to screen.
            EPSG_string (str): The EPSG code of the UTM coordinates you want (326XX) with zone XX is for north, 327XX is for south.
            Latitude_string (str): The name of the latitude column you want
            Longitude_string (str): The name of the longitude column you want.

        Return:
            float: Two lists containing easting and northing

        Author: SMM
        """
        print("Yo, getting this stuff: "+EPSG_string)
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
        """This removes data from a point function that is below a threshold value

        Args:
            data_name (str): The name of the data member to select
            Threshold_value (float): Below this threshold points will be removed.

        Returns:
            None removes data from the object (not reversible!!)

        Author: SMM

        """
        print("I am thinning the data for you!")

        # Get the data for thinning
        if data_name not in self.VariableList:
            print("The data " + data_name + " is not one of the data elements in this point data")
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
## Data manipulation
##==============================================================================
##==============================================================================
    def ThinDataSelection(self,data_name,data_for_selection_list):
        """This function takes a list of values and retains the members in data name corresponding to that selection

        Args:
            data_name (str): The name of the data member to select
            data_for_selection_list (int): A list of values to retain. Useful for things like selecting basins or sources.

        Returns:
            None removes data from the object (not reversible!!)

        Author: SMM

        """


        print("I am thinning the data for you from a list!")

        # Get the data for thinning
        if data_name not in self.VariableList:
            print("The data " + data_name + " is not one of the data elements in this point data")
        else:
            this_data = self.PointData[data_name]

        this_data = [int(x) for x in this_data]
        #print("The original data I need to thin is: ")
        #print(this_data)


        # Start a new data dict
        NewDataDict = {}
        NewLat = []
        NewLon = []
        for name in self.VariableList:
            NewDataDict[name] = []


        # Get all the data to be deleted
        print("The data I am keeping is: ")
        print(data_for_selection_list)
        delete_indices = []
        for index, data in enumerate(this_data):
            #print("Data: "+str(data))
            if data not in data_for_selection_list:
                #print("I'm not keeping it")
                delete_indices.append(index)
            else:
                #print("I'll have that one. ")
                NewLat.append(self.Latitude[index])
                NewLon.append(self.Longitude[index])
                for name in self.VariableList:
                    this_element = self.PointData[name][index]
                    NewDataDict[name].append(this_element)


        # Now reset the data dict
        self.PointData = NewDataDict
        self.Latitude = NewLat
        self.Longitude = NewLon

        #print("The updated data is:")
        #print(self.PointData[data_name])


##==============================================================================
##==============================================================================
## Format conversion
##==============================================================================
##==============================================================================
    # This translates the CRNData object to an Esri shapefile
    def TranslateToReducedShapefile(self,FileName):
        """This converts the point data to a shapefile

        Args:
            FileName (str): the name of the file to be printed. The code strips the extension and turns it into .shp, so you can give it the name of the csv file and ti will still work.

        Return:
            None, but prints a new shapefile

        Author: SMM
        """

        import osgeo.ogr as ogr

        #  set up the shapefile driver
        driver = ogr.GetDriverByName("ESRI Shapefile")

        # Get the path to the file
        this_path = LSDOst.GetPath(FileName)
        DataName = self.FilePrefix

        FileOut = this_path+DataName+".shp"

        print("The filename will be: " + FileOut)

        # delete the existing file
        if os.path.exists(FileOut):
            driver.DeleteDataSource(FileOut)
            print("That file exists, I am deleting it in order to start again.")
        else:
            print("I am making a new shapefile for you")

        # create the data source
        data_source = driver.CreateDataSource(FileOut)

        # create the spatial reference, in this case WGS84 (which is ESPG 4326)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)

        print("Creating the layer")

        # create the layer
        layer = data_source.CreateLayer(DataName, srs, ogr.wkbPoint)

        # Add the field names
        for index,name in enumerate(self.VariableList):
            print("The variable name is " + name + " and the type is: " + str(self.DataTypes[index]))


            if self.DataTypes[index] is int:
                layer.CreateField(ogr.FieldDefn(name, ogr.OFTInteger))
            elif self.DataTypes[index] is float:
                layer.CreateField(ogr.FieldDefn(name, ogr.OFTReal))
            elif self.DataTypes[index] is str:
                print("Making a sting layer for layer " + name)
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
        """This converts the point data to a GeoJSON

        Args:
            FileName (str): the name of the file to be printed. The code strips the extension and turns it into .geojson, so you can give it the name of the csv file and ti will still work.

        Return:
            None, but prints a new GeoJSON

        Author: SMM
        """

        # Parse a delimited text file of volcano data and create a shapefile

        import osgeo.ogr as ogr
        import osgeo.osr as osr


        #  set up the shapefile driver
        driver = ogr.GetDriverByName("GeoJSON")

        # Get the path to the file
        this_path = LSDOst.GetPath(FileName)
        DataName = self.FilePrefix

        FileOut = this_path+DataName+".geojson"

        print("The filename will be: " + FileOut)

        # delete the existing file
        if os.path.exists(FileOut):
            driver.DeleteDataSource(FileOut)

        # create the data source
        data_source = driver.CreateDataSource(FileOut)

        # create the spatial reference,  in this case WGS84 (which is ESPG 4326)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)

        print("Creating the layer")

        # create the layer
        layer = data_source.CreateLayer("PointData", srs, ogr.wkbPoint)

        #print "Adding the field names"

        # Add the field names
        for index,name in enumerate(self.VariableList):
            print("The variable name is " + name + " and the type is: " + str(self.DataTypes[index]))


            if self.DataTypes[index] is int:
                layer.CreateField(ogr.FieldDefn(name, ogr.OFTInteger))
            elif self.DataTypes[index] is float:
                layer.CreateField(ogr.FieldDefn(name, ogr.OFTReal))
            elif self.DataTypes[index] is str:
                print("Making a sting layer for layer " + name)
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
