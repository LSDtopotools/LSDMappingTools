## LSDMap_VectorTools.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## These functions are tools to deal with vector data using shapely
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## FJC
## 26/06/17
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from . import LSDMap_GDALIO as LSDMap_IO
from shapely.geometry import Point, Polygon
import os
from os.path import exists
from osgeo import ogr, osr
import LSDPlottingTools as LSDPT
import gdal as gdal
from osgeo.gdalconst import GA_ReadOnly
from LSDMapFigure import PlottingHelpers as Helper

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
# BASIN FUNCTIONS
# These functions do various operations on basin polygons
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
def GetBasinOutlines(DataDirectory, basins_fname):
	"""
	This function takes in the raster of basins and gets a dict of basin polygons,
	where the key is the basin key and the value is a shapely polygon of the basin.

	IMPORTANT: In this case the "basin key" is usually the junction number:
		this function will use the raster values as keys and in general
		the basin rasters are output based on junction indices rather than keys

	Args:
		DataDirectory (str): the data directory with the basin raster
		basins_fname (str): the basin raster

	Returns:
		list of shapely polygons with the basins

	Author: FJC
	"""
	# read in the basins raster
	this_fname = basins_fname.split('.')
	print(basins_fname)
	OutputShapefile = this_fname[0]+'.shp'

	# polygonise the raster
	BasinDict = LSDMap_IO.PolygoniseRaster(DataDirectory, basins_fname, OutputShapefile)
	return BasinDict

def GetMultipleBasinOutlines(DataDirectory):
  """
	This function takes in multiple rasters of basins and gets a dict of basin polygons,
	where the key is the basin key derived from the file name and the value is a shapely polygon of the basin.

	IMPORTANT: In this case the "basin key" is usually the junction number:
		this function will use the raster values as keys and in general
		the basin rasters are output based on junction indices rather than keys

	Args:
		DataDirectory (str): the data directory with the basin raster

	Returns:
		list of shapely polygons with the basins

	Author: MDH
	"""
  # get a list of basins and declare the dictionary to populate
  basin_dict = Helper.MapBasinsToKeys(DataDirectory)
  BasinsDict = {}

  #loop across the basins
  for outlet_jn, basin_key in basin_dict.iteritems():
    this_fname = "basin"+str(outlet_jn)+"_AllBasins.bil"

    TempBasins = GetBasinOutlines(DataDirectory,this_fname)

    for temp_outlet, temp_basin_key in TempBasins.iteritems():
      if len(TempBasins) > 1:
        print("WARNING: MULTIPLE BASINS IN basin #", outlet_jn)
      TempBasins[int(outlet_jn)] = TempBasins.pop(temp_outlet)

    BasinsDict.update(TempBasins)

  return BasinsDict

def GetBasinCentroids(DataDirectory, basins_fname):
	"""
	This function takes in the raster of basins and returns a dict where the
	key is the basin key and the value is the shapely point of the centroid

	In most cases the "basin key" is actually the junction index: it comes
	from the basins labeled within the basin raster, which is output with
	junction indices rather than junction keys

	Args:
		DataDirectory (str): the data directory with the basin raster
		fname_prefix (str): the prefix for the DEM

	Returns:
		dict of centroid points

	Author: FJC
	"""
	# get the basin polygons
	BasinDict = GetBasinOutlines(DataDirectory, basins_fname)

	# get the centroids
	CentroidDict = {}
	for basin_key, basin in BasinDict.iteritems():
		CentroidDict[basin_key] = Point(basin.centroid)

	return CentroidDict

def GetPointWithinBasins(DataDirectory,basins_fname):
	"""
	This function takes in the raster of basin and returns a dict where the
	key is the basin key and the value is a shapely point that is representative
	of the basin (guaranteed to be within the polygon)

	In most cases the "basin key" is actually the junction index: it comes
	from the basins labeled within the basin raster, which is output with
	junction indices rather than junction keys

	Args:
		DataDirectory (str): the data directory with the basin raster
		fname_prefix (str): the prefix for the DEM

	Returns:
		dict of representative points

	Author: FJC
	"""
	# get the basin polygons
	BasinDict = GetBasinOutlines(DataDirectory, basins_fname)

	# get the centroids
	PointDict = {}
	for basin_key, basin in BasinDict.iteritems():
		PointDict[basin_key] = Point(basin.representative_point())

	return PointDict

def GetPointsWithinMultipleBasins(DataDirectory,basins_fname):
  """
  This function takes in rasters of basins and returns a dict where the
  key is the basin key and the value is a shapely point that is representative
  of the basin (guaranteed to be within the polygon)

  In most cases the "basin key" is actually the junction index: it comes
  from the basins labeled within the basin raster, which is output with
  junction indices rather than junction keys

  Args:
	  DataDirectory (str): the data directory with the basin raster
	  fname_prefix (str): the prefix for the DEM

  Returns:
	  dict of representative points

  Author: FJC
  """
  # get the basin polygons
  BasinDict = GetMultipleBasinOutlines(DataDirectory)
  print("BASIN DICT IS")
  print(BasinDict)

  # get the centroids
  PointDict = {}
  for basin_key, basin in BasinDict.iteritems():
    PointDict[basin_key] = Point(basin.representative_point())

  print("POINT DICT IS")
  print(PointDict)

  return PointDict

def GetPointWithinBasinsBuffered(DataDirectory,basins_fname, basin_list = [], buffer_frac=0.1):
	"""
	This function takes in the raster of basins, and buffers each basin
	(makes each one smaller). It then gets the centroid of each buffered
	basin and returns as a dict where the key is the basin key and the value
	is a shapely point that is the centroid of the buffered basin.

	In most cases the "basin key" is actually the junction index: it comes
	from the basins labeled within the basin raster, which is output with
	junction indices rather than junction keys

	This doesn't work at the moment - need to think of a way to specify the buffer
	distance appropriately.

	Args:
		DataDirectory (str): the data directory with the basin raster
		fname_prefix (str): the prefix for the DEM
		buffer_frac (float): the fraction of the basin to be removed by the
		buffer, default = 0.1

	Returns:
		dict of representative points

	Author: FJC
	"""
	# get the basin polygons
	BasinDict = GetBasinOutlines(DataDirectory, basins_fname)

	# buffer and get the centre of the buffered polygons
	PointDict = {}
	for basin_key, basin in BasinDict.iteritems():
		# get the x and y lengths of the basin and append to list
		print("This basin key is: "+str(basin_key))
		lengths = []
		bounds = basin.bounds
		lengths.append(bounds[2] - bounds[0])
		lengths.append(bounds[3] - bounds[1])
		print(min(lengths))

		# buffer with a fraction of the minimum length
		new_basin = Polygon(basin.buffer(min(lengths)*buffer_frac*-1))

		# get the centroid of the buffered basin
		PointDict[basin_key] = Point(new_basin.centroid)

	return PointDict


##### This part is copied from the LSD_GeologyTools.py file to make the functions accessible from another scripts and thus easier to ingest, it will be cleaned up at some points.
##### the aim of those functions is to raterize a lithologic raster

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
	format = "GTiff".encode('utf-8')
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
	rasterDS =  gdal.GetDriverByName('GTiff'.encode('utf-8')).Create(outraster, xRes, yRes, 1,  gdal.GDT_Byte)

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
	geol_field = geol_field.encode('utf-8')
	for feature in daLayer:
		ID = feature.GetField(geol_field)
		GEOL = feature.GetField("GEOL_CODE".encode('utf-8'))

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
	new_field = ogr.FieldDefn("GEOL_CODE".encode('utf-8'), ogr.OFTInteger)
	daLayer.CreateField(new_field)

	# lets see what the layers are
	print("Let me tell you what the names of the fields are after I added one!")
	layerDefinition = daLayer.GetLayerDefn()
	for i in range(layerDefinition.GetFieldCount()):
		print(layerDefinition.GetFieldDefn(i).GetName())


	# Make a key for the bedrock
	geol_dict = dict()
	geol_iterator = 0
	geol_field = geol_field.encode('utf-8')
	for feature in daLayer:
		GEOL = feature.GetField(geol_field)

		if GEOL not in geol_dict:
			geol_iterator = geol_iterator+1
			print("I found a new rock type, GEOL: "+ str(GEOL)+ " and rock type: " + str(geol_iterator))
			geol_dict[GEOL] = geol_iterator

		# now get the geol code
		this_geol_code = geol_dict[GEOL]
		# set the feature
		feature.SetField("GEOL_CODE".encode('utf-8'), this_geol_code)

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
	cliffbaggu = "ESRI shapefile"
	cliffbaggu = cliffbaggu.encode('utf-8')
	driver = ogr.GetDriverByName(cliffbaggu)
	#src.Destroy()


	# Now write to the the outfile
	out_ds = driver.CreateDataSource(new_shapefile_name)
	# create the output layer
	#out_lyr = out_ds.CreateLayer("yo",srs = daLayer.GetSpatialRef(),geom_type=ogr.wkbPolygon)
	out_lyr = out_ds.CreateLayer("yo".encode('utf-8'),srs = daLayer.GetSpatialRef(),geom_type=geom_type)

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


def rasterize_shapefile(path_to_shp, res = 30, field = ""):
	"""
		I am going to lead the rasterization of the shapefile into a tif raster.
		I'll work on the bil version at some points.

		@param:
			path_to_shapefile (str) the path and name of the shapefile

		@returns: Nothing but write a raster
		@Author: BG
		@date: 28/09/2017

	"""

	print("I will raterize your shapefile:")
	shapefile_name = path_to_shp
	print(shapefile_name)

	#launching the  rasterization
	new_shapefile_name, geol_dict = geologic_maps_modify_shapefile(shapefile_name, geol_field = field)
	tifname = Rasterize_geologic_maps_pythonic(new_shapefile_name,raster_resolution = res, geol_field = field)
	Correct_Raterized_GLIM_map(tifname)

	print("Now removing the temporary files")
	os.remove(new_shapefile_name)
	os.remove(tifname)
	print("done with the rasterization")
