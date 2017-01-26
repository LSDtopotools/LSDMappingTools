# ==================================
# Extracting basin-averaged values from a raster (.flt) given the junction number at the outlet of the basin (shape file). 
# Creates a csv file with the values.
# Configured to manipulate the zone directories and files for the global denudation rates compilation carried out by Marie-Alice Harel 
# ***** DO NOT LAUNCH THIS WITH SPYDER ! *****
# *** USE A COMMAND LINE IN A TERMINAL INSTEAD ***
# ==================================

"""
Created on Thu Mar 26 17:07:35 2015

@author: mharel
"""

#!/usr/bin/python
import gdal, ogr, osr, numpy
import csv
import pandas
#from shapely.wkb import loads
#from shapely.ops import cascaded_union

#--------------------------------------------------------------------------
# Enter zone number :
zonenb = 67   
newzonenb = 67   # New zone number if zone name belongs to older referencing

# Choose a parameter to be averaged over the zone basins :
paramch = 'seism'
#  vege
#  seism
#--------------------------------------------------------------------------

DataDirectory ="/exports/csce/datastore/geos/users/mharel/Topo_Data/general/New_zone_ref/zone"+str(zonenb)+"/"
print DataDirectory

def proc():
                
        # Each process needs its own pointer.
        shape = ogr.Open(input_zone_polygon)
        lyr = shape.GetLayer()
        raster = gdal.Open(input_value_raster)
        xmin, xmax, ymin, ymax = lyr.GetExtent()
        #print "Extent : ", lyr.GetExtent()
#        # Get extent of feat
#        geom = feat.GetGeometryRef()
#        if (geom.GetGeometryName() == 'MULTIPOLYGON'):
#            count = 0
#            pointsX = []; pointsY = []
#            for polygon in geom:
#                geomInner = geom.GetGeometryRef(count)    
#                ring = geomInner.GetGeometryRef(0)
#                numpoints = ring.GetPointCount()
#                for p in range(numpoints):
#                        lon, lat, z = ring.GetPoint(p)
#                        pointsX.append(lon)
#                        pointsY.append(lat)    
#                count += 1
#        elif (geom.GetGeometryName() == 'POLYGON'):
#            ring = geom.GetGeometryRef(0)
#            numpoints = ring.GetPointCount()
#            pointsX = []; pointsY = []
#            for p in range(numpoints):
#                    lon, lat, z = ring.GetPoint(p)
#                    pointsX.append(lon)
#                    pointsY.append(lat)
#        else:
#            sys.exit()
#    
#        xmin = min(pointsX)
#        xmax = max(pointsX)
#        ymin = min(pointsY)
#        ymax = max(pointsY)
        
#        print geom.GetGeometryName()
        # Specify offset and rows and columns to read
        xoff = int((xmin - xOrigin)/pixelWidth)
        yoff = int((yOrigin - ymax)/pixelWidth)
        xcount = int((xmax - xmin)/pixelWidth)+1
        ycount = int((ymax - ymin)/pixelWidth)+1
                
        # Create memory target raster
        target_ds = gdal.GetDriverByName('MEM').Create('', xcount, ycount, gdal.GDT_Byte)
        target_ds.SetGeoTransform((
            xmin, pixelWidth, 0,
            ymax, 0, pixelHeight,
        ))
        # Create for target raster the same projection as for the value raster
        raster_srs = osr.SpatialReference()
        raster_srs.ImportFromWkt(raster.GetProjectionRef())
        target_ds.SetProjection(raster_srs.ExportToWkt())
    
        # Rasterize zone polygon to raster
        gdal.RasterizeLayer(target_ds, [1], lyr, burn_values=[1])
        # Read raster as arrays
        banddataraster = raster.GetRasterBand(1)
        # error here        
        dataraster = banddataraster.ReadAsArray(xoff, yoff, xcount, ycount).astype(numpy.float)
        bandmask = target_ds.GetRasterBand(1)
        datamask = bandmask.ReadAsArray(0, 0, xcount, ycount).astype(numpy.float)

        # Mask zone of raster
        zoneraster = numpy.ma.masked_array(dataraster,  numpy.logical_not(datamask))

        # Calculate statistics of zonal raster
        value = numpy.mean(zoneraster)

        print value
        return value


resu = []  # Empty list to store the results 
# Extract junction numbers for this zone
#data = pandas.read_csv('/home/mharel/LSDVisu_work/compil-data-MAH.csv')
data = pandas.read_csv('compil-data-MAH.csv')
datazon = data[(data.zone == newzonenb)].copy()
#lenzon = len(datazon['Junction'])

for junct in datazon['Junction']: 
    junction = int(junct)
    print "Junctin is " +str(junction)
    
    # Raster dataset
    if paramch == 'eleva':
        input_value_raster = DataDirectory+"zone"+str(zonenb)+".flt"
    else:
        input_value_raster = DataDirectory+paramch+"_zone"+str(zonenb)+".flt"

    # Vector dataset(zones)
    input_zone_polygon = DataDirectory+"shape_"+str(junction)+".shp"
     
    # Open data
    rast = gdal.Open(input_value_raster)
    print rast.RasterXSize, rast.RasterYSize
    shp = ogr.Open(input_zone_polygon)

    if shp is None:
        raise IOError('Could not open file ' + str(input_zone_polygon))

    # Get raster georeference info
    transform = rast.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = transform[5]
             
    # Start the processes
    layer = shp.GetLayer()
    featList = range(layer.GetFeatureCount())
    #print "Number of features = " +str(featList)
               
    #pool = Pool(processes=24)   
    #value = pool.map(proc,0,8)
    value = proc()
    resu.append([value])
    print "value = " +str(value)
    
with open ('/exports/csce/datastore/geos/users/mharel/Topo_Data/general/'+'zone'+str(newzonenb)+'_resufile_'+paramch+'.csv', 'w') as csvfile:
    g = csv.writer(csvfile, delimiter = ',')
    g.writerows(resu)
    
print "Done."
