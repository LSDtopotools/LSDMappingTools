#!/usr/bin/python
#coding=utf-8

"""
Description : Crop and project an image in the same grid as a reference image

Author : Amaury Dehecq
"""

#Python libraries
from osgeo import gdal, osr
import sys, os

# =============================================================================
def Usage():
    print('*** Crop and project an image in the same grid as a reference image ***')
    print()
    print('Usage : crop2image.py im_ref im2crop outfile')
    print()
    print('input parameters :')
    print('  im_ref : str, path to the reference image')
    print('  im2crop : str, path to the image to crop')
    print('  outfile : str, name of the output file')
    print() 
# =============================================================================
   
# example of command line :
#-------------------------
# python crop2image.py /exports/csce/datastore/geos/users/mharel/Topo_Data/general/Old_zone_ref/zone69/zone69.flt /exports/csce/datastore/geos/users/mharel/World_data/seism_globe.flt /exports/csce/datastore/geos/users/mharel/Topo_Data/general/Old_zone_ref/zone69/seism_zone69.flt

if len(sys.argv)<4:
    Usage()
    sys.exit(1)


#Read arguments
im_ref = sys.argv[1]
im2crop = sys.argv[2]
outfile = sys.argv[3]


#Read reference image spatial reference system
ds = gdal.Open(im_ref)
srs_dest=osr.SpatialReference()
wkt = ds.GetProjectionRef()
srs_dest.ImportFromWkt(wkt)
proj=srs_dest.ExportToProj4()


#Read reference image size and corners coordinates
trans = ds.GetGeoTransform()
pixelWidth = trans[1]
pixelHeight = trans[5]
xmin = trans[0]
xmax = trans[0] + ds.RasterXSize*trans[1] + ds.RasterYSize*trans[2]
ymin = trans[3] + ds.RasterXSize*trans[4] + ds.RasterYSize*trans[5]
ymax = trans[3]

#Crop and reproject
cmd = "gdalwarp -te %.8f %.8f %.8f %.8f -tr %.8f %.8f -t_srs '%s' %s %s -overwrite" %(xmin,ymin,xmax,ymax,pixelWidth,pixelHeight,proj,im2crop,outfile)
print(cmd); os.system(cmd)
