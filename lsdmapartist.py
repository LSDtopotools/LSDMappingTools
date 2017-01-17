#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 21:50:53 2017

LSDMapArtist

@author: dav

Object-oriented plotting module for constructing
drape maps in a reusable, generic way.

(This is also partly an exercise in writing more complex classes in python...)

"""

import LSDPlottingTools as LSDP

class BaseRaster(object):
    """ 
    Class BaseRaster represents the data associated with the basic raster
    used to create to image to be plotted. It also contains the methods
    associated with performing any basic analyses to that raster.
    """
    def __init__(self, RasterName, Directory):
        
        self._RasterFileName = RasterName
        self._RasterDirectory = Directory
        self._FullPathRaster = self._RasterDirectory + self._RasterFileName
        
        # I think the BaseRaster should contain a numpy array of the Raster
        self._RasterArray = LSDP.ReadRasterArrayBlocks(self._FullPathRaster)
        
        # Get the extents as a list
        self._RasterExtents = LSDP.GetRasterExtent(self._FullPathRaster)
    
    @property
    def extents(self):
        return self._RasterExtents
        
    @property
    def fullpath_to_raster(self):
        return self._FullPathRaster


class BackgroundRaster(BaseRaster):
    """
    Class Background raster represents the background image over which the
    drape data can be overlain.
    """
    def __init__(self, RasterName, Directory):
        #BaseRaster.__init__(RasterName, Directory)
        super(BackgroundRaster, self).__init__(RasterName, Directory)
        
        self._hillshade_raster = LSDP.Hillshade(self.fullpath_to_raster)

class DrapeRaster(BaseRaster):
    pass        

class DrapeMap:
    pass

class MultiDrapeMap:
    pass


DataDirectory = "/mnt/SCRATCH/Dev/LSDMappingTools/test_rasters/peak_flow_rasters/"
RasterName = "Elevations0.asc"

raster = BaseRaster(RasterName, DataDirectory)

