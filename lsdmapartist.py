# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 21:50:53 2017

LSDMapArtist

@author: dav

Object-oriented plotting module for constructing
drape maps in a reusable, generic way.

Experimental. Use at your own risk.

This software is realsed under the Artistic Licence v2.0

"""

import LSDPlottingTools as LSDP
import matplotlib.pyplot as plt


class BaseRaster(object):
    """ 
    Class BaseRaster represents the data associated with the basic rasters
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
    drape data can be overlain. (overlaid? overlayed? superimposed?)
    """
    def __init__(self, RasterName, Directory, backgroundtype="Hillshade"):
        
        # We need to initialise our base class though first!
        super(BackgroundRaster, self).__init__(RasterName, Directory)
        
        # Could also have done it this way:
        # BaseRaster.__init__(RasterName, Directory)
        
        self._backgroundtype = backgroundtype
        
        self._render_background(self.fullpath_to_raster)
        #self._hillshade_raster = LSDP.Hillshade(self.fullpath_to_raster)
        
    def _render_background(self, fullpath_to_raster):
        """
        Renders the background image that 
        will form the drape plot, e.g. a hillshade
        """
        if self._backgroundtype == "Hillshade":
          self.Hillshade = LSDP.Hillshade(self.fullpath_to_raster)
        else:
          print ("That background style is not yet supported. Currently only " \
                " 'Hillshade' is supported")
        
    def show_hillshade(self):
        """
        For testing the hillshade
        """
        plt.imshow(self.Hillshade,
                   cmap="gray",
                   extent=self.extents)
        
class DrapeRaster(BaseRaster):
    """
    This class is for constructing the drape raster which gets overlaid on
    top of the BackgroundRaster
    """
    def __init__(self, RasterName, Directory):
        super(DrapeRaster, self).__init__(RasterName, Directory)
        
        self._render_drape(self.fullpath_to_raster)      
        
    def _render_drape(self, fullpath_to_raster):
        pass

class DrapePlot(object):
	"""
	Class DrapeMap contains the methods for creating a drape map, i.e
	a BackgroundRaster overlain with a DrapeRaster.
	"""
    def __init__(self, DrapeRasterName, BackgroundRasterName, Directory):
		
		self._backgroundRaster = BackgroundRaster(BackgroundRasterName, Directory)
		self._drapeRaster = DrapeRaster(DrapeRasterName, Directory)
		
	def make_drape_plot():
		"""Creates a matplotlib Axes object with the drape map.
		
		Returns:
		    A matplotlib.Axes instance containing the drape plot.
		"""
		return drape_plot_axes

class MultiDrapePlot(object):
	"""
	Class MultiDrapeMap creates multiplot drape maps.
	
	Composed of multiple DrapePlot
	"""
    pass


DataDirectory = "/mnt/SCRATCH/Dev/LSDMappingTools/test_rasters/peak_flow_rasters/"
RasterName = "Elevations0.asc"

raster = BaseRaster(RasterName, DataDirectory)

