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
import matplotlib.cm as _cm
import numpy as np


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
    
    @property
    def xmin(self):
        self.x_min = self._RasterExtents[0]

    @property
    def ymin(self):
        self.y_min = self._RasterExtents[1]
        
    @property
    def xmax(self):
        self.x_max = self._RasterExtents[2]
        
    @property
    def ymax(self):
        self.x_max = self._RasterExtents[3]

        
class BackgroundRaster(BaseRaster):
    """
    Class Background raster represents the background image over which the
    drape data can be overlain. (overlaid? overlayed? superimposed?)
    """
    def __init__(self, RasterName, Directory, background_type="Hillshade"):
        
        # We need to initialise our base class though first!
        super(BackgroundRaster, self).__init__(RasterName, Directory)
        
        # Could also have done it this way:
        # BaseRaster.__init__(RasterName, Directory)
        
        self._backgroundtype = background_type
        
        self._render_background(self.fullpath_to_raster)
        #self._hillshade_raster = LSDP.Hillshade(self.fullpath_to_raster)
        
    def _render_background(self, fullpath_to_raster):
        """
        Renders the background image that 
        will form the drape plot, e.g. a hillshade
        """
        if self._backgroundtype == "Hillshade":
            self.Hillshade = LSDP.Hillshade(self.fullpath_to_raster)
            self.colourmap = "gray"
          
        elif self._backgroundtype == "Terrain":
            self.Hillshade = LSDP.ReadRasterArrayBlocks(self.fullpath_to_raster)
            self.colourmap = "terrain"
        else:
            print ("That background style is not yet supported. Currently only " \
                   " 'Hillshade' and 'Terrain' are supported.")
        
    def show_hillshade(self):
        """
        For testing the hillshade
        """
        plt.imshow(self.Hillshade,
                   cmap=self.colourmap,
                   extent=self.extents)
        
class DrapeRaster(BaseRaster):
    """
    This class is for constructing the drape raster which gets overlaid on
    top of the BackgroundRaster
    """
    def __init__(self, RasterName, Directory,
                 drape_min_threshold, drape_max_threshold,
                 middle_mask_range):
        
        # Initialise the super-class
        super(DrapeRaster, self).__init__(RasterName, Directory)
        
        self._render_drape(self.fullpath_to_raster) 

        self._drapeminthreshold = drape_min_threshold
        self._drapemaxthreshold = drape_max_threshold
        self._middlemaskrange = middle_mask_range
        
        self._initialise_masks()
        
    def _initialise_masks(self):
        if self._drapeminthreshold is not None:
            self.mask_low_values()
        if self._drapemaxthreshold is not None:
            self.mask_high_values()
        if self._middlemaskrange is not None:
            self.mask_middle_values()
        
    def mask_low_values(self):
        low_values_index = self._RasterArray < self._drapeminthreshold
        self._RasterArray[low_values_index] = np.nan
                         
    def mask_high_values(self):
        high_values_index = self._RasterArray < self._drapemaxthreshold
        self._RasterArray[high_values_index] = np.nan
                         
    def mask_middle_values(self):
        """Masks a centre range of values."""
        masked_mid_values_index = (np.logical_and(self._RasterArray > self._middlemaskrange[0], 
                                   self._RasterArray < self._middlemaskrange[1]))
        self._RasterArray[masked_mid_values_index] = np.nan
        
    def _render_drape(self, fullpath_to_raster):
        pass
    
    def set_drape_min_threshold(self):
        pass
    
    def set_drape_max_threshold(self):
        pass

class DrapePlot(object):
    """
    Class DrapeMap contains the methods for creating a drape map, i.e
    a BackgroundRaster overlain with a DrapeRaster.
    """
    def __init__(self, DrapeRasterName, BackgroundRasterName, Directory,
                 drape_colourmap, background_type="Hillshade",
                 drape_min_threshold=None, drape_max_threshold=None,
                 vmin=None, vmax = None, middle_mask_range=None):	
        
        self.Background = BackgroundRaster(BackgroundRasterName, Directory,
                                           background_type)
        self.Drape = DrapeRaster(DrapeRasterName, Directory,
                                 drape_min_threshold=drape_min_threshold,
                                 drape_max_threshold=drape_max_threshold,
                                 middle_mask_range=middle_mask_range)
        self._drape_colourmap = drape_colourmap
        
        self._vmin = vmin
        self._vmax = vmax

        
        self.make_drape_plot()
		
    def make_drape_plot(self):
        """Creates a matplotlib Axes object with the drape map."""
        self.fig, self.ax = plt.subplots()
        
        # Plot the background
        self.im = self.ax.imshow(self.Background.Hillshade,
                                 self.Background.colourmap,
                                 extent=self.Background.extents,
                                 interpolation="nearest")
        # Plot the drape (overlay data) on top.
        self.im = self.ax.imshow(self.Drape._RasterArray,
                                 self._drape_colourmap,
                                 extent=self.Drape.extents,
                                 interpolation="nearest",
                                 vmin=self._vmin, vmax=self._vmax)

    def make_drape_colourbar(self, cbar_label="Test"):
        """Adds a colourbar for the drape"""
        self.fig.subplots_adjust(right=0.85)
        self.cax = self.fig.add_axes([0.9, 0.1, 0.03, 0.8])
        
        self.cbar = self.fig.colorbar(self.im, cax=self.cax) 
        self.cbar.set_label(cbar_label)
        # This is needed to update the canvas figure with the colour bar.
        self.fig.canvas.draw()
        
    def set_coordinate_labels_axis(self, x_axis_label='Easting (m)',
                                   y_axis_label='Northing (m)'):
        """Sets the x and y axis labels for the entire figure"""
        self.fig.text(0.5, 0.04, x_axis_label, ha='center')
        self.fig.text(0.04, 0.5, y_axis_label, va='center', rotation='vertical')
        self.fig.canvas.draw()
        

    def show_plot(self):
        self.fig.show()
        

class MultiDrapePlot(object):
    """
    Class MultiDrapeMap creates multiplot drape maps.
    
    Composed of multiple DrapePlot
    """
    pass




