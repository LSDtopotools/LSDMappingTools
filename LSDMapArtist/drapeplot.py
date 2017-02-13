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

# LSDPlottingTools must be in your pythonpath
import LSDPlottingTools as LSDP
import matplotlib.pyplot as plt
import matplotlib.cm as _cm
import matplotlib.colors as _mcolors
import matplotlib.axes
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
    def raster_filename(self):
        return self._RasterFileName
    
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
            #self.colourmap = LSDP.colours.UsefulColourmaps.niceterrain
            self.colourmap = LSDP.colours.UsefulColourmaps.darkearth
            #self.colourmap = "terrain"
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

class DrapeAxes(object):
    """
    Class DrapeAxes contains the methods for creating a drape map, i.e
    a BackgroundRaster overlain with a DrapeRaster.
    """
    def __init__(self, DrapeRasterName, BackgroundRasterName, Directory,
                 drape_colourmap="jet", 
                 background_type="Hillshade", 
                 show_background_colourbar=False,
                 colourbar_label="",
                 drape_min_threshold=None, 
                 drape_max_threshold=None,
                 colourbar_norm_type=None,
                 vmin=None, 
                 vmax=None, 
                 middle_mask_range=None,
                 drape_alpha=1.0,
                 clip_to_background=False,
                 coord_type="UTM", *args, **kwargs):	
        
        
#        for key, value in kwargs.items():
#            setattr(self, '_'+key, value)
        self.fig, self.ax = plt.subplots()
        
        self.Background = BackgroundRaster(BackgroundRasterName, 
                                           Directory,
                                           background_type)

        self.Drape = DrapeRaster(DrapeRasterName, 
                                 Directory,
                                 drape_min_threshold=drape_min_threshold,
                                 drape_max_threshold=drape_max_threshold,
                                 middle_mask_range=middle_mask_range)

        self._drape_colourmap = drape_colourmap
        self._colourbar_label = colourbar_label
        self._show_background_colourbar = show_background_colourbar
        self._vmin = vmin
        self._vmax = vmax
        
        self._colourbar_normalisation = self.set_colourbar_norm_type(colourbar_norm_type)
        
        self._drape_alpha = drape_alpha
        
        self._set_coord_type(coord_type)
        
        self._num_drapes = 0  # Number of drapes in the image.
        # We will increase this everytime we call ax.imshow.
        
        # Stores the Image instances generated from imshow()
        self._drape_list = []
        
        if clip_to_background:
            self.mask_by_background_nodatas()

        self.make_drape_plot()
    
    #def __call__(self):
    #    return self.ax
	
    def set_colourbar_norm_type(self, colourbar_norm_type):
        if colourbar_norm_type == "SymLogNorm":
            norm = _mcolors.SymLogNorm(linthresh=0.3, vmin=self._vmin, vmax=self._vmax)
        elif colourbar_norm_type is None:
            norm = None
        else:
            raise NotImplementedError("That type of normalisation is not coded in yet..")
            norm = None
        
        return norm
        
	
    def make_drape_plot(self):
        """Creates a matplotlib Axes object with the drape map."""
        
        # but what if fig XOR ax is None??
#        if fig is None and ax is None:
#            self.fig, self.ax = plt.subplots()
#        else:
#            self.fig=fig, self.ax=ax

        
        
        # Plot the background
        self.im_background = self.ax.imshow(self.Background.Hillshade,
                                 self.Background.colourmap,
                                 extent=self.Background.extents,
                                 interpolation="nearest")
        self._num_drapes += 1
        self._drape_list.append(self.im_background)
        
        if self._show_background_colourbar:
            # Plot the background image colour bar
            self._generic_colourbar_plotter(self.im_background, "Elevation (m)")
        
        # Plot the drape (overlay data) on top.
        # Should be separate function really...
        self.im = self.ax.imshow(self.Drape._RasterArray,
                                 self._drape_colourmap,
                                 extent=self.Drape.extents,
                                 interpolation="nearest",
                                 vmin=self._vmin, vmax=self._vmax,
                                 norm=self._colourbar_normalisation,
                                 alpha=self._drape_alpha
                                 )
                                 #norm=_mcolors.PowerNorm(gamma=0.2))
        
        self._drape_list.append(self.im)
        self._num_drapes += 1
        
        self._set_axis_labels(self._xaxis_label, self._yaxis_label)
        
        # Add the colourbar for the drape
        self._generic_colourbar_plotter(self.im, self._colourbar_label)
        
        # Add a title
        self._set_subplot_autolabel()
        
        #return self.fig, self.ax
    
    def _set_coord_type(self, coord_type):
        """Sets the coordinate type"""
        if coord_type == "UTM":
            self._coord_type = "UTM"
            self._xaxis_label = "Easting (m)"
            self._yaxis_label = "Northing (m)"
        
        # Example, do not actually use...
        elif coord_type == "Kruskal–Szekeres":
            self._coord_type = "Kruskal–Szekeres"
            self._xaxis_label = "X"
            self._yaxis_label = "T"
        
        else:
            raise NotImplementedError("Sorry, the coordinate type: ", coord_type, 
                             "is not yet supported")

    def _generic_colourbar_plotter(self, mappable, cbar_label):
        """A generic colourbar plotter"""
        self.cbar = self.fig.colorbar(mappable, fraction=0.046, pad=0.04)
        self.cbar.set_label(cbar_label)
        self.fig.canvas.draw()
        
    def make_drape_colourbar(self, cbar_label="Test"):
        """Adds a colourbar for the drape"""
        #self.fig.subplots_adjust(right=0.85)
        #self.cax = self.fig.add_axes([0.9, 0.1, 0.03, 0.8])
        
        self.cbar = self.fig.colorbar(self.im)#, cax=self.cax) 
        self.cbar.set_label(cbar_label)
        # This is needed to update the canvas figure with the colour bar.
        self.fig.canvas.draw()
      
    def _set_axis_labels(self, x_axis_label="", y_axis_label=""):
        self.ax.set_xlabel(x_axis_label)
        self.ax.set_ylabel(y_axis_label)
        
    def set_fig_axis_labels(self, x_axis_label='Easting (m)',
                                   y_axis_label='Northing (m)'):
        """Sets the x and y axis labels for the entire figure
        
        This is a hacky-method since figures don't technically,
        have axis labels, so you are actually setting a text box
        value and manually positioning it on the figure
        """
        self.fig.text(0.5, 0.04, x_axis_label, ha='center')
        self.fig.text(0.04, 0.5, y_axis_label, va='center', rotation='vertical')
        self.fig.canvas.draw()

    def _set_subplot_autolabel(self):
        subplotlabel = LSDP.labels.make_line_label(self.Drape.raster_filename)
        self.ax.set_title(subplotlabel)
        
    def set_subplot_labels(self, text):
        self.ax.set_title(text)
        
    def mask_by_background_nodatas(self):
        """Uses the Background raster mask to mask the drape extents"""
        masked_background_index = np.isnan(self.Background._RasterArray) #== np.nan or -9999.0
        self.Drape._RasterArray[masked_background_index] = np.nan

    def show_plot(self):
        self.fig.show()
    
    @property    
    def num_drapes(self):
        return self._num_drapes
    
    @property
    def coord_type(self):
        return self._coord_type
        

class DrapeFig(object):
    """
    Class MultiDrapeMap creates multiplot drape maps.
    
    Composed of multiple DrapeAxes
    """
    pass




