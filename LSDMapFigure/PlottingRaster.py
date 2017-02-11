# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 21:50:53 2017

LSDPlottingRaster

@author: DAV and SMM

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
        
        # set the default colourmap
        self._colourmap = "gray"
    
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

    # The default colormap is gray
    @property
    def colourmap(self):
        self._colourmap = "gray"
    
    # This controls the zorder of the raster
    @property
    def zorder(self):
        self._zorder = 1
        
        
    def set_raster_type(self, rastertype):
        """
        Renders the background image that 
        will form the drape plot, e.g. a hillshade
        """
        self._rastertype = rastertype
        if self._rastertype == "Hillshade":
            print("I'm using a hillshade colour scheme")
            self._colourmap = "gray"
          
        elif self._rastertype == "Terrain":
            print("I'm using a terrain colour scheme")
            self._colourmap = LSDP.colours.UsefulColourmaps.darkearth
        else:
            print ("That raster style is not yet supported. Currently only " \
                   " 'Hillshade' and 'Terrain' are supported.")

    def set_colourmap(self, cmap):
        """
        There must be a more pythonic way to do this!
        """
        self._colourmap = cmap
            
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

    def show_raster(self):
        plt.imshow(self._RasterArray,
                   cmap=self._colourmap,
                   extent=self.extents)
        plt.show()        
        
class MapFigure(object):        
    def __init__(self, BaseRasterName, Directory, 
                 coord_type="UTM", *args, **kwargs):	        

        # A map figure has one figure
        #self.fig = plt.figure(1, facecolor='white',figsize=(6,3))
        #self.fig
        
        # There can be mulitple axes in the figure. These are maptlotlib artists
        # that can be used to place plotting elements
        #self.ax = self.fig.add_axes([0.1,0.1,0.7,0.7])
        #self.ax = plt.subplots()
        
        #self.fig, self.ax = plt.subplots()
        fig = plt.figure()
        self.ax_list = []
        ax = fig.add_axes([0,0,1,1])
        self.ax_list.append(ax)
        

        # Names of the directory and the base raster        
        self._Directory = Directory
        self._BaseRasterName = BaseRasterName
        self._BaseRasterFullName = Directory+BaseRasterName
 
        
        self.FigFileName = self._Directory+"TestFig.png"
        self.FigFormat = "png"        
        
        
        # The way this is going to work is that you can have many rasters in the
        # plot that get appended into a list. Each one has its own colourmap
        # and properties
        self._RasterList = []
        self._RasterList.append(BaseRaster(BaseRasterName,Directory))
        
        # The coordinate type. UTM and UTM with tick in km are supported at the moment         
        self._set_coord_type(coord_type)
        
        # Get the tick properties 
        self._xmin = self._RasterList[0].xmin
        self._ymin = self._RasterList[0].ymin
        self._xmax = self._RasterList[0].xmax
        self._ymax = self._RasterList[0].ymax
        self._n_target_ticks = 5 
        self.make_ticks()
        
        self._num_drapes = 0  # Number of drapes in the image.
        # We will increase this everytime we call ax.imshow.
        
        # Stores the Image instances generated from imshow()
        self._drape_list = []
        
        self.ax_list = self.make_base_image(self.ax_list)
        print(self.ax_list[0])


    def make_ticks(self):
        if self._coord_type == "UTM":
            self.tick_xlocs,self.tick_ylocs,self.tick_x_labels,self.tick_y_labels = LSDP.GetTicksForUTMNoInversion(self._BaseRasterFullName,self._xmax,self._xmin,
                             self._ymax,self._ymin,self._n_target_ticks)
        elif self._coord_type == "UTM_km":
            self.tick_xlocs,self.tick_ylocs,self.tick_x_labels,self.tick_y_labels = LSDP.GetTicksForUTMNoInversion(self._BaseRasterFullName,self._xmax,self._xmin,
                             self._ymax,self._ymin,self._n_target_ticks)
            n_hacked_digits = 3
            self.tick_x_labels = LSDP.TickLabelShortenizer(self.tick_x_labels,n_hacked_digits)
            self.tick_y_labels = LSDP.TickLabelShortenizer(self.tick_y_labels,n_hacked_digits)
        else:
            raise ValueError("Sorry, the coordinate type: ", self._coord_type, 
                             "is not yet supported") 
            
        print("I made the ticks.")
        print("x labels are: ")
        print(self.tick_x_labels)
        print("x locations are:")
        print(self.tick_xlocs)
        print("y labels are: ")
        print(self.tick_y_labels)
        print("y locations are:")
        print(self.tick_ylocs)

    def add_ticks_to_axis(self,ax):
        ax.set_xticklabels(self.tick_x_labels)
        ax.set_yticklabels(self.tick_y_labels)
        ax.set_xticks(self.tick_xlocs)
        ax.set_yticks(self.tick_ylocs)
        ax.set_xlabel(self._xaxis_label)
        ax.set_ylabel(self._yaxis_label)
        
        return ax
        
    def make_base_image(self,ax_list):
        
        # We need to initiate with a figure
        #self.ax = self.fig.add_axes([0.1,0.1,0.7,0.7])

        print("This colourmap is: "+ self._RasterList[0]._colourmap)
        im = self.ax_list[0].imshow(self._RasterList[0]._RasterArray, self._RasterList[0]._colourmap, extent = self._RasterList[0].extents, interpolation="nearest")
        
        # This affects all axes because we set share_all = True.
        #ax.set_xlim(self._xmin,self._xmax)
        #ax.set_ylim(self._ymin,self._ymax)
        self.ax_list[0] = self.add_ticks_to_axis(self.ax_list[0])       
        self._drape_list.append(im)
        
        print("The number of axes are: "+str(len(self._drape_list)))
        
        print(self.ax_list[0])
        return self.ax_list
        
    def add_drape_image(self,RasterName,Directory,colourmap = "gray",
                        alpha=0.5,
                        show_colourbar = False,
                        colorbarlabel = "Colourbar"):
        
        print("N axes are: "+str(len(self.ax_list)))
        print(self.ax_list[0])
        
        self.ax_list = self._add_drape_image(self.ax_list,RasterName,Directory,colourmap,alpha,show_colourbar,colorbarlabel)

    def _add_drape_image(self,ax_list,RasterName,Directory,
                         colourmap = "gray",
                         alpha=0.5,
                         show_colourbar = False,
                         colorbarlabel = "Colourbar"):
         
        self._RasterList.append(BaseRaster(RasterName,Directory))
        self._RasterList[-1].set_colourmap(colourmap)
        
        # We need to initiate with a figure
        #self.ax = self.fig.add_axes([0.1,0.1,0.7,0.7])

        im = self.ax_list[0].imshow(self._RasterList[-1]._RasterArray, self._RasterList[-1]._colourmap, extent = self._RasterList[0].extents, interpolation="nearest",alpha = alpha)
        
        # This affects all axes because we set share_all = True.
        #ax.set_xlim(self._xmin,self._xmax)
        #ax.set_ylim(self._ymin,self._ymax)
        self.ax_list[0] = self.add_ticks_to_axis(self.ax_list[0])       
        self._drape_list.append(im)
        
        print("The number of axes are: "+str(len(self._drape_list)))
        
        if show_colourbar:
            self.ax_list = self.add_colourbar(self.ax_list,im,self._RasterList[-1],
                                              colorbarlabel = colorbarlabel)
            
        
        return self.ax_list
        
    def add_colourbar(self,ax_list,im,BaseRaster,colorbarlabel = "Colourbar"):
        fig = matplotlib.pyplot.gcf()
        ax_list.append(fig.add_axes([0.1,0.8,0.2,0.5]))
        cbar = plt.colorbar(im,cmap=BaseRaster._colourmap,spacing='uniform', orientation='horizontal',cax=ax_list[-1])
        #cbar.set_label(colorbarlabel, fontsize=10)
        ax_list[-1].set_xlabel(colorbarlabel, fontname='Arial',labelpad=-35) 
        
        return ax_list

        
        
        
    def _set_coord_type(self, coord_type):
        """Sets the coordinate type"""
        if coord_type == "UTM":
            self._coord_type = "UTM"
            self._xaxis_label = "Easting (m)"
            self._yaxis_label = "Northing (m)"

        elif coord_type == "UTM_km":
            self._coord_type = "UTM_km"
            self._xaxis_label = "Easting (km)"
            self._yaxis_label = "Northing (km)"            
            
        # Example, do not actually use...
        elif coord_type == "Kruskal–Szekeres":
            self._coord_type = "Kruskal–Szekeres"
            self._xaxis_label = "X"
            self._yaxis_label = "T"
        
        else:
            raise ValueError("Sorry, the coordinate type: ", coord_type, 
                             "is not yet supported")        

    def show_plot(self):
        
        self.fig.show()            
 
    def save_fig(self):
        
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(4, 3)
        fig.savefig('test2png.png', dpi=100)        
        
        #self.fig.show()   
        #print("The figure format is: " + self.FigFormat)
        #plt.savefig(self.FigFileName,format=self.FigFormat)
        #self.fig.clf()           