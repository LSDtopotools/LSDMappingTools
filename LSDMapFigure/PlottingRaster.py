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
from . import PlottingHelpers as phelp
import matplotlib.pyplot as plt
import matplotlib.cm as _cm
import matplotlib.colors as _mcolors
import matplotlib.axes
import numpy as np
from matplotlib import rcParams


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
        self._RasterAspectRatio = (self._RasterExtents[1]-self._RasterExtents[0])/(self._RasterExtents[3]-self._RasterExtents[2])

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
        self.y_min = self._RasterExtents[2]

    @property
    def xmax(self):
        self.x_max = self._RasterExtents[1]

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

    #def get_apect_ratio(self):
    #    # Get the aspect ratio
    #    self._RasterAspectRatio = (self.xmax-self.xmin)/(self.ymax-self.ymin)
    #    return self._RasterAspectRatio

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

    def axis_styler(self,ax_list,axis_style="Normal"):
        """This sets the line width and fonts on the axes. 
        
        Args:
            ax_list (axes objects): the list of axis objects
            axis_style (string): The syle of the axis. See options below. 
        """

        if axis_style == "Normal":
            lw = 1              # line width
            ftsz = 10           # Size of tick label font
            tpd = 2             # Tick padding
            label_ftsz = 12     # Fontsize of axis label
        elif axis_style == "Thick":
            lw = 2
            ftsz = 10
            tpd = 2
            label_ftsz = 12            
        elif axis_style == "Thin":
            lw = 0.5
            ftsz = 8
            tpd = 1
            label_ftsz = 10   
        elif axis_style == "Big":
            lw = 2
            ftsz = 12
            tpd = 3 
            label_ftsz = 14 
        elif axis_style == "Madhouse":
            # This is just a crazy style to test if the figure is actually recieving these instructions
            lw = 4
            ftsz = 20
            tpd = 3 
            label_ftsz = 6             
        else:
            print("Using the default axis styling")
            lw = 1
            ftsz = 10
            tpd = 2
            label_ftsz = 12             
           
        for ax in ax_list:
            # Now to fix up the axes
            ax.spines['top'].set_linewidth(lw)
            ax.spines['left'].set_linewidth(lw)
            ax.spines['right'].set_linewidth(lw)
            ax.spines['bottom'].set_linewidth(lw)
            
            ax.xaxis.label.set_size(label_ftsz)
            ax.yaxis.label.set_size(label_ftsz)
                        
            # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
            ax.tick_params(axis='both', width=lw, pad = tpd, labelsize = ftsz )
            
        return ax_list


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
                        colorbarlabel = "Colourbar",
                        colourbar_orientation = "horizontal"):

        print("N axes are: "+str(len(self.ax_list)))
        print(self.ax_list[0])

        self.ax_list = self._add_drape_image(self.ax_list,RasterName,Directory,colourmap,alpha,show_colourbar,colorbarlabel,colourbar_orientation)

    def _add_drape_image(self,ax_list,RasterName,Directory,
                         colourmap = "gray",
                         alpha=0.5,
                         show_colourbar = False,
                         colorbarlabel = "Colourbar",
                         colourbar_orientation = "horizontal"):

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
                                              colorbarlabel = colorbarlabel, cbar_orientation = colourbar_orientation)


        return self.ax_list

    def add_colourbar(self,ax_list,im,BaseRaster,colorbarlabel = "Colourbar",cbar_orientation = 'horizontal'):
        fig = matplotlib.pyplot.gcf()
        ax_list.append(fig.add_axes([0.1,0.8,0.2,0.5]))
        cbar = plt.colorbar(im,cmap=BaseRaster._colourmap,spacing='uniform', orientation=cbar_orientation,cax=ax_list[-1])
        #cbar.set_label(colorbarlabel, fontsize=10)
        
        if cbar_orientation == 'horizontal': 
            ax_list[-1].set_xlabel(colorbarlabel, fontname='Arial',labelpad=-35)
        else:
            ax_list[-1].set_ylabel(colorbarlabel, fontname='Arial',labelpad=-35)

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
 
    def save_fig(self,fig_width_inches = 4,cbar_location = "Top",FigFileName = 'TestFig.png',FigFormat = 'png',Fig_dpi = 100, axis_style = "Normal"):

        
        self.ax_list = self.axis_styler(self.ax_list,axis_style)
        
        map_aspect_ratio = self._RasterList[0]._RasterAspectRatio
        print("The aspect ratio is: "+str(map_aspect_ratio))

        fig = matplotlib.pyplot.gcf()
        fig_size_inches, map_axes, cbar_axes = phelp.MapFigureSizer(fig_width_inches,
                                                              map_aspect_ratio,cbar_loc = cbar_location)

        fig.set_size_inches(fig_size_inches[0], fig_size_inches[1])
        self.ax_list[0].set_position(map_axes)

        print("Number of axes are: " + str(len(self.ax_list)))


        if cbar_axes == None:
            del self.ax_list[-1]
        else:
            self.ax_list[-1].set_position(cbar_axes)

        fig.savefig(FigFileName, format=FigFormat, dpi=Fig_dpi)

        #self.fig.show()
        #print("The figure format is: " + self.FigFormat)
        #plt.savefig(self.FigFileName,format=self.FigFormat)
        #self.fig.clf()

    def SetRCParams(self,label_size):
        """This sets a load of RC params
        """
        print("I am setting the font size to: "+str(label_size))
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['arial']
        rcParams['font.size'] = label_size
        rcParams['lines.linewidth']  = 1.5
