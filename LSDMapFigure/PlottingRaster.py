# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 21:50:53 2017

LSDPlottingRaster

@author: DAV, SMM

Object-oriented plotting module for constructing
drape maps in a reusable, generic way.

Experimental. Use at your own risk.

This software is released under the Artistic Licence v2.0

"""

# LSDPlottingTools must be in your pythonpath
import LSDPlottingTools as LSDP
from . import PlottingHelpers as phelp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as _cm
import matplotlib.colors as _mcolors
from matplotlib import colors
import matplotlib.axes
import numpy as np
from matplotlib import ticker
from matplotlib import rcParams


class BaseRaster(object):
    """
    Class BaseRaster represents the data associated with the basic rasters
    used to create to image to be plotted. It also contains the methods
    associated with performing any basic analyses to that raster.

    Args:
        RasterName (str): The name of the rasters (with extension). It is read by gdal so should cope with mulitple formats
        Directory (str): The path to the raster. Needs to have the trailing slash
        NFF_opti (bool): experimental test of reading raster using numpy.fromfile() which a super efficient binary reader

    Author: DAV and SMM
    """
    def __init__(self, RasterName, Directory, NFF_opti = False, alpha = 1):

        self._RasterFileName = RasterName
        self._RasterDirectory = Directory
        self._FullPathRaster = self._RasterDirectory + self._RasterFileName

        # I think the BaseRaster should contain a numpy array of the Raster
        if(NFF_opti):
            self._RasterArray = LSDP.ReadRasterArrayBlocks_numpy(self._FullPathRaster)
        else:
            self._RasterArray = LSDP.ReadRasterArrayBlocks(self._FullPathRaster)

        # Get the extents as a list
        self._RasterExtents = LSDP.GetRasterExtent(self._FullPathRaster)
        self._RasterAspectRatio = (self._RasterExtents[1]-self._RasterExtents[0])/(self._RasterExtents[3]-self._RasterExtents[2])

        # set the default colourmap
        self._colourmap = "gray"

        #alpha_value
        self._alpha = alpha

        # get the EPSG string
        self._EPSGString = LSDP.LSDMap_IO.GetUTMEPSG(self._FullPathRaster)
        #print("The EPSGString is: "+ self._EPSGString)

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

        Args:
            rastertype (str): The type of raster. Not many supported, but basically just changes the colourmap

        Author: DAV
        """
        self._rastertype = rastertype
        if self._rastertype == "Hillshade":
            print("I'm using a hillshade colour scheme")
            self._colourmap = "gray"

        elif self._rastertype == "Terrain":
            print("I'm using a terrain colour scheme")
            self._colourmap = LSDP.colours.UsefulColourmaps.darkearth
        else:
            print("That raster style is not yet supported. Currently only " \
                   " 'Hillshade' and 'Terrain' are supported.")

    def set_colourmap(self, cmap):
        """
        There must be a more pythonic way to do this! But sets the colourmap

        Args:
            cmap (list or str): the colourmap

        Author: DAV
        """
        self._colourmap = cmap

    def _initialise_masks(self):
        """
        Some logic Declan wrote to mask values.

        Author: DAV
        """
        if self._drapeminthreshold is not None:
            self.mask_low_values()
        if self._drapemaxthreshold is not None:
            self.mask_high_values()
        if self._middlemaskrange is not None:
            self.mask_middle_values()

    def mask_low_values(self):#
        """
        Reads from the self._drapeminthreshold to mask low values.

        Author: DAV
        """
        low_values_index = self._RasterArray < self._drapeminthreshold
        self._RasterArray[low_values_index] = np.nan

    def mask_high_values(self):
        """
        Reads from the self._drapeminthreshold to mask high values.

        Author: DAV
        """
        high_values_index = self._RasterArray < self._drapemaxthreshold
        self._RasterArray[high_values_index] = np.nan

    def mask_middle_values(self):
        """
        Masks a centre range of values.

        Author: DAV
        """
        masked_mid_values_index = (np.logical_and(self._RasterArray > self._middlemaskrange[0],
                                   self._RasterArray < self._middlemaskrange[1]))
        self._RasterArray[masked_mid_values_index] = np.nan

    def show_raster(self):
        """
        Low level show function. Only really used for debugging since it contains
        no sizing, labelling etc.

        Author: DAV
        """
        plt.imshow(self._RasterArray,
                   cmap=self._colourmap,
                   extent=self.extents)
        plt.show()

    def replace_raster_values(self, old_values = [], new_values = []):
        """
        Function to take a list of raster values and replace it with
        a new list. Can be used to overwrite basin junction IDs with other
        information about the basin, for example.

        Args:
            old_values (list): The old values in the raster
            new_values (list): The replacement values. Needs to be the same size
                as the old_values list

        Author: FJC

        Date: 17/06/17
        """
        for idx, value in enumerate(old_values):
            old_values_index = self._RasterArray == value
            self._RasterArray[old_values_index] = new_values[idx]
        #print self._RasterArray


class MapFigure(object):
    """
    This is the main object used for plotting. It contains the underlying axes of the figures.
    At the moment the4 axes contain typically a colourbar and a main image axis.
    The function also contains routines for sizing, draping, adding point data,
    etc.
    """
    def __init__(self, BaseRasterName, Directory,
                 coord_type="UTM", colourbar_location = "None", basemap_colourmap = "gray", plot_title = "None", NFF_opti = False,alpha = 1,*args, **kwargs):
        """
        Initiates the object.

        Args:
            BaseRasterName (string): The name of the raster (no directory, but need extension)
            Directory (string): directory of the data.
            coord_type (str): The type of coordinate system. At the moment we only support UTM and UTM_km (the latter makes the tick labels easier to handle).
            colourbar_location (string): Can be none, top bottom left or right. Controls where the colourbar is located.
            basemap_colourmap (string or colormap): The colourmap of the base raster.
            plot_title (string): The title of the plot, if "None" then will not be plotted.
            NFF_opti (bool): If true, use a fast python native file loading. Much faster but not completely tested.

        Author: SMM and DAV

        """
        fig = plt.figure()
        self.ax_list = []
        ax = fig.add_axes([0,0,1,1])
        self.ax_list.append(ax)

        # check the colourbar location
        if colourbar_location == "Top" or colourbar_location == "TOP":
            colourbar_location = "top"
        if colourbar_location == "Bottom" or colourbar_location == "BOTTOM":
            colourbar_location = "bottom"
        if colourbar_location == "Left" or colourbar_location == "LEFT":
            colourbar_location = "left"
        if colourbar_location == "Right" or colourbar_location == "RIGHT":
            colourbar_location = "right"

        print("Your colourbar will be located: "+ colourbar_location)
        if colourbar_location == "top" or colourbar_location == "bottom":
            self.colourbar_location = colourbar_location
            self.colourbar_orientation = "horizontal"
        elif colourbar_location == "left" or colourbar_location == "right":
            self.colourbar_location = colourbar_location
            self.colourbar_orientation = "vertical"
        elif colourbar_location == "None":
            self.colourbar_location = "None"
            self.colourbar_orientation = "None"
        else:
            print("You did not choose a valid colourbar location. Options are top, left, right, bottom and None. Defaulting to None.")
            self.colourbar_location = "None"
            self.colourbar_orientation = "None"

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
        if basemap_colourmap == "gray":
            self._RasterList.append(BaseRaster(BaseRasterName,Directory, NFF_opti = NFF_opti, alpha = alpha))
        else:
            self._RasterList.append(BaseRaster(BaseRasterName,Directory, NFF_opti = NFF_opti, alpha = alpha))
            self._RasterList[-1].set_colourmap(basemap_colourmap)

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

        # A title if needed
        self.title = plot_title

    def SetCustomExtent(self,xmin,xmax,ymin,ymax):
        """
        This function sets the plot extent in map coordinates and remakes the axis ticks

        Args:
          xmin: the minimum extent in easting
          xmax: the maximum extent in easting
          ymin: the minimum extent in northing
          ymax: the maximum extent in northing

        Author: MDH
        """
        # Get the tick properties
        self._xmin = xmin
        self._ymin = ymin
        self._xmax = xmax
        self._ymax = ymax
        self.make_ticks()

        # Annoying but the scatter plot resets the extents so you need to reassert them
        self.ax_list[0].set_xlim(self._xmin,self._xmax)
        self.ax_list[0].set_ylim(self._ymin,self._ymax)
        self.ax_list = self.make_base_image(self.ax_list)


    def make_ticks(self):
        """
        This function makes the tick marks and the tick labels.
        It has been optimised so you get nice looking ticks, so you shouldn't have to mess with them after this is called.

        Author: SMM
        """

        if self._coord_type == "UTM":
            self.tick_xlocs,self.tick_ylocs,self.tick_x_labels,self.tick_y_labels = LSDP.GetTicksForUTMNoInversion(self._BaseRasterFullName,self._xmax,self._xmin,
                             self._ymax,self._ymin,self._n_target_ticks,self.min_tick_spacing)
        elif self._coord_type == "UTM_km":
            self.tick_xlocs,self.tick_ylocs,self.tick_x_labels,self.tick_y_labels = LSDP.GetTicksForUTMNoInversion(self._BaseRasterFullName,self._xmax,self._xmin,
                             self._ymax,self._ymin,self._n_target_ticks,minimum_tick_spacing=1000)
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
        """
        This is a low level function for placing the ticks on the active image axis.

        Args:
            Ax (object): The axis object

        Author: SMM
        """
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

        Author: SMM
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
        elif axis_style == "Ultra_Thin":
            lw = 0.4
            ftsz = 4
            tpd = 0.3
            label_ftsz = 6
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
        """
        This function creates the base image. It creates the axis for the base image,
        further drapes and point data are placed upon this image.

        Args:
            ax_list: A list of axes, we append the base raster to the [0] element
                of the axis

        Author: DAV and SMM
        """

        # We need to initiate with a figure
        #self.ax = self.fig.add_axes([0.1,0.1,0.7,0.7])

        print("This colourmap is: "+ self._RasterList[0]._colourmap)
        im = self.ax_list[0].imshow(self._RasterList[0]._RasterArray, self._RasterList[0]._colourmap, extent = self._RasterList[0].extents, interpolation="nearest", alpha = self._RasterList[0]._alpha)

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
                        colorbarlabel = "Colourbar", discrete_cmap=False, n_colours=10,
                        norm = "None",
                        colour_min_max = [],
                        modify_raster_values=False,
                        old_values=[], new_values=[], cbar_type=float,
                        NFF_opti = False, custom_min_max = [], zorder=1):
        """
        This function adds a drape over the base raster.

        Args:
            RasterName (string): The name of the raster (no directory, but need extension)
            Directory (string): directory of the data
            colourmap (string or colourmap): The colourmap. Can be a strong for default colourmaps
            alpha (float): The transparency of the drape (1 is opaque, 0 totally transparent)
            show_colourbar (bool): True to show colourbar
            colourbarlabel (string): The label of the colourbar
            discrete_cmap (bool): If true, make discrete values for colours, otherwise a gradient.
            n_colours (int): number of colours in discrete colourbar
            norm (string): Normalisation of colourbar. I don't understand this so don't change
            colour_min_max( list of int/float): if it contains two elements, map the colourbar between these two values.
            modify_raster_values (bool): If true, it takes old_values in list and replaces them with new_values
            old_values (list): A list of values to be replaced in raster. Useful for masking and renaming
            new_values (list): A list of the new values. This probably should be done with a map: TODO
            cbar_type (type): Sets the type of the colourbar (if you want int labels, set to int)
            NFF_opti (bool): If true, uses the new file loading functions. It is faster but hasn't been completely tested.
            custom_min_max (list of int/float): if it contains two elements, recast the raster to [min,max] values for display.

        Author: SMM
        """
        print("N axes are: "+str(len(self.ax_list)))
        print(self.ax_list[0])

        self.ax_list = self._add_drape_image(self.ax_list,RasterName,Directory,colourmap,alpha,
                                             colorbarlabel,discrete_cmap,n_colours, norm,
                                             colour_min_max,modify_raster_values,old_values,
                                             new_values,cbar_type, NFF_opti, custom_min_max, zorder=zorder)
        #print("Getting axis limits in drape function: ")
        #print(self.ax_list[0].get_xlim())


    def _add_drape_image(self,ax_list,RasterName,Directory,
                         colourmap = "gray",
                         alpha=0.5,
                         colorbarlabel = "Colourbar", discrete_cmap=False,
                         n_colours=10, nroma = "None",
                         colour_min_max = [],
                         modify_raster_values = False,
                         old_values=[], new_values = [], cbar_type=float,
                         NFF_opti = False, custom_min_max = [],zorder=1):
        """
        This function adds a drape over the base raster. It does all the dirty work
        I can't quite remember why I did it in two steps but I vaguely recall trying it in one step and it didn't work.

        Args:
            RasterName (string): The name of the raster (no directory, but need extension)
            Directory (string): directory of the data
            colourmap (string or colourmap): The colourmap. Can be a strong for default colourmaps
            alpha (float): The transparency of the drape (1 is opaque, 0 totally transparent)
            show_colourbar (bool): True to show colourbar
            colourbarlabel (string): The label of the colourbar
            discrete_cmap (bool): If true, make discrete values for colours, otherwise a gradient.
            n_colours (int): number of colours in discrete colourbar
            norm (string): Normalisation of colourbar. I don't understand this so don't change
            colour_min_max( list of int/float): if it contains two elements, map the colourbar between these two values.
            modify_raster_values (bool): If true, it takes old_values in list and replaces them with new_values
            old_values (list): A list of values to be replaced in raster. Useful for masking and renaming
            new_values (list): A list of the new values. This probably should be done with a map: TODO
            cbar_type (type): Sets the type of the colourbar (if you want int labels, set to int)
            NFF_opti (bool): If true, uses the new file loading functions. It is faster but hasn't been completely tested.
            custom_min_max (list of int/float): if it contains two elements, recast the raster to [min,max] values for display.

        Author: SMM
        """
        Raster = BaseRaster(RasterName,Directory, NFF_opti = NFF_opti)
        if modify_raster_values == True:
            Raster.replace_raster_values(old_values, new_values)

        if discrete_cmap == True:
            print("N colours: "+str(n_colours))
            colourmap = self.cmap_discretize(colourmap, n_colours)




        self._RasterList.append(Raster)
        self._RasterList[-1].set_colourmap(colourmap)
        # I am recasting the raster to custom extents
        if len(custom_min_max)!=0:
            if len(custom_min_max)== 2:
                print("I am setting customisable minimum and maximum values: %s ¦¦ %s" %(custom_min_max[0],custom_min_max[1]))
                self._RasterList[-1]._RasterArray[self._RasterList[-1]._RasterArray<custom_min_max[0]] = custom_min_max[0]
                self._RasterList[-1]._RasterArray[self._RasterList[-1]._RasterArray>custom_min_max[1]] = custom_min_max[1]
            else:
                print("I cannot customize your minimum and maximum because I don't understand your input. It should be [min,max] with min max as integers or floats")

        # We need to initiate with a figure
        #self.ax = self.fig.add_axes([0.1,0.1,0.7,0.7])
        if len(colour_min_max)!=0:
            if len(colour_min_max)== 2:
                print("I am setting customisable colourbar minimum and maximum values: %s ¦¦ %s" %(colour_min_max[0],colour_min_max[1]))
                im = self.ax_list[0].imshow(self._RasterList[-1]._RasterArray, self._RasterList[-1]._colourmap, extent = self._RasterList[0].extents,
                                 interpolation="nearest",alpha = alpha, norm = mpl.colors.Normalize(vmin=colour_min_max[0], vmax=colour_min_max[1]),zorder=zorder)
            else:
                print("I cannot customize your colour minimum and maximum because I don't understand your input. It should be [min,max] with min max as integers or floats")
        else:
            if(nroma != "None"):
                im = self.ax_list[0].imshow(self._RasterList[-1]._RasterArray, self._RasterList[-1]._colourmap, extent = self._RasterList[0].extents,
                                 interpolation="nearest",alpha = alpha, norm = nroma, zorder=zorder)
            else:
                im = self.ax_list[0].imshow(self._RasterList[-1]._RasterArray, self._RasterList[-1]._colourmap,
                                 extent = self._RasterList[0].extents, interpolation="nearest",alpha = alpha, zorder=zorder)
        # This affects all axes because we set share_all = True.
        #ax.set_xlim(self._xmin,self._xmax)
        #ax.set_ylim(self._ymin,self._ymax)
        self.ax_list[0] = self.add_ticks_to_axis(self.ax_list[0])
        self._drape_list.append(im)

        print("The number of axes are: "+str(len(self._drape_list)))

        if self.colourbar_orientation != "None":
            self.ax_list = self.add_colourbar(self.ax_list,im,self._RasterList[-1],
                                              colorbarlabel = colorbarlabel, discrete=discrete_cmap,
                                              n_colours=n_colours, cbar_type=cbar_type)


        return self.ax_list



    def add_basin_plot(self,RasterName,BasinInfoPrefix,Directory,
                         colourmap = "gray",alpha=0.5,
                         show_colourbar = "False", colorbarlabel = "Colourbar",
                         discrete_cmap=False, n_colours=10, cbar_type=float,
                         use_keys_not_junctions = True,
                         label_basins = True, adjust_text = False, rename_dict = {},
                         value_dict = {}, mask_list = [],
                         edgecolour='black', linewidth=1, cbar_dict = {}, parallel=False,
                         outlines_only = False,zorder = 1):
        """
        This is a basin plotting routine. It plots basins as polygons which
        can be coloured and labelled in various ways.

        Args:
            RasterName (string): The name of the raster (no directory, but need extension)
            BasinInfoPrefix (string): The prefix (before "_BasinInfo.csv") of the basin info file produced by the chi_mapping_tool
            Directory (string): directory of the data
            colourmap (string or colourmap): The colourmap. Can be a strong for default colourmaps
            alpha (float): The transparency of the drape (1 is opaque, 0 totally transparent)
            show_colourbar (bool): True to show colourbar
            colourbarlabel (string): The label of the colourbar
            discrete_cmap (bool): If true, make discrete values for colours, otherwise a gradient.
            n_colours (int): number of colours in discrete colourbar
            cbar_type (type): Sets the type of the colourbar (if you want int labels, set to int)
            use_keys_not_junctions (bool): If true, the basin keys rather than the junction indices are used to map to the basins. f false the junction indices are used.
            label_basins (bool): If true, add text labels to basins.
            adjust_text (bool): If true calls the text adjustment routine. Takes a long time!
            rename_dict (dict): a dictionary where the key is the basin to rename (either key or junc, depending on use_keys_not_junctions) and the value is a string of the new name.
            value_dict (dict): the key is the basin (either key or junc, depending on use_keys_not_junctions) and the value is a new value that is used as a colour for the basin.
            mask_list (list of ints): Any basin named in this list (can be either a key or junction index depending on use_keys_not_junctions) is removed from the polgons and not plotted.
            edgecolour (string): colour of the lines around the basins.
            linewidth(float): width of the line around the basins.
            cbar_dict (dict): an optional dictionary to set the min and max of the colourbars, where the key is the
            min and the max, and the values are what you want to set the colourbar to. Leave empty if you just want this
            to be the same as the value dict.
            parallel (bool): option flag for processing multiple basin raster files triggered by parallel chi mapping tool.
            outlines_only (bool): If true, only plot the outlines.

        Author: SMM
        """

        from shapely.geometry import Polygon, Point
        from descartes import PolygonPatch
        from matplotlib.collections import PatchCollection

        # Get the basin outlines
        # Basins referes to a dict where the key is the junction index and the
        # value is a shapely polygon object
        if not parallel:
          Basins = LSDP.GetBasinOutlines(Directory, RasterName)
        else:
          Basins = LSDP.GetMultipleBasinOutlines(Directory)

        # Now check if you want to mask the basins
        # get the basin IDs to make a discrete colourmap for each ID
        #load the file
        if not parallel:
          BasinInfoDF = phelp.ReadBasinInfoCSV(Directory, BasinInfoPrefix)
        else:
          BasinInfoDF = phelp.AppendBasinInfoCSVs(Directory,BasinInfoPrefix)

        # Extract the basin keys
        basin_keys = list(BasinInfoDF['basin_key'])
        basin_keys = [int(x) for x in basin_keys]

        # Extract the junction indices
        basin_junctions = list(BasinInfoDF['outlet_junction'])
        basin_junctions = [int(x) for x in basin_junctions]

        # we need dicts for refering to each of these
        key_to_junction_dict = dict(zip(basin_keys,basin_junctions))
        junction_to_key_dict= dict(zip(basin_junctions,basin_keys))

        # Now mask some basins.
        # This has a lot of tedious logical statements to ensure there are
        # no errors associated with looking for a key in a dict that doens't exist
        for basin in mask_list:
            if use_keys_not_junctions:
                if basin in key_to_junction_dict:
                    this_junc = key_to_junction_dict[basin]
                    if this_junc in Basins:
                        del Basins[this_junc]
                    else:
                        print("I'm trying to mask a key, " +str(basin)+ ", which has junction " +str(this_junc)+ " that isn't there.")
            else:
                if basin in Basins:
                    del Basins[basin]
                else:
                    print("I'm trying to mask a basin, " +str(basin)+ " that isn't there.")

        # Now label the basins
        if label_basins:
            # This will hold the labels. Need to initiate here to ensure it lives outside control statements
            texts = []

            # First get the points
            Points = {}
            print("The number of basins are: "+str(len(Basins)))
            for basin_key, basin in Basins.items():
                Points[basin_key] = Point(basin.representative_point())
            print("The number of points are: "+str(len(Points)))

            # Now check if there is a renaming dictionary
            if len(rename_dict) == 0:
                if use_keys_not_junctions:
                    texts = self.add_text_annotation_from_shapely_points_v2(Points, text_colour='k', label_dict=junction_to_key_dict,zorder = zorder +1)
                else:
                    texts = self.add_text_annotation_from_shapely_points_v2(Points, text_colour='k',zorder = zorder +1)
            else:
                # Okay so the way tyhis is going to work is that we ware going to
                # look for renamed basins but basins that haven't been renamed are
                # going to get their old names.

                # First we see if the renamed basins are in the lists
                if use_keys_not_junctions:
                    # We start with a clean dict that will be used for labelling
                    new_label_dict = junction_to_key_dict

                    # We are using keys. We need to replace the label dict with
                    # with the strings from the renamed basins
                    for key in rename_dict:
                        print("I am renaming. The key is: " +str(key))
                        # get the junction number of this key
                        if key in key_to_junction_dict:
                            this_junc = key_to_junction_dict[key]
                            print("The junction is: "+ str(this_junc))
                            new_label_dict[this_junc] = rename_dict[key]
                        else:
                            print("I am missing this key")

                    # Use this new label dict to rename the junctions
                    texts = self.add_text_annotation_from_shapely_points_v2(Points, text_colour='k', label_dict=new_label_dict,zorder = zorder +1)

                else:
                    # Now we need a new junction dict for this
                    new_label_dict = {}
                    for junction in junction_to_key_dict:
                        new_label_dict[junction] = junction

                    for key in rename_dict:
                        # get the junction number of this key
                        if key in new_label_dict:
                            new_label_dict[key] = rename_dict[key]

                    # Use this new label dict to rename the junctions
                    texts = self.add_text_annotation_from_shapely_points_v2(Points, text_colour='k', label_dict=new_label_dict,zorder = zorder +1)

            if adjust_text == True:
                print("I am adjusting the text for you. Warning: this takes a long time!")
                LSDP.adjust_text(texts)
                print("Finished adjusting text.")


        # Get the appropriate colourmap
        if type(colourmap) == str:
            this_cmap=plt.get_cmap(colourmap)
        else:
            this_cmap=colourmap

        """
        Some notes on colormapping:
            If the value dict is present, you colour by the value dict.
            The value dict has
            If the value dict is empty, you colour using a discrete cmap
        """

        min_value = 0
        max_value = n_colours-1
        if len(value_dict) == 0:
            # There are no values so we ensure there is a discrete colourmap
            discrete_cmap = True
            this_cmap = self.cmap_discretize(this_cmap, n_colours)

            # The colours go from 0 to number of colours minus 1
            min_value = 0
            max_value = n_colours-1
            cNorm  = colors.Normalize(min_value, max_value)
            new_colours = plt.cm.ScalarMappable(norm=cNorm, cmap=this_cmap)

            # now plot the polygons
            print('Plotting the polygons, colouring by basin...')
            for key, poly in Basins.items():
                print(key, int(key))
                colourkey = int(key) % n_colours
                # We need two patches since we don't want the edges transparent
                if not outlines_only:
                    this_patch = PolygonPatch(poly, fc=new_colours.to_rgba(colourkey), ec="none", alpha=alpha,zorder = zorder)
                    self.ax_list[0].add_patch(this_patch)
                this_patch = PolygonPatch(poly, fc="none", ec=edgecolour, alpha=1,zorder = zorder)
                self.ax_list[0].add_patch(this_patch)
        else:
            if discrete_cmap:
                this_cmap = self.cmap_discretize(this_cmap, n_colours)

            # Now we need to normalise from the values, so we go from the minimum
            # To the maximum
            # check if you want to renormalise the colourbar based on specific values
            if len(cbar_dict) != 0:
                min_value = cbar_dict.get('min')
                max_value = cbar_dict.get('max')
            else:
                max_value = max([i for i in value_dict.values()])
                min_value = min([i for i in value_dict.values()])

            # Normalise the colourmap
            cNorm  = colors.Normalize(min_value, max_value)
            new_colours = plt.cm.ScalarMappable(norm=cNorm, cmap=this_cmap)

            # we need a grayed out value for basins that don't have a value
            gray_colour = "#a9a9a9"

            # now plot the polygons
            print(Basins)
            print(junction_to_key_dict)
            print('Plotting the polygons, colouring by value...')
            for junc, poly in Basins.items():

                # If we are using keys, we need to check to see if the key referred to by
                # this junction is in the value dict
                if not outlines_only:
                    if use_keys_not_junctions:
                        print(junc)
                        this_key = junction_to_key_dict[int(junc)]
                        #print ("This key is: "+str(this_key)+", and this value is: "+str(value_dict[this_key]))
                        if this_key in value_dict:
                            this_patch = PolygonPatch(poly, fc=new_colours.to_rgba( value_dict[this_key] ), ec="none", alpha=alpha, zorder = zorder)
                            self.ax_list[0].add_patch(this_patch)
                        else:
                            this_patch = PolygonPatch(poly, fc=gray_colour, ec="none", alpha=alpha, zorder = zorder)
                            self.ax_list[0].add_patch(this_patch)
                    else:
                        # We are using junction indices so these link directly in to the polygon keys
                        if junc in value_dict:
                            this_patch = PolygonPatch(poly, fc=new_colours.to_rgba( value_dict[junc] ), ec="none", alpha=alpha, zorder = zorder)
                            self.ax_list[0].add_patch(this_patch)
                        else:
                            this_patch = PolygonPatch(poly, fc=gray_colour, ec="none", alpha=alpha, zorder = zorder)
                            self.ax_list[0].add_patch(this_patch)

                # We need to add the outline seperately because we don't want it to be transparent
                this_patch = PolygonPatch(poly, fc="none", ec=edgecolour, alpha=1, zorder = zorder)
                self.ax_list[0].add_patch(this_patch)



        # Now plot the colourbar
        if show_colourbar:
            print("The colourbar orientation for basin plotting is: "+self.colourbar_orientation)
            if self.colourbar_orientation != "None":
                print("Let me add a colourbar for your basin data")
                print("CBAR TYPE IS", cbar_type)
                this_cbar_type = cbar_type
                this_cbarlabel = colorbarlabel
                self.ax_list = self.add_objectless_colourbar(self.ax_list,
                                                         min_value, max_value,
                                                         cmap = this_cmap,colorbarlabel = this_cbarlabel,
                                                         discrete=discrete_cmap, n_colours=n_colours, cbar_type=this_cbar_type)

            else:
                print("You have asked me to plot a colourbar but the location is none. Give the mapfigure object a colourbar location.")



    def cmap_discretize(self, cmap, N):
        """
        Return a discrete colormap from the continuous colormap cmap.
        From http://scipy.github.io/old-wiki/pages/Cookbook/Matplotlib/ColormapTransformations

        Args:
             cmap: colormap instance, eg. cm.jet.
             N: number of colors.

        Returns:
            discrete colourmap

        Author: FJC

        """

        if type(cmap) == str:
         cmap = plt.get_cmap(cmap)

        colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
        colors_rgba = cmap(colors_i)
        indices = np.linspace(0, 1., N+1)
        cdict = {}

        for ki,key in enumerate(('red','green','blue')):
         cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in xrange(N+1) ]

        # Return colormap object.
        return _mcolors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)

    def add_colourbar(self,ax_list,im,BaseRaster,colorbarlabel = "Colourbar",discrete=False, n_colours=10, cbar_type=float):
        """
        This adds the colourbar to the image.
        IMPORTANT: It assumes the colourbar occupies the last axis element

        Args:
            ax_list: The list of axes objects. Assumes colourbar is in axis_list[-1]
            im: The image object
            BaseRaster: The base raster.
            colorbarlabel (string): The label of the colourbar
            discrete_cmap (bool): If true, make discrete values for colours, otherwise a gradient.
            n_colours (int): number of colours in discrete colourbar.
            cbar_type (type): Sets the type of the colourbar (if you want int labels, set to int).

        Author: SMM
        """
        fig = matplotlib.pyplot.gcf()
        ax_list.append(fig.add_axes([0.1,0.8,0.05,0.2]))
        cbar = plt.colorbar(im,cmap=BaseRaster._colourmap,spacing='uniform', orientation=self.colourbar_orientation,cax=ax_list[-1])

        if discrete==True:
            # change ticks
            self.fix_colourbar_ticks(BaseRaster, cbar, n_colours, cbar_type)

        #Will's changes:
        # Changed rotation of colourbar text to 90 and the labelpad to -75 for "left"

        if self.colourbar_location == 'top':
            ax_list[-1].set_xlabel(colorbarlabel, fontname='Arial',labelpad=5)
        elif self.colourbar_location == 'bottom':
            ax_list[-1].set_xlabel(colorbarlabel, fontname='Arial',labelpad=5)
        elif self.colourbar_location == 'left':
            ax_list[-1].set_ylabel(colorbarlabel, fontname='Arial',labelpad=-75,rotation=90)
        elif self.colourbar_location == 'right':
            ax_list[-1].set_ylabel(colorbarlabel, fontname='Arial',labelpad=10,rotation=270)
        return ax_list

    def fix_colourbar_ticks(self, BaseRaster, cbar,n_colours, cbar_type=float,
                            use_baseraster = True, min_value = 0, max_value = 0, cbar_label_rotation=30):
        """
        This function takes a discrete colourbar and fixes the ticks so they are
        in the middle of each colour

        Update: SMM, made this more flexible so you can set the minimum and maximum
         values. Required becase basin plotting and point plotting doesn't use
         raster values
        Update 2: BG, somehow, writing "update" in upper cases make my text editors (Atom and brackets) completely panic.
        so I changed it, sorry.

        Args:
            BaseRaster: the base raster object
            cbar: the colourbar
            n_colours (int): number of colours in discrete colourbar.
            cbar_type (type): Sets the type of the colourbar (if you want int labels, set to int).
            use_baseraster (bool): True if you want to use the baseraster, otherwise uses the min_value and max_value
            min_value: the minimum value on the colourbar
            max_value: the maximum value on the colourbar
            cbar_label_rotation (float): rotate the tick labels

        Returns:
            None but fixes ticks

        Author: FJC
        """

        # get the min and the max of the colourbar
        if use_baseraster:
            vmin = np.nanmin(BaseRaster._RasterArray)
            vmax = np.nanmax(BaseRaster._RasterArray)
        else:
            print("I'm fixing the ticks, but won't use a base raster, since you told me not to.")
            vmin = min_value
            vmax = max_value
        print("The min and max for the colourbar are:")
        print((vmin, vmax))
        print(n_colours)
        print(cbar_type)

        # get the additional end spacing for colourbar
        tick_spacing = float(vmax-vmin)/float(n_colours)
        print(tick_spacing)
        new_vmin = vmin-(tick_spacing/2)
        new_vmax = vmax+(tick_spacing/2)+tick_spacing

        #get list of tick locations
        tick_locs = np.arange(new_vmin, new_vmax, step=tick_spacing)
        print(tick_locs)

        # update ticks
        tick_locator = ticker.FixedLocator(tick_locs)
        cbar.locator = tick_locator
        cbar.update_ticks()

        #print BaseRaster._RasterArray[np.isnan(BaseRaster._RasterArray) == False]

        # get tick labels
        tick_labels = np.linspace(vmin, vmax, n_colours)
        if cbar_type == int:
            tick_labels = [str(int(x)) for x in tick_labels]
        else:
            tick_labels = [str(x) for x in tick_labels]
        print(tick_labels)

        if self.colourbar_orientation == "horizontal":
            cbar.ax.set_xticklabels(tick_labels, rotation=cbar_label_rotation)
        else:
            cbar.ax.set_yticklabels(tick_labels, rotation=cbar_label_rotation)

    def add_point_colourbar(self,ax_list,sc,cmap = "cubehelix",colorbarlabel = "Colourbar",
                            discrete=False, n_colours=10, cbar_type=float):
        """
        This adds a colourbar for any points that are on the DEM.

        Args:
            ax_list: The list of axes objects. Assumes colourbar is in axis_list[-1]
            sc: The scatterplot object. Generated by plt.scatter
            cmap (string or colourmap): The colourmap.
            colorbarlabel (string): The label of the colourbar

        Author: SMM
        """
        fig = matplotlib.pyplot.gcf()
        ax_list.append(fig.add_axes([0.1,0.8,0.2,0.5]))
        cbar = plt.colorbar(sc,cmap=cmap, orientation=self.colourbar_orientation,cax=ax_list[-1])

        if self.colourbar_location == 'top':
            ax_list[-1].set_xlabel(colorbarlabel, fontname='Arial',labelpad=5)
        elif self.colourbar_location == 'bottom':
            ax_list[-1].set_xlabel(colorbarlabel, fontname='Arial',labelpad=5)
        elif self.colourbar_location == 'left':
            ax_list[-1].set_ylabel(colorbarlabel, fontname='Arial',labelpad=-75,rotation=90)
        elif self.colourbar_location == 'right':
            ax_list[-1].set_ylabel(colorbarlabel, fontname='Arial',labelpad=10,rotation=270)
        return ax_list


    def add_objectless_colourbar(self,ax_list,minimum_value, maximum_value,
                                 cmap = "cubehelix",colorbarlabel = "Colourbar",
                                 discrete=False, n_colours=10, cbar_type=float):
        """
        This adds a colourbar that is not attached to any particular object

        Args:
            ax_list: The list of axes objects. Assumes colourbar is in axis_list[-1]
            minimum_value (float or int): minimum value on colourbar
            maximum_value (float or int): maximum value on colourbar
            cmap (string or colourmap): The colourmap.
            colorbarlabel (string): The label of the colourbar

        Returns:
            The axis list

        Author: SMM
        """
        fig = matplotlib.pyplot.gcf()
        ax_list.append(fig.add_axes([0.1,0.8,0.2,0.5]))
        cnorm = colors.Normalize(minimum_value, maximum_value)
        this_cmap = cmap
        cbar = mpl.colorbar.ColorbarBase(ax_list[-1], cmap=this_cmap,
                                norm=cnorm,
                                orientation=self.colourbar_orientation)

        if discrete==True:
            # change ticks
            self.fix_colourbar_ticks(BaseRaster, cbar, n_colours, cbar_type, False, minimum_value, maximum_value)
            cbar_label_rotation=30


        #Will's changes:
        # Changed rotation of colourbar text to 90 and the labelpad to -75 for "left"

        if self.colourbar_location == 'top':
            ax_list[-1].set_xlabel(colorbarlabel, fontname='Arial',labelpad=5)
        elif self.colourbar_location == 'bottom':
            ax_list[-1].set_xlabel(colorbarlabel, fontname='Arial',labelpad=5)
        elif self.colourbar_location == 'left':
            ax_list[-1].set_ylabel(colorbarlabel, fontname='Arial',labelpad=-75,rotation=90)
        elif self.colourbar_location == 'right':
            ax_list[-1].set_ylabel(colorbarlabel, fontname='Arial',labelpad=10,rotation=270)
        return ax_list


    def add_point_data(self, thisPointData,column_for_plotting = "None",
                       this_colourmap = "cubehelix", show_colourbar = False, colourbar_location = "bottom",
                       colorbarlabel = "Colourbar",
                       scale_points = False,column_for_scaling = "None",
                       scaled_data_in_log = False,
                       max_point_size = 5, min_point_size = 0.5,
                       colour_log = False, colour_manual_scale = [],
                       manual_size = 0.5, alpha = 1, minimum_log_scale_cut_off = -10, label_field = "None",
                       font_size = 6, offset = 100, zorder=1, marker = "o", discrete_colours = False, NColours = 10,scale_in_absolute = False, color_abs =False, unicolor = "blue",
                       recast_scale_min_max = [], scale_in_abs_after_recasting = False):

        """
        This add point data to the map.

        Args:
            thisPointData (object): an LSDMap_PointData object.
            column_for_plotting (string): The column in the pint data that is used for plotting the points.
            this_colourmap (string or colourmap): The colourmap.
            show_colourbar (bool): If true, will display a colourbar for the points
            colourbar_location (string): The location of the colourbar, can be either "top", "bottom", "left", or "right"
            colorbarlabel (string): The label of the colourbar
            scale_point (bool): If true, point size is scaled by the point value.
            column_for_scaling (string): The column name that is used to scale the point size
            scaled_data_in_log (bool): If true, the points are scaled in proportion to the logarithm of their value.
            max_point_size (float): Maximum size in points of the symbols.
            min_point_size (float): Minumum size in points of the symbols.
            colour_log (bool): If the colours are scaled by logarithm.
            colour_manual_scale (list): A two element list containing the minimum and maximum values of colourbar (i.e. if you want to cut off the colurbar at a certain value).
            manual_size (float): If scale_points is false then this is the size of the points.
            alpha (float): transparency (between 0 and 1).
            minimum_log_scale_cut_off (float): If the log of the value is less than this the point is not plotted.
            label_field (str): text annotation below the point contained in this column in the csv file
            offset (int/float): offset of the text below the point
            font_size (int): everything is in the title
            zorder (int): priority for plotting
            marker (str): the marker used in the plots
            discrete_colours (bool): If true, the colourmap will be discrete
            NColours (int) The number of colours n the colourmap
            scale_in_absolute (bool): scale the data using absolute values
            abs (bool): color the data using absolute values
            unicolor (str): set a unique color in case of no plotting column - Default: blue
            color_abs: get the absolute data for scale
            recast_scale_min_max: recast the min and max of the array before scaling
            scale_in_abs_after_recasting: give the abolute value of scaling after recasting the data

        Author: SMM, BG
        """


        # Get the axis limits to assert after
        this_xlim = self.ax_list[0].get_xlim()
        this_ylim = self.ax_list[0].get_ylim()

        EPSG_string = self._RasterList[0]._EPSGString
        print("I am going to plot some points for you. The EPSG string is:"+EPSG_string)


        # convert to easting and northing or pull easting northing from file
        # I had an old file that didnt report lat/long so pull directly if lat/lon not found
        # MDH 1/3/18
        try:
            [easting,northing] = thisPointData.GetUTMEastingNorthing(EPSG_string)
        except:
            # check to see if easting and northing data already exists
            easting = thisPointData.QueryData("easting").as_matrix().astype(float)
            northing = thisPointData.QueryData("northing").as_matrix().astype(float)

        print("I got the easting and northing")

        # check if the column for plotting exists
        # BG - 16/01/2018 - Adding some exception management. Sometimes, a list can be returned by QueryData and crash here
        try:
            this_data = thisPointData.QueryData(column_for_plotting).values
        except AttributeError:
            this_data = thisPointData.QueryData(column_for_plotting)

        print("I got the data column you wanted")
        if(color_abs):
            print("I will color your data using its absolute value")
            this_data = np.abs(this_data)
        # Log the color if required
        if(colour_log):
            this_data = np.log10(this_data)
            print("I have taken the log your colour data, the minimum is %s and the maximum is %s" %(np.nanmin(this_data), np.nanmax(this_data)))

        # Now the data for scaling. Point size will be scaled by these data
        scale_data = thisPointData.QueryData(column_for_scaling)
        print("I also got the data for scaling, which is in column "+column_for_scaling)
        scale_data = np.asarray(scale_data)
        if(scale_in_absolute):
            scale_data = np.abs(scale_data)
        #scale_data = scale_data.flatten()
        print("The size of the array is: ")
        print(scale_data.shape)

        # If there is scaled data, convert to log if that option is selected
        if scaled_data_in_log:
            print("I am going to convert data to log for point scaling.")
            if len(scale_data) == 0 or len(scale_data) != len(easting):
                scale_data = [0.5]
            else:
                # We need this logic since we can get nans and -Infs from 0 and negative numbers
                scale_data = np.log10(scale_data)
                print("I logged (is it a verb?) your scaled data, the minimum is %s and the maximum is %s but all the values inferior to %s will be %s" %(np.nanmin(scale_data), np.nanmax(scale_data), minimum_log_scale_cut_off, minimum_log_scale_cut_off))
                scale_data[scale_data < minimum_log_scale_cut_off] = minimum_log_scale_cut_off
        else:
            print("You are not going to use a log scale to scale the size of the points")


        # scale the points if you want
        if scale_points == True:
            print("I am scaling your points for you")
            if len(scale_data) == 0 or len(scale_data) != len(easting):
                print("There doesn't seem to be any scaling data. Reverting to manual size.")
                point_scale = manual_size
            else:

                if(scale_in_abs_after_recasting):
                    scale_data = np.abs(scale_data)

                if(len(recast_scale_min_max) == 2):
                    scale_data[scale_data < recast_scale_min_max[0]] = recast_scale_min_max[0]
                    scale_data[scale_data > recast_scale_min_max[1]] = recast_scale_min_max[1]



                max_sd = np.nanmax(scale_data)
                min_sd = np.nanmin(scale_data)

                print("max is: "+str(max_sd)+ " and min is: "+ str(min_sd))

                # now rescale the data. Always a linear scaling.
                new_scale = []
                size_range = max_point_size-min_point_size
                for datum in scale_data:
                    frac = (datum-min_sd)/(max_sd-min_sd)
                    new_size = frac*size_range+min_point_size
                    new_scale.append(new_size)
                ns_array = np.asarray(new_scale)
                point_scale = ns_array
                print("I have got a scaled point array,")
        else:
            print("I will not scale your points.")
            point_scale = manual_size



        print("I will plot the points now.")
        if len(this_data) == 0 or len(this_data) != len(easting):
            print("I am only plotting the points.")
            unicolor = unicolor
            sc = self.ax_list[0].scatter(easting,northing,s=point_scale, c= unicolor,cmap=this_colourmap,edgecolors='none', alpha = alpha,zorder=zorder, marker = marker)
        else:
            print("I will colour by the points")
            if(colour_manual_scale != []):
                print("let me rescale the colour using your array")
                if(len(colour_manual_scale) == 2):
                    cNorm  = _mcolors.Normalize(vmin=colour_manual_scale[0], vmax=colour_manual_scale[1])
                    #scalarMap = _cm.ScalarMappable(norm = cNorm, cmap= this_colourmap)
                    #tps_color = scalarMap.to_rgba(this_data)
                    #scalarMap.set_array(tps_color)
                    #this_colourmap = scalarMap
                    #sc = self.ax_list[0].scatter(easting,northing,s=point_scale, c=tps_color,cmap=this_colourmap,edgecolors='none', alpha = alpha)
                    sc = self.ax_list[0].scatter(easting,northing,s=point_scale, c=this_data,cmap=this_colourmap,norm=cNorm,edgecolors='none', alpha = alpha,zorder=zorder, marker = marker)

                else:
                    print("Your colour_log_manual_scale should be something like [min,max], aborting")
                    quit()
            else:
                #sc = self.ax_list[0].scatter(easting,northing,s=point_scale, c=this_data,cmap=this_colourmap,edgecolors='none', alpha = alpha,zorder=zorder, marker = marker)
                if discrete_colours:
                    # make a color map of fixed colors
                    NUM_COLORS = NColours

                    this_cmap = this_colourmap
                    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
                    plt.cm.ScalarMappable(norm=cNorm, cmap=this_colourmap)
                    channel_data = [x % NUM_COLORS for x in this_data]

                    sc = self.ax_list[0].scatter(easting,northing,s=point_scale, c=channel_data,cmap=this_colourmap,norm=cNorm,edgecolors='none', alpha = alpha,zorder=zorder)
                else:
                    sc = self.ax_list[0].scatter(easting,northing,s=point_scale, c=this_data,cmap=this_colourmap,edgecolors='none', alpha = alpha,zorder=zorder, marker = marker)

        # Setting the labelling
        if(label_field != "None"):
            print("labelling from this tool is not available yet, Boris is working on it")
            tg = thisPointData.QueryData(label_field)
            print(tg)
            for i in range(len(easting)):
                print(str(tg[i]))
                sc =self.ax_list[0].text(easting[i]-offset,northing[i]-offset,str(tg[i]),fontsize = font_size)

        # Annoying but the scatter plot resets the extents so you need to reassert them
        self.ax_list[0].set_xlim(this_xlim)
        self.ax_list[0].set_ylim(this_ylim)

        if show_colourbar == True:
            print("Your colourbar will be located: "+ colourbar_location)
            if colourbar_location == "top" or colourbar_location == "bottom":
                self.colourbar_location = colourbar_location
                self.colourbar_orientation = "horizontal"
            elif colourbar_location == "left" or colourbar_location == "right":
                self.colourbar_location = colourbar_location
                self.colourbar_orientation = "vertical"
            elif colourbar_location == "None":
                self.colourbar_location = "None"
                self.colourbar_orientation = "None"
            if self.colourbar_location != "None":
                print("Let me add a colourbar for your point data")
                self.ax_list = self.add_point_colourbar(self.ax_list,sc,cmap=this_colourmap, colorbarlabel = colorbarlabel)

    def add_line_data(self, ThisLineFile, linestyle = '-', edgecolour = "k", linewidth=0.5, zorder = 1, alpha=0):
        """
        This adds line data from a named shapefile to the map.

        Args:
            ThisLineFile (object): a shapefile for plotting lines.
            linestyle: matplotlib line style
            edgecolour (string): colour of the lines around the basins.
            linewidth(float): width of the line around the basins.
            zorder (int): priority for plotting
            alpha (float): transparency (between 0 and 1).

        Author: MDH
        """

        # import linestring
        from shapely.geometry import shape, LineString
        import fiona

        # Get the axis limits to assert after
        this_xlim = self.ax_list[0].get_xlim()
        this_ylim = self.ax_list[0].get_ylim()

        EPSG_string = self._RasterList[0]._EPSGString
        print("I am going to plot some lines for you. The EPSG string is:"+EPSG_string)

        # load the data
        with fiona.open(ThisLineFile) as input:
            for feature in input:
                geom = shape(feature['geometry'])
                x, y = geom.xy
                self.ax_list[0].plot(x,y,linestyle,color=edgecolour,lw=linewidth, zorder=zorder, alpha=alpha)

        # Annoying but the scatter plot resets the extents so you need to reassert them
        self.ax_list[0].set_xlim(this_xlim)
        self.ax_list[0].set_ylim(this_ylim)

    def plot_segment_of_knickzone(self, thisPointData, color = "k", lw = 1):
        # Get the axis limits to assert after
        this_xlim = self.ax_list[0].get_xlim()
        this_ylim = self.ax_list[0].get_ylim()

        EPSG_string = self._RasterList[0]._EPSGString
        print("I am going to plot some points for you. The EPSG string is:"+EPSG_string)

        # convert to easting and northing
        [easting,northing] = thisPointData.GetUTMEastingNorthing(EPSG_string)
        print("I got the easting and northing")

        if(len(easting)>1):
            self.ax_list[0].plot(easting,northing, c = color, lw = lw)
        elif (len(easting) == 1):

            self.ax_list[0].scatter(easting,northing, c = color, s = lw * 0.2)

        # Annoying but the scatter plot resets the extents so you need to reassert them
        self.ax_list[0].set_xlim(this_xlim)
        self.ax_list[0].set_ylim(this_ylim)


    def add_channel_network_from_points(self, thisPointData, colour='k', alpha = 0.7, zorder=1):
        """
        This function plots the channel network from the map figure, you must pass in the
        channel network as an LSDMap_PointData object.

        Args:
            thisPointData (object): an LSDMap_PointData object
            colour (string): colour you want the channel network to be, default = black
            alpha (float): transparency, 1 = opaque. Default = 0.7
            zorder (float): priority for layering of plots

        Returns:
            plots the channel network

        Author: FJC
        """
        # Get the axis limits to assert after
        this_xlim = self.ax_list[0].get_xlim()
        this_ylim = self.ax_list[0].get_ylim()

        EPSG_string = self._RasterList[0]._EPSGString
        print("I am going to plot the channel network for you. The EPSG string is:"+EPSG_string)

        # convert to easting and northing
        [easting,northing] = thisPointData.GetUTMEastingNorthing(EPSG_string)
        print("I got the easting and northing")

        sc = self.ax_list[0].scatter(easting,northing,s=0.1, c=colour, facecolor=colour, alpha = alpha, zorder=zorder)


    def add_text_annotation_from_points(self, thisPointData,column_for_plotting = "None",
                                        selection_criteria = [], PANDEX=False, border_colour='k', text_colour='r', alpha=1):
        """
        This adds annotations to points. Used for annotating basins or sources, for example.

        Args:
            thisPointData (object): an LSDMap_PointData object.
            column_for_plotting (string): The column in the pint data that is used for plotting the points.
            selection_critera (list): This selects given values for plotting.
            PANDEX (bool): If true uses pandas data loading. Much faster but not fully tested.
            border_colour (string or colour): The colour of the edges around the textboxes.
            text_colour (string or colour): The colour of the text.
            alpha (float): The transparency of the text.

        Returns:
            A text annotation object

        Author: SMM
        """

        # A list of text objects
        texts = []

        # Get the axis limits to assert after
        this_xlim = self.ax_list[0].get_xlim()
        this_ylim = self.ax_list[0].get_ylim()

        # Format the bounding box
        bbox_props = dict(boxstyle="circle,pad=0.1", fc="w", ec=border_colour, lw=0.5,alpha = alpha)

        # see if the data column exists
        test_data = thisPointData.QueryData(column_for_plotting, PANDEX=PANDEX)

        if len(test_data) == 0:
            print("There is no data with the column name: "+column_for_plotting)
        else:

            # Thin the data
            if len(selection_criteria) == 1:
                thisPointData.ThinData(column_for_plotting,selection_criteria)
            elif len(selection_criteria) >1:
                thisPointData.ThinDataSelection(column_for_plotting,selection_criteria)

            # get the easting and northing
            EPSG_string = self._RasterList[0]._EPSGString
            print("EPSG_string is: "+EPSG_string)
            [this_easting,this_northing] = thisPointData.GetUTMEastingNorthing(EPSG_string)
            print((this_easting, this_northing))
            thinned_data = thisPointData.QueryData(column_for_plotting, PANDEX=PANDEX)
            print(thinned_data)

            for idx, data in enumerate(thinned_data):
                texts.append(self.ax_list[0].text(this_easting[idx],this_northing[idx], str(data),fontsize = 8, color= text_colour,alpha=alpha,bbox=bbox_props))
            #print ("I'm adding the text, yo")

        # Annoying but the scatter plot resets the extents so you need to reassert them
        self.ax_list[0].set_xlim(this_xlim)
        self.ax_list[0].set_ylim(this_ylim)

        return texts

    def add_text_annotation_from_shapely_points(self, points, label_dict={}, border_colour='k', text_colour='r', alpha=1,zorder=10):
        """
        This adds annotations from a dictionary of shapely points, for annotating basins or sources.


        Args:
            points: This is a dictionary the keys are the raster values to annotate and the values are the point objects.
            label_dict: The labels are also stored in a dictionary, where the key is the original value (e.g. basin
                junction, and the value is a string that you want to label with (e.g. the basin key).

        Returns:
            A text annotation object

        FJC 24/06/17
        """
        from shapely.geometry import Point

        # rewrite with new values if you need to (for basins)
        print(points)
        print(label_dict)

        new_points = {}
        if label_dict:
            for key, label in label_dict.items():
                # get the point for this key
                new_points[label] = points.get(key)
                print((key, label, new_points[label]))
            points = new_points

        # A list of text objects
        texts = []

        # Get the axis limits to assert after
        this_xlim = self.ax_list[0].get_xlim()
        this_ylim = self.ax_list[0].get_ylim()

        # Format the bounding box
        bbox_props = dict(boxstyle="Round4,pad=0.1", fc="w", ec=border_colour, lw=0.5,alpha = alpha)

        print(points)
        for key, point in points.items():
            x = point.x
            y = point.y
            texts.append(self.ax_list[0].text(point.x, point.y, str(key), fontsize=8, color=text_colour,alpha=alpha,bbox=bbox_props, ha= 'center',zorder=zorder))
            #print ("I'm adding the text, yo")

        # Annoying but the scatter plot resets the extents so you need to reassert them
        self.ax_list[0].set_xlim(this_xlim)
        self.ax_list[0].set_ylim(this_ylim)

        return texts

    def add_text_annotation_from_shapely_points_v2(self, points, label_dict={}, border_colour='k', text_colour='r', alpha=1,zorder=10):
        """
        This adds annotations from a dictionary of shapely points, for annotating basins or sources.

        THIS WORKS IN SMM's USE CASE
        I'm retaining the old one so not to break fionas code



        Args:
            points: This is a dictionary the keys are the raster values to annotate and the values are the point objects.
            label_dict: The labels are also stored in a dictionary, where the key is the original value (e.g. basin
                junction, and the value is a string that you want to label with (e.g. the basin key).

        Returns:
            A text annotation object

        SMM 24/06/17
        """
        from shapely.geometry import Point

        # A list of text objects
        texts = []

        # Get the axis limits to assert after
        this_xlim = self.ax_list[0].get_xlim()
        this_ylim = self.ax_list[0].get_ylim()

        # Format the bounding box
        bbox_props = dict(boxstyle="Round4,pad=0.1", fc="w", ec=border_colour, lw=0.5,alpha = alpha)

        for key, point in points.items():
            x = point.x
            y = point.y

            # If there is no label dict, just append with the text
            if len(label_dict) == 0:
                texts.append(self.ax_list[0].text(point.x, point.y, str(key), fontsize=8, color=text_colour,alpha=alpha,bbox=bbox_props, ha= 'center',zorder=zorder))
            else:
                if key in label_dict:
                    this_label = str(label_dict[key])
                    texts.append(self.ax_list[0].text(point.x, point.y, this_label, fontsize=8, color=text_colour,alpha=alpha,bbox=bbox_props, ha= 'center',zorder=zorder))

        # Annoying but the scatter plot resets the extents so you need to reassert them
        self.ax_list[0].set_xlim(this_xlim)
        self.ax_list[0].set_ylim(this_ylim)

        return texts

    def add_arrows_from_points(self, arrow_df, azimuth_header='azimuth', arrow_length=0.1, colour='black', linewidth=1, alpha = 1):
        """
        This function adds arrows at locations and azimuths specified by the
        pandas dataframe.  The X and Y locations should have column headers
        called 'X' and 'Y', the user specifies the column header of the azimuths
        column (default = 'azimuth')

        Args:
            arrow_df: pandas dataframe with the X and Y locations and the azimuths
            azimuth_header (str): the name of the column header with the azimuth locations
            colour: arrow colour, can pass in the name of one of the column headers and arrows will be coloured according to this.

        Author: FJC

        """
        import math
        import matplotlib.patches as patches

        # convert azimuths to radians
        azimuths = list(arrow_df[azimuth_header])
        az_radians = np.radians(azimuths)
        print(az_radians)

        # get points
        X = np.array(arrow_df['X'])
        Y = np.array(arrow_df['Y'])
        dx = np.array([ arrow_length*np.sin(i) for i in az_radians ])
        dy = np.array([ arrow_length*np.cos(i) for i in az_radians ])
        new_X = X - dx/2
        new_Y = Y - dy/2

        print((dx,dy))

        self.ax_list[0].quiver(new_X,new_Y,dx,dy,angles='xy',scale_units='xy',scale=1, width=0.002)

    def add_strike_and_dip_symbols(self,dip_df,colour='black',linewidth=1,alpha=1, symbol_length=10):
        """
        This function adds strike and dip symbols to the map from a pandas
        dataframe with the dip, dip direction, and strike
        Will be labelled with the dip.

        Args:
            dip_df: pandas dataframe with the dip and dip direction info
            colour: colour of the symbols
            linewidth: linewidth of the symbols
            alpha: transparency of the symbols.
            symbol_length: length of the strike element of the strike/dip symbol

        Author: FJC
        """
        # first we need get the strike to radians
        strikes = np.array(dip_df['strike'])
        strike_radians = np.radians(strikes)

        # get points
        X = np.array(dip_df['X'])
        Y = np.array(dip_df['Y'])
        # work out the dx and dy of the strike line
        dx = np.array([ symbol_length*np.sin(i) for i in strike_radians ])
        dy = np.array([ symbol_length*np.cos(i) for i in strike_radians ])
        # we want the centre of the line to be at the X/Y point
        start_X = X - dx/2
        start_Y = Y - dy/2

        # now we need to convert dip dirs to radians
        dip_dirs = np.array(dip_df['dip_azimuth'])
        dd_radians = np.radians(dip_dirs)
        # get the dx and dy of the dip dir line
        dx_dd = np.array([ (symbol_length/5)*np.sin(i) for i in dd_radians ])
        dy_dd = np.array([ (symbol_length/5)*np.cos(i) for i in dd_radians ])

        # get the dips as strings
        dips = list(dip_df['dip'])
        dips = [str(np.round(x,1)) for x in dips]

        bbox_props = dict(boxstyle="Round4,pad=0.1", fc="none", ec="none", lw=0.5,alpha = alpha)

        # now plot each symbol
        for i, angle in enumerate(strike_radians):
            # start and end for the strike
            end_X = start_X[i] + dx[i]
            end_Y = start_Y[i] + dy[i]
            # start and end for the dip dir
            end_X_dd = X[i] + dx_dd[i]
            end_Y_dd = Y[i] + dy_dd[i]
            # plot the lines
            self.ax_list[0].plot([start_X[i], end_X],[start_Y[i], end_Y],c=colour, lw=linewidth,alpha=alpha)
            self.ax_list[0].plot([X[i], end_X_dd], [Y[i], end_Y_dd], c=colour, lw=linewidth,alpha=alpha)
            # add the dip labelling
            self.ax_list[0].text(X[i], Y[i]-symbol_length, dips[i], fontsize=4, color=colour,alpha=alpha,bbox=bbox_props, ha= 'center',zorder=2)


    def plot_polygon_outlines(self,polygons, colour='black', linewidth=1, alpha = 1):
        """
        This function plots an outline of a series of shapely polygons
        Modified to also plot shapely Multipolygons if passed.

        Args:
            ax_list: list of axes
            polygons: dict of shapely polygons

        Author: FJC
        """
        from shapely.geometry import Polygon, MultiPolygon

        print('Plotting the polygon outlines...')

        for key, poly in polygons.items():
            if poly.geom_type == 'Polygon':
                x,y = poly.exterior.xy
                self.ax_list[0].plot(x,y, c=colour, lw = linewidth, alpha = alpha)
            elif poly.geom_type == 'MultiPolygon':
                for singlepoly in poly:
                    x,y = singlepoly.exterior.xy
                    self.ax_list[0].plot(x,y, c=colour, lw = linewidth, alpha = alpha)

    def plot_filled_polygons(self,polygons, facecolour='green', edgecolour='black', linewidth=1, alpha=0.5):
        """
        This function plots a series of shapely polygons but fills them in

        Args:
            ax_list: list of axes
            polygons: list of shapely polygons

        Author: FJC
        """
        from shapely.geometry import Polygon
        from descartes import PolygonPatch
        from matplotlib.collections import PatchCollection

        print('Plotting the polygons...')

        #patches = []
        for key, poly in polygons.items():
            this_patch = PolygonPatch(poly, fc=facecolour, ec=edgecolour, alpha=alpha)
            self.ax_list[0].add_patch(this_patch)

    def _set_coord_type(self, coord_type):
        """Sets the coordinate type

        Args:
            coord_type (string): The coordinate type. See options below

        Author: SMM
        """
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
        """
        Shows the plot. Duh.
        """
        self.fig.show()

    def save_fig(self,fig_width_inches = 4,FigFileName = 'TestFig.png',
                 FigFormat = 'png',Fig_dpi = 100,
                 axis_style = "Normal", transparent=False,
                 adjust_cbar_characters=True,
                 fixed_cbar_characters=4, return_fig = False):
        """
        This saves the figure to file.

        Args:
            fig_width_inches (float): The figure width in inches. It is in inches because yanks wrote matplotlib and that is the unit of measurement.
            FigFilename (str): The filename. You will need to add the directory here if you want it somewhere other than your python working directory.
            FigFormat (str): Figure format. Can be png, svg, pdf, gif, etc.
            Fig_dpi (int): dots per inch
            axis_stlye (string): This sets the axis style. Options are "Normal","Thick","Thin","Big", and "Madhouse"
            transparent (bool): If true the background is transparent (i.e., you don't get a white rectangle in the background)
            adjust_cbar_characters (bool): If true, adjust the spacing of the colourbar to account for the characters in the cbar label
            fixed_cbar_characters (int): ONLY used if adjust_cbar_characters=False. The number of characters to pad the cbar for.
            return_fig (bool): return the figure rather than saving a plot. In case you want some personnalisation. CAreful, if your personalisation may be useful for everyone, just code it for everyone.

        Author: SMM
        """

        self.ax_list = self.axis_styler(self.ax_list,axis_style)

        map_aspect_ratio = self._RasterList[0]._RasterAspectRatio
        print("The aspect ratio is: "+str(map_aspect_ratio))

        # We have to make some adjustments for the colourbar labels since if
        # they have lots of digits we need more space.
        # This is a fantastic thing you have to do when you are hard coding the
        # layout of the figure. Don't worry, this will all be worth it when we have awsome
        # figures every time.
        max_cbar_characters = 0
        if adjust_cbar_characters == True:
            print("I need to adjust the spacing of the colourbar.")
            if self.colourbar_location == "left" or self.colourbar_location == "right":
                print("You have a colourbar on the left or right, I need to check the number of characters in the labels.")
                labels = [item.get_text() for item in self.ax_list[-1].get_yticklabels()]
                print(labels)
                for label in labels:
                    if len(label) > max_cbar_characters:
                        max_cbar_characters = len(label)
                print("The longest colourbar label has "+str(max_cbar_characters)+" characters.")
        else:
            print("I am fixing the colourbar labels to have: "+str(fixed_cbar_characters)+ " characters")
            max_cbar_characters = fixed_cbar_characters



        fig = matplotlib.pyplot.gcf()

        # Now we size the figure. This requires some finessing since we are
        # hard coding the widths of everything

        if max_cbar_characters <= 3:
            cbar_width = 0.2
            cbar_text_width = 0.4   # This is in inches. Because yanks wrote matplotlib.
        else:
            cbar_width = 0.2
            cbar_text_width = 0.4+0.15*(max_cbar_characters-3)   # This is in inches. Because yanks wrote matplotlib.
            print("I'm adjusting the colourbar text width to "+str(cbar_text_width)+" inches")

        print("The cbar characters are: "+str(max_cbar_characters)+" and the cbar text width is: "+str(cbar_text_width))
        fig_size_inches, map_axes, cbar_axes = phelp.MapFigureSizer(fig_width_inches,
                                                              map_aspect_ratio,
                                                              self.colourbar_location,
                                                              self.title,
                                                              cbar_width,
                                                              cbar_text_width)

        fig.set_size_inches(fig_size_inches[0], fig_size_inches[1])
        self.ax_list[0].set_position(map_axes)

        # Annoying but the scatter plot resets the extents so you need to reassert them
        self.ax_list[0].set_xlim(self._xmin,self._xmax)
        #self.ax_list[0].set_ylim(self._ymax,self._ymin)
        self.ax_list[0].set_ylim(self._ymin,self._ymax)

        # add the title
        if self.title != "None":
            self.ax_list[0].set_title(self.title)

        print("Number of axes are: " + str(len(self.ax_list)))


        if cbar_axes == None:
            del self.ax_list[-1]
        else:
            self.ax_list[-1].set_position(cbar_axes)
        # if cbar_axes != None:
        #     self.ax_list[-1].set_position(cbar_axes)

        # I am returning the figure if wanted, otherwise I am saving the figure and clearing it
        if(return_fig):
            return fig
        else:
            # saving and closing
            fig.savefig(FigFileName, format=FigFormat, dpi=Fig_dpi, transparent=transparent)
            fig.clf()
            plt.close(fig)

    def SetRCParams(self,label_size):
        """
        This sets some RC params.

        Args:
            label_size(int): Font size of the labels

        Author: DAV

        """
        print("I am setting the font size to: "+str(label_size))
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['arial']
        rcParams['font.size'] = label_size
        rcParams['lines.linewidth']  = 1.5
