## LSDMap_Subplots.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## These functions are tools to deal with creating nice subplots from multiple
## rasters
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## FJC
## 22/12/2016
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from __future__ import absolute_import, division, print_function, unicode_literals

from glob import glob
import os.path
import numpy as np
import matplotlib.pyplot as pp
import string
import matplotlib.image as mpimg
import matplotlib.cm as cmx
from matplotlib import rcParams
from . import LSDMap_GDALIO as LSDMap_IO
from . import LSDMap_BasicPlotting as LSDMap_BP
from . import labels as lsdlabels

#==============================================================================
# Convert cm to inch for figure sizing
#------------------------------------------------------------------------------
def cm2inch(value):
    """
    Convert cm to inch for figure sizing
    """
    return value/2.54

#==============================================================================
# Function to create nice field sites figure from all the hillshades in a folder
# N_HSFiles = number of field sites
# Also reads in a series of map files to show locations of each site. At the moment
# you have to manually indicate on these maps where the field site is - would
# be nice to automate this when I have some time to mess around with it.
# FJC 22/12/16
#------------------------------------------------------------------------------

def field_sites(DataDirectory, N_HSFiles, NRows, NCols, n_target_ticks):

    # Set up fonts
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = 6

    # Read in the files for each site
    HSFiles = sorted(glob(DataDirectory+'*_HS.bil'), key=str)
    MapFiles = sorted(glob(DataDirectory+'*_map.eps'), key=str)
    print(MapFiles)

    n_files = len(HSFiles)+len(MapFiles)
    print("Number of files = ", n_files)

    # Now make the subplots
    fig, ax = pp.subplots(NRows,NCols, figsize=(cm2inch(12),cm2inch(15)), frameon=False)
    ax = ax.ravel()

    #get a list to label the subfigures
    alphabet = list(string.ascii_lowercase)
    file_counter = 0

    for i in range (n_files):

        if i < N_HSFiles:
        # first deal with the hillshade files. If there are more hillshade files simply
        # change the value of N_HSFiles
            print("The hillshade file name is: ", HSFiles[i])

            hillshade_raster = LSDMap_IO.ReadRasterArrayBlocks(HSFiles[i])
            hillshade_raster = np.ma.masked_where(hillshade_raster == -9999, hillshade_raster)

            # now get the extent
            extent_raster = LSDMap_IO.GetRasterExtent(HSFiles[i])

            # get DEM info
            CellSize,XMin,XMax,YMin,YMax = LSDMap_IO.GetUTMMaxMin(HSFiles[i])
            print(YMin, YMax)

            #plot the rasters
            ax[i].imshow(hillshade_raster, extent = extent_raster, cmap=cmx.gray)
            ax[i].text(0.05,0.97, alphabet[i], horizontalalignment='left', verticalalignment='top', bbox=dict(facecolor='white', edgecolor='k', pad=3), fontsize = 8, transform=ax[i].transAxes)

            #change ticks
            # xlocs = ax[i].xaxis.get_ticklocs()
            # ylocs = ax[i].yaxis.get_ticklocs()

            # new_xlocs,new_ylocs,new_x_labels,new_y_labels = LSDMap_BP.GetTicksForUTM(HSFiles[i],xlocs.max(),xlocs.min(),ylocs.max(),ylocs.min(),n_target_ticks)
            #
            # # change the location of the ticks depending on subplot placement
            # if i < 2:
            #     ax[i].xaxis.tick_top()
            # if i % 2 != 0:
            #     ax[i].yaxis.tick_right()

            ax[i].set_xticklabels([])
            ax[i].set_yticklabels([])
            #ax[i].tick_params(axis='x', pad=7)
            #ax[i].tick_params(axis='y', pad=7)

        if i >= N_HSFiles:
            print("The map file name is: ", MapFiles[file_counter])
            # add in the location maps for the sites (e.g. USA, UK)
            img = mpimg.imread(MapFiles[file_counter])
            ax[i].imshow(img)
            ax[i].text(0.02,1.00, alphabet[i], horizontalalignment='left', verticalalignment='top', fontsize = 8, transform=ax[i].transAxes)
            file_counter=file_counter+1

            #remove border
            ax[i].axis('off')

    # Save figure

    # Add a big subplot to get a common x and y label for the subplots
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    pp.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

    #pp.xlabel('Easting (m)', fontsize=8, labelpad=-122)
    #pp.ylabel('Northing (m)', fontsize=8, labelpad=10, position=(0.0,0.67))
    pp.tight_layout(pad=0.1, w_pad = 0.1, h_pad = 0.2)
    OutputFigureName = "field_sites"
    OutputFigureFormat = 'pdf'
    pp.savefig(DataDirectory+OutputFigureName + '.' + OutputFigureFormat, format=OutputFigureFormat, dpi=300)
    #pp.show()

#==============================================================================
# Function to create comparison plots for floodplain mapping between the published
# flood maps and the geometric method
# FJC 22/12/16
#------------------------------------------------------------------------------

def multiple_flood_maps(DataDirectory):
    """
    Make nice subplots of floodplain rasters for different field sites

    """
    import seaborn as sns

    # Set up fonts
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = 12

    # Read in the files for each site
    FPFiles = sorted(glob(DataDirectory+'*_FP*.bil'), key=str)

    n_files = len(FPFiles)
    print("Number of files = ", n_files)

    # Now make the subplots
    fig, ax = pp.subplots(2,3, figsize=(cm2inch(15),cm2inch(11)))
    ax = ax.ravel()

    #use seaborn to get a nice color palette
    cmap_oranges = sns.light_palette("#ff8f66", input="hex", as_cmap=True, reverse=True)
    cmap_ice = sns.light_palette("#00ffff", input="hex", as_cmap=True, reverse=True)

    alphabet = list(string.ascii_lowercase)

    for i in range (n_files):
        print("The floodplain file name is: ", FPFiles[i])

        # get the name of the field site
        fname = FPFiles[i].split('/')
        print(fname)
        split_fname = fname[-1].split('_')
        print(split_fname)
        HSFile = DataDirectory+split_fname[1]+'_'+split_fname[2]+"_HS.bil"
        print(HSFile)

        hillshade_raster = LSDMap_IO.ReadRasterArrayBlocks(HSFile)
        FP_raster = LSDMap_IO.ReadRasterArrayBlocks(FPFiles[i])
        FP_raster = np.ma.masked_where(FP_raster <= 0, FP_raster)

        # now get the extent
        extent_raster = LSDMap_IO.GetRasterExtent(HSFile)

        # get DEM info
        CellSize,XMin,XMax,YMin,YMax = LSDMap_IO.GetUTMMaxMin(HSFile)
        print(YMin, YMax)

        #plot the rasters
        ax[i].imshow(hillshade_raster, extent = extent_raster, cmap=cmx.gray)

        if i < 3:
            ax[i].imshow(FP_raster, extent = extent_raster, cmap=cmap_oranges, alpha=0.8)
        else:
            ax[i].imshow(FP_raster, extent = extent_raster, cmap=cmap_ice, alpha=0.6)

        ax[i].text(0.03,0.97, alphabet[i], bbox=dict(facecolor='white', edgecolor='k', pad=5), horizontalalignment='left', verticalalignment='top', transform=ax[i].transAxes)
        #scalebars.add_scalebar(ax[i], matchx=False, sizex=500.0, loc=3, borderpad =1, lw=3, matchy=False, hidex=False, hidey=False)

        ax[i].set_xticklabels([])
        ax[i].set_yticklabels([])
        ax[i].set_xticks([])
        ax[i].set_yticks([])

        ax[i].tick_params(axis='x', pad=3)
        ax[i].tick_params(axis='y', pad=3)

    # Save figure
    pp.tight_layout(pad=0.5, h_pad=0, w_pad=0.1)
    pp.subplots_adjust(wspace=0.05,hspace=0)
    OutputFigureName = "Comparison_published_maps"
    OutputFigureFormat = 'pdf'
    pp.savefig(DataDirectory+OutputFigureName + '.' + OutputFigureFormat, format=OutputFigureFormat, transparent=True, dpi=300)

#==============================================================================
#    Make subplots showing the difference between the mapped and predicted
#    floodplain initiation points. Uses Fiona (yaaay) to read in the shapefiles
#    FJC 05/01/17
#------------------------------------------------------------------------------

def flood_maps_with_shapefile(DataDirectory):

    from fiona import collection
    from descartes import PolygonPatch

    # Set up fonts
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = 8

    # Read in the files for each site
    HSFiles = sorted(glob(DataDirectory+'*_HS*.bil'), key=str)
    FPFiles = sorted(glob(DataDirectory+'*_FIPs_FP*.shp'), key=str)
    PointFiles = sorted(glob(DataDirectory+'*_FIPs_MP*.shp'), key=str)

    n_files = len(FPFiles)
    print("Number of files = ", n_files)

    # Now make the subplots
    fig, ax = pp.subplots(1,2, figsize=(cm2inch(12),cm2inch(7)))
    ax = ax.ravel()

    #get a list with the figure letterings
    figure_letter = ["a", "b"]
    titles = ["Mid Bailey Run, OH", "Coweeta, NC"]

    for i in range (n_files):

        print("The hillshade file name is: ", HSFiles[i])
        print("The floodplain file name is: ", FPFiles[i])
        print("The shapefile name is: ", PointFiles[i])

        # get the name of the field site
        split_fname = FPFiles[i].split('_FP')
        print(split_fname[0])

        hillshade_raster = LSDMap_IO.ReadRasterArrayBlocks(HSFiles[i])

        # now get the extent
        extent_raster = LSDMap_IO.GetRasterExtent(HSFiles[i])

        # get DEM info
        CellSize,XMin,XMax,YMin,YMax = LSDMap_IO.GetUTMMaxMin(HSFiles[i])
        print(YMin, YMax)

        # plot the raster
        ax[i].imshow(hillshade_raster, extent = extent_raster, cmap=cmx.gray)

        # plot the floodplain shapefile using fiona and descartes
        with collection(FPFiles[i], 'r') as input:
            for f in input:
                ax[i].add_patch(PolygonPatch(f['geometry'], fc='blue', ec='blue', lw=0.1, alpha=0.8))

        #plot the mapped points
        with collection(PointFiles[i],'r') as input:
            for point in input:
                x = point['geometry']['coordinates'][0]
                y = point['geometry']['coordinates'][1]
                ax[i].scatter(x,y, c="red", s=15, zorder=100)

        #ax[i].text(0.03,0.97, figure_letter[i], bbox=dict(facecolor='white', edgecolor='k', pad=3), horizontalalignment='left', verticalalignment='top', fontsize = 8, transform=ax[i].transAxes)

        #change ticks
        xlocs = ax[i].xaxis.get_ticklocs()
        ylocs = ax[i].yaxis.get_ticklocs()

        n_target_tics = 10
        new_xlocs,new_ylocs,new_x_labels,new_y_labels = LSDMap_BP.GetTicksForUTM(HSFiles[i],xlocs.max(),xlocs.min(),ylocs.max(),ylocs.min(),n_target_tics)

        ax[i].set_xticklabels(new_x_labels, rotation=30)
        ax[i].set_yticklabels(new_y_labels)

        #set axis limits
        ax[i].set_xlim(XMin,XMax)
        ax[i].set_ylim(YMin,YMax)

        if i == 0:
            ax[i].set_xlabel('Easting (m)', position=(1,0))
            ax[i].set_ylabel('Northing (m)')
        if i > 0:
            ax[i].yaxis.tick_right()

        ax[i].tick_params(axis='x', pad=2)
        ax[i].tick_params(axis='y', pad=2)

        ax[i].set_title(titles[i], fontsize=9)

    # Save figure
    pp.tight_layout(pad=0.5)
    OutputFigureName = "Comparison_with_mapped_FIPs"
    OutputFigureFormat = 'pdf'
    pp.savefig(DataDirectory+OutputFigureName + '.' + OutputFigureFormat, format=OutputFigureFormat, dpi=300)


def findmaxval_multirasters(FileList):
    """
    Loops through a list or array of rasters (np arrays)
    and finds the maximum single value in the set of arrays.
    """
    overall_max_val = 0

    for i in range (len(FileList)):

        raster_as_array = LSDMap_IO.ReadRasterArrayBlocks(FileList[i])
        this_max_val = np.nanmax(raster_as_array)

        if this_max_val > overall_max_val:
            overall_max_val = this_max_val
            print(overall_max_val)

    return overall_max_val

def findminval_multirasters(FileList):
    """
    Loops through a list or array of rasters (np arrays)
    and finds the minimum single value in the set of arrays.
    """
    overall_min_val = 0

    for i in range (len(FileList)):

        raster_as_array = LSDMap_IO.ReadRasterArrayBlocks(FileList[i])
        this_min_val = np.nanmin(raster_as_array)

        if this_min_val > overall_min_val:
            overall_min_val = this_min_val
            print(overall_min_val)

    return overall_min_val


def MultiDrapeFloodMaps(DataDir, ElevationRaster, DrapeRasterWild, cmap,
                        drape_min_threshold=None, drape_max=None, cbar_label=None):
    """Creates a figure with multiple drape maps over a hillshade.

    Plots flood extents from water depth rasters
    draped over the catchment elevation raster
    in a series of subplots

    Takes a wildcard for the drapes
    Expexts a fixed elevation raster, but this could be
    modified in future.

    Parameters:
	DataDir (str): Path to the directory containing the data files
	ElevationRaster (str): Name of the elevation raster used to create the hillshade
	DrapeRasterWild (str): Wildcard string used to find all the drape files in the
				directory.
	cmap: Can be the string name of a colourmap, or a Colourmap object
	drape_min (float, optional): Minimum value for the drape raster, i.e. values
		below this threshold will be masked and not plotted.
	drape_max (float, optional): Maximum value for the drape raster, i.e. values
		above this value will be masked and not plotted.
	cbar_label (str, optional): Label for the colourbar on the figure. This
		is the colourbar for the drape colourmap.

    Notes:
        Consider, if plotting multiple datasets, how you
        are going to deal with min a max values in the colur range.
        imshow will automatically set vmin and vmax and stretch the colour bar
        over this - which can be visually misleading. Ideally, you
        want to have the same colour map used for *all* subplots, and
        this is not default behaviour.

    Note: If `drape_max` is not set, the function searches for the maximum value
	in the range of rasters found by expanding the `DrapeRasterWild` argument
	and searching for the maximum value out of all rasters found.

    Raises:
	Exception: If the maximum value in the drape maps could not be found.

    """


    f, ax_arr = pp.subplots(2, 2, figsize=(10, 5), sharex=True, sharey=True)
    ax_arr = ax_arr.ravel()

    FPFiles = sorted(glob(DataDir+DrapeRasterWild), key=str)
    n_files = len(FPFiles)
    print("Number of files = ", n_files)

    elev_raster_file = DataDir + ElevationRaster

    hillshade = LSDMap_BP.Hillshade(elev_raster_file)
    #hillshade_array = LSDP.ReadRasterArrayBlocks(elev_raster_file)

    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(elev_raster_file)

    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # now get the tick marks
    n_target_tics = 5
    xlocs,ylocs,new_x_labels,new_y_labels = LSDMap_BP.GetTicksForUTM(elev_raster_file,x_max,x_min,y_max,y_min,n_target_tics)

    print("xmax: " + str(x_max))
    print("xmin: " + str(x_min))
    print("ymax: " + str(y_max))
    print("ymin: " + str(y_min))

    """
    Find the maximum water depth in all rasters.
    You need this to normalize the colourscale accross
    all plots when teh imshow is done later.
    """

    try:
        print("Calculating max drape raster value by scanning rasters...")
        max_water_depth = findmaxval_multirasters(FPFiles)
        drape_max = max_water_depth

    except:
        print("Something went wrong trying to obtain the max value in \
                your drape raster file list.")
    finally:
        print("The drape(s) max value is set to: ", drape_max)


    #im = mpimg.AxesImage()

    for i in range(n_files):

        print("The floodplain file name is: ", FPFiles[i])
        FP_raster = LSDMap_IO.ReadRasterArrayBlocks(FPFiles[i])
        #FP_raster = np.ma.masked_where(FP_raster <= 0, FP_raster)

        filename = os.path.basename(FPFiles[i])
        title = lsdlabels.make_line_label(filename)
        print(title)

        low_values_index = FP_raster < drape_min_threshold
        FP_raster[low_values_index] = np.nan

        im = ax_arr[i].imshow(hillshade, "gray", extent=extent_raster, interpolation="nearest")
        """
        Now we can set vmax to be the maximum water depth we calcualted earlier, making our separate
        subplots all have the same colourscale
        """
        im = ax_arr[i].imshow(FP_raster, cmap, extent=extent_raster,
                                alpha=1.0, interpolation="nearest",
                                vmin=drape_min_threshold,
                                vmax=drape_max)
        ax_arr[i].set_title(title)
        pp.setp( ax_arr[i].xaxis.get_majorticklabels(), rotation=70 )

    f.subplots_adjust(right=0.85)
    cax = f.add_axes([0.9, 0.1, 0.03, 0.8])

    cbar = f.colorbar(im, cax=cax)
    cbar.set_label(cbar_label)
    #cbar.set_ticks(np.linspace(0, 8, 8))
    #cbar = LSDP.colours.colorbar_index(f, cax, 8, cmap,
    #                                     drape_min_threshold, drape_max)

    #tick_locator = ticker.MaxNLocator(nbins=8)
    #cbar.locator = tick_locator
    #cbar.update_ticks()

    f.text(0.5, 0.04, 'Easting (m)', ha='center', fontsize=17)
    f.text(0.04, 0.5, 'Northing (m)', va='center', rotation='vertical', fontsize=17)


def MultiDrapeErodeDiffMaps(DataDir, ElevationRaster, DrapeRasterWild, cmap,
                        drape_min_threshold=None, cbar_label=None,
                        drape_max_threshold=None,
                        middle_mask_range=None):
    """Plots multiple drape maps of erosion/deposition (a DEM of difference)
       over a hillshade raster of the basin.

    Takes a wildcard for the drapes
    Expexts a single elevation raster for the background hillshade,
    but this could be modified in future.

    Parameters:
	DataDir (str): Path to the directory containing the data files
	ElevationRaster (str): Name of the elevation raster used to create the hillshade
	DrapeRasterWild (str): Wildcard string used to find all the drape files in the
				directory.
	cmap: Can be the string name of a colourmap, or a Colourmap object
	drape_min_threshold (float, optional): Minimum value for the drape raster, i.e. values
		below this threshold will be masked and not plotted.
	drape_max_threshold (float, optional): Maximum value for the drape raster, i.e. values
		above this value will be masked and not plotted.
	cbar_label (str, optional): Label for the colourbar on the figure. This
		is the colourbar for the drape colourmap.
     middle_mask_range (tuple, optional): A tuple or list of two values, used
                       to mask an inner range of values in the drape raster.
                       e.g. if you pass `(-0.1, 0.1)` then all the values
                       in the range -0.1 to 0.1 will be masked and not plotted
                       on the final map. Use for masking very small values
                       either side of zero.

    Notes:
        Consider, if plotting multiple datasets, how you
        are going to deal with min a max values in the colur range.
        imshow will automatically set vmin and vmax and stretch the colour bar
        over this - which can be visually misleading. Ideally, you
        want to have the same colour map used for *all* subplots, and
        this is not default behaviour.

    Note: If `drape_max_threshold` is not set, the function searches for the maximum value
	in the range of rasters found by expanding the `DrapeRasterWild` argument
	and searching for the maximum value out of all rasters found.

    Raises:
	Exception: If the maximum value in the drape maps could not be found.

    Author: DAV & FJC
    """
    #import lsdmatplotlibextensions as mplext

    f, ax_arr = pp.subplots(2, 2, figsize=(10, 5), sharex=True, sharey=True)
    ax_arr = ax_arr.ravel()

    FPFiles = sorted(glob(DataDir+DrapeRasterWild), key=str)
    n_files = len(FPFiles)
    print("Number of files = ", n_files)

    elev_raster_file = DataDir + ElevationRaster

    hillshade = LSDMap_BP.Hillshade(elev_raster_file)
    #hillshade_array = LSDP.ReadRasterArrayBlocks(elev_raster_file)

    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(elev_raster_file)

    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # now get the tick marks
    n_target_tics = 5
    xlocs, ylocs, new_x_labels, new_y_labels = LSDMap_BP.GetTicksForUTM(
                                                            elev_raster_file,
                                                            x_max,
                                                            x_min,
                                                            y_max,
                                                            y_min,
                                                            n_target_tics)

    print("xmax: " + str(x_max))
    print("xmin: " + str(x_min))
    print("ymax: " + str(y_max))
    print("ymin: " + str(y_min))

    """
    Find the maximum water depth in all rasters.
    You need this to normalize the colourscale accross
    all plots when teh imshow is done later.
    """
    if drape_max_threshold is None:
        try:
            print("Calculating max drape raster value by scanning rasters...")
            max_water_depth = findmaxval_multirasters(FPFiles)
            drape_max_threshold = max_water_depth

        except ValueError:
            print("Something went wrong trying to obtain the max value in \
                    your drape raster file list.")
        finally:
            print("The drape(s) max value is set to: ", drape_max_threshold)

    if drape_min_threshold is None:
        try:
            print("Calculating min drape raster value by scanning rasters...")
            min_water_depth = findminval_multirasters(FPFiles)
            drape_min_threshold = min_water_depth

        except ValueError:
            print("Something went wrong trying to obtain the min value in \
                    your drape raster file list.")
        finally:
            print("The drape(s) min value is set to: ", drape_min_threshold)


    for i in range(n_files):

        print("The floodplain file name is: ", FPFiles[i])
        FP_raster = LSDMap_IO.ReadRasterArrayBlocks(FPFiles[i])
        #FP_raster = np.ma.masked_where(FP_raster <= 0, FP_raster)

        filename = os.path.basename(FPFiles[i])
        title = lsdlabels.make_line_label(filename)
        print(title)

        # Mask the extreme high values
        hi_values_index = FP_raster > drape_max_threshold
        FP_raster[hi_values_index] = np.nan

        # Mask the extreme low values
        lo_values_index = FP_raster < drape_min_threshold
        FP_raster[lo_values_index] = np.nan

        # Mask the middle values that are really close to zero (i.e. if you
        # have negative and positive values in the raster, such as in a DEM
        # of difference with both erosion and deposition.)
        if middle_mask_range is not None:
            masked_mid_values_index = (np.logical_and(FP_raster > middle_mask_range[0],
                                                      FP_raster < middle_mask_range[1]))
            FP_raster[masked_mid_values_index] = np.nan

        im = ax_arr[i].imshow(hillshade, "gray", extent=extent_raster, interpolation="nearest")
        """
        Now we can set vmax to be the maximum water depth we calcualted earlier, making our separate
        subplots all have the same colourscale
        """
        im = ax_arr[i].imshow(FP_raster, cmap, extent=extent_raster,
                                alpha=1.0, interpolation="nearest",
                                vmin=drape_min_threshold,
                                vmax=drape_max_threshold)
        ax_arr[i].set_title(title)
        pp.setp( ax_arr[i].xaxis.get_majorticklabels(), rotation=70 )

    f.subplots_adjust(right=0.85)
    cax = f.add_axes([0.9, 0.1, 0.03, 0.8])

    cbar = f.colorbar(im, cax=cax)
    cbar.set_label(cbar_label)
    #cbar.set_ticks(np.linspace(0, 8, 8))
    #cbar = mplext.colours.colorbar_index(f, cax, 8, cmap,
    #                                     drape_min_threshold, drape_max)

    #tick_locator = ticker.MaxNLocator(nbins=8)
    #cbar.locator = tick_locator
    #cbar.update_ticks()

    f.text(0.5, 0.04, 'Easting (m)', ha='center', fontsize=17)
    f.text(0.04, 0.5, 'Northing (m)', va='center', rotation='vertical', fontsize=17)
