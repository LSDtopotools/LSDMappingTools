## LSDMap_ChiPlotting.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## These functions are tools to deal with chi maps
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## SMM
## 14/12/2016
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import os
from . import cubehelix
import matplotlib.pyplot as plt
import matplotlib as mpl
#from cycler import cycler
from matplotlib import rcParams
from matplotlib import colors
import LSDPlottingTools.LSDMap_GDALIO as LSDMap_IO
#import LSDMap_BasicManipulation as LSDMap_BM
#import LSDMap_OSystemTools as LSDOst
import LSDPlottingTools.LSDMap_BasicPlotting as LSDMap_BP
import LSDPlottingTools.LSDMap_PointTools as LSDMap_PD
import LSDPlottingTools.LSDMap_BasicManipulation as LSDMap_BM
import LSDPlottingTools.statsutilities as LSDStats
from LSDMapFigure.PlottingRaster import MapFigure
from LSDMapFigure.PlottingRaster import BaseRaster
from LSDMapFigure import PlottingHelpers as Helper
from LSDPlottingTools import colours as lsdcolours
from LSDPlottingTools import statsutilities as SUT
from LSDPlottingTools import init_plotting_DV
from LSDPlottingTools import adjust_text
import LSDPlottingTools as LSDP
import LSDMapFigure.PlottingHelpers as Helper


def ConvertBasinIndexToJunction(BasinPointData,BasinIndexList):
    """This transforms a basin index list (simply the order of the basins, starting from low to high junction number) to a junction list.

    This allows users to go between basin rasters (with junctions listed) and the simpler basin indexing system (which is sequential)

    Args:
        BasinPointData (LSDMap_PointData): a point data object
        BasinIndexList (list of ints): The basin indices to be converted to junctions

    Returns:
        A list on ints with the basin junctions

    Author: SMM
    """



    these_data = BasinPointData.QueryData("outlet_junction")
    these_data = [int(x) for x in these_data]

    #print("The junctions are: ")




    basin_junction_list = []
    for basinindex in BasinIndexList:
        basin_junction_list.append(these_data[basinindex])


    return basin_junction_list


##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def FindSourceInformation(thisPointData):
    """This function finds the source locations, with chi elevation, flow distance, etc.

    Args:
        thisPointData (LSDMap_PointData) A LSDMap_PointData object that is derived from the Chi_mapping_tool component of *LSDTopoTools*.

    Returns:
        A dict with key of the source node that returns a dict that has the FlowDistance, Chi, and Elevation of each source.
        Used for plotting source numbers on profile plots.

    Author: SMM
    """

    # Get the chi, m_chi, basin number, and source ID code
    chi = thisPointData.QueryData('chi')
    chi = [float(x) for x in chi]
    elevation = thisPointData.QueryData('elevation')
    elevation = [float(x) for x in elevation]
    fdist = thisPointData.QueryData('flow distance')
    fdist = [float(x) for x in fdist]
    source = thisPointData.QueryData('source_key')
    source = [int(x) for x in source]
    latitude = thisPointData.GetLatitude()
    longitude = thisPointData.GetLongitude()




    Chi = np.asarray(chi)
    Elevation = np.asarray(elevation)
    Fdist = np.asarray(fdist)
    Source = np.asarray(source)
    Latitude = np.asarray(latitude)
    Longitude = np.asarray(longitude)

    n_sources = Source.max()+1
    print("N sources is: "+str(n_sources))

    # This loops through all the source indices, and then picks out the
    # Elevation, chi coordinate and flow distance of each node
    # Then it returns a dictionary containing the elements of the node
    these_source_nodes = {}
    for src_idx in range(0,n_sources):
        m = np.ma.masked_where(Source!=src_idx, Source)

        # Mask the unwanted values
        maskX = np.ma.masked_where(np.ma.getmask(m), Chi)
        maskElevation = np.ma.masked_where(np.ma.getmask(m), Elevation)
        maskFlowDistance = np.ma.masked_where(np.ma.getmask(m), Fdist)
        maskLatitude = np.ma.masked_where(np.ma.getmask(m), Latitude)
        maskLongitude= np.ma.masked_where(np.ma.getmask(m), Longitude)

        # get the locations of the source
        this_dict = {}
        idx_of_max_FD = maskX.argmax()
        this_dict["FlowDistance"]=maskFlowDistance[idx_of_max_FD]
        this_dict["Chi"]=maskX[idx_of_max_FD]
        this_dict["Elevation"]=maskElevation[idx_of_max_FD]
        this_dict["Latitude"]=maskLatitude[idx_of_max_FD]
        this_dict["Longitude"]=maskLongitude[idx_of_max_FD]

        # get the minimum of the source
        idx_of_min_Chi = maskX.argmin()
        chi_length = maskX[idx_of_max_FD]-maskX[idx_of_min_Chi]
        this_dict["SourceLength"]=chi_length

        these_source_nodes[src_idx] = this_dict

    return these_source_nodes

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def FindShortSourceChannels(these_source_nodes,threshold_length):
    """This function gets the list of sources that are shorter than a threshold value

    Args:
        these_source_nodes (dict): A dict from the FindSourceInformation module
        threshold_length (float): The threshold of chi lenght of the source segment

    Return:
        long_sources: A list of integers of source with the appropriate length

    Author: SMM
    """
    long_sources = []
    for key in these_source_nodes:
        if these_source_nodes[key]["SourceLength"] > threshold_length:
            long_sources.append(key)

    return long_sources


##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## This function plots the chi slope on a shaded relief map
## It uses the Kirby and Whipple colour scheme
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def BasicChiPlotGridPlotKirby(FileName, DrapeName, chi_csv_fname, thiscmap='gray',drape_cmap='gray',
                            colorbarlabel='Elevation in meters',clim_val = (0,0),
                            drape_alpha = 0.6,FigFileName = 'Image.pdf',FigFormat = 'show',
                            elevation_threshold = 0, size_format = "ESURF", dpi_save = 500):

    """This function plots the chi slope on a shaded relief map. It uses the Kirby and Whipple colour scheme.

    Args:
        FileName (str): The name (with full path and extension) of the DEM.
        DrapenName (str): The name (with full path and extension) of the drape file (usually a hillshade, but could be anything)
        chi_csv_fname (str): The name (with full path and extension) of the cdv file with chi, chi slope, etc information. This file is produced by the chi_mapping_tool.
        thiscmap (colormap): The colourmap for the elevation raster
        drape_cmap (colormap):  The colourmap for the drape raster
        colorbarlabel (str): the text label on the colourbar.
        clim_val  (float,float): The colour limits for the drape file. If (0,0) it uses the minimum and maximum values of the drape file. Users can assign numbers to get consistent colourmaps between plots.
        drape_alpha (float): The alpha value of the drape
        FigFileName (str): The name of the figure file
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        elevation_threshold (float): elevation_threshold chi points below this elevation are removed from plotting.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).

    Returns:
        Does not return anything but makes a plot.

    Author: SMM
    """

    from matplotlib import colors

    label_size = 10

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size
    #plt.rc('text', usetex=True)

    # get the data
    raster = LSDMap_IO.ReadRasterArrayBlocks(FileName)
    raster_drape = LSDMap_IO.ReadRasterArrayBlocks(DrapeName)

    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(FileName)

    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # make a figure,
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    gs = plt.GridSpec(100,100,bottom=0.25,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])

    gs = plt.GridSpec(100,100,bottom=0.25,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])

    # This is the axis for the colorbar
    ax2 = fig.add_subplot(gs[10:15,15:70])

    # now get the tick marks
    n_target_tics = 5
    xlocs,ylocs,new_x_labels,new_y_labels = LSDMap_BP.GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics)

    im1 = ax.imshow(raster[::-1], thiscmap, extent = extent_raster, interpolation="nearest")

    # set the colour limits
    print("Setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1]))
    if (clim_val == (0,0)):
        print("Im setting colour limits based on minimum and maximum values")
        im1.set_clim(0, np.nanmax(raster))
    else:
        print("Now setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1]))
        im1.set_clim(clim_val[0],clim_val[1])

    plt.hold(True)

    # Now for the drape: it is in grayscale
    #print "drape_cmap is: "+drape_cmap
    im3 = ax.imshow(raster_drape[::-1], drape_cmap, extent = extent_raster, alpha = drape_alpha, interpolation="nearest")

    # Set the colour limits of the drape
    im3.set_clim(0,np.nanmax(raster_drape))

    ax.spines['top'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)

    ax.set_xticklabels(new_x_labels,rotation=60)
    ax.set_yticklabels(new_y_labels)

    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
    ax.tick_params(axis='both', width=1, pad = 2)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(2)

    # Now we get the chi points
    EPSG_string = LSDMap_IO.GetUTMEPSG(FileName)
    print("EPSG string is: " + EPSG_string)

    thisPointData = LSDMap_PD.LSDMap_PointData(chi_csv_fname)
    thisPointData.ThinData('elevation',elevation_threshold)

    # convert to easting and northing
    [easting,northing] = thisPointData.GetUTMEastingNorthing(EPSG_string)


    # The image is inverted so we have to invert the northing coordinate
    Ncoord = np.asarray(northing)
    Ncoord = np.subtract(extent_raster[3],Ncoord)
    Ncoord = np.add(Ncoord,extent_raster[2])

    M_chi = thisPointData.QueryData('m_chi')
    #print M_chi
    M_chi = [float(x) for x in M_chi]


    # make a color map of fixed colors
    this_cmap = colors.ListedColormap(['#2c7bb6','#abd9e9','#ffffbf','#fdae61','#d7191c'])
    bounds=[0,50,100,175,250,1205]
    norm = colors.BoundaryNorm(bounds, this_cmap.N)

    sc = ax.scatter(easting,Ncoord,s=0.5, c=M_chi,cmap=this_cmap,norm=norm,edgecolors='none')

    # This affects all axes because we set share_all = True.
    ax.set_xlim(x_min,x_max)
    ax.set_ylim(y_max,y_min)

    ax.set_xticks(xlocs)
    ax.set_yticks(ylocs)

    cbar = plt.colorbar(sc,cmap=this_cmap,norm=norm,spacing='uniform', ticks=bounds, boundaries=bounds,orientation='horizontal',cax=ax2)
    cbar.set_label(colorbarlabel, fontsize=10)
    ax2.set_xlabel(colorbarlabel, fontname='Arial',labelpad=l_pad)

    print("The figure format is: " + FigFormat)
    if FigFormat == 'show':
        plt.show()
    elif FigFormat == 'return':
        return fig
    else:
        plt.savefig(FigFileName,format=FigFormat,dpi=dpi_save)
        fig.clf()




##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## This function plots the chi slope on a shaded relief map
## It uses a cubehelix colourmap over the log 10 of the channel steepness
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def BasicChiPlotGridPlot(FileName, DrapeName, chi_csv_fname, thisPointData, thiscmap='gray',drape_cmap='gray',
                            colorbarlabel='log$_{10}k_{sn}$',clim_val = (0,0),
                            drape_alpha = 0.6,FigFileName = 'Image.pdf',FigFormat = 'show',
                            elevation_threshold = 0, size_format = "ESURF", dpi_save = 750):

    """This is the main chi plotting script that prints a chi steepness map over the hillshade. Note that the colour scale for the chi slope values are always cubehelix

    Args:
        FileName (str): The name (with full path and extension) of the DEM.
        DrapenName (str): The name (with full path and extension) of the drape file (usually a hillshade, but could be anything)
        chi_csv_fname (str): The name (with full path and extension) of the cdv file with chi, chi slope, etc information. This file is produced by the chi_mapping_tool.
        thiscmap (colormap): The colourmap for the elevation raster
        drape_cmap (colormap):  The colourmap for the drape raster
        colorbarlabel (str): the text label on the colourbar.
        clim_val  (float,float): The colour limits for the drape file. If (0,0) it uses the minimum and maximum values of the drape file. Users can assign numbers to get consistent colourmaps between plots.
        drape_alpha (float): The alpha value of the drape
        FigFileName (str): The name of the figure file
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        elevation_threshold (float): elevation_threshold chi points below this elevation are removed from plotting.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).

    Returns:
        Prints a plot to file.

    Author:
        Simon M. Mudd

    """

    import math
    label_size = 10

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # get the data
    raster = LSDMap_IO.ReadRasterArrayBlocks(FileName)
    raster_drape = LSDMap_IO.ReadRasterArrayBlocks(DrapeName)

    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(FileName)
    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # make a figure,
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    gs = plt.GridSpec(100,100,bottom=0.25,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])

    # now get the tick marks
    n_target_tics = 5
    xlocs,ylocs,new_x_labels,new_y_labels = LSDMap_BP.GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics)

    im1 = ax.imshow(raster[::-1], thiscmap, extent = extent_raster, interpolation="nearest")

    # set the colour limits
    print("Setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1]))
    if (clim_val == (0,0)):
        print("Im setting colour limits based on minimum and maximum values")
        im1.set_clim(0, np.nanmax(raster))
    else:
        print("Now setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1]))
        im1.set_clim(clim_val[0],clim_val[1])

    plt.hold(True)

    # Now for the drape: it is in grayscale
    #print "drape_cmap is: "+drape_cmap
    im3 = ax.imshow(raster_drape[::-1], drape_cmap, extent = extent_raster, alpha = drape_alpha, interpolation="nearest")

    # Set the colour limits of the drape
    im3.set_clim(0,np.nanmax(raster_drape))

    # Set up axes and ticks
    ax.spines['top'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)
    ax.set_xticklabels(new_x_labels,rotation=60)
    ax.set_yticklabels(new_y_labels)
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
    ax.tick_params(axis='both', width=1, pad = 2)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(2)

    # Now we get the chi points
    EPSG_string = LSDMap_IO.GetUTMEPSG(FileName)
    print("EPSG string is: " + EPSG_string)

    #thisPointData = LSDMap_PD.LSDMap_PointData(chi_csv_fname)
    thisPointData.ThinData('elevation',elevation_threshold)

    # convert to easting and northing
    [easting,northing] = thisPointData.GetUTMEastingNorthing(EPSG_string)

    # The image is inverted so we have to invert the northing coordinate
    Ncoord = np.asarray(northing)
    Ncoord = np.subtract(extent_raster[3],Ncoord)
    Ncoord = np.add(Ncoord,extent_raster[2])

    M_chi = thisPointData.QueryData('m_chi')
    M_chi = [float(x) for x in M_chi]


    log_m_chi = []
    for value in M_chi:
        if value < 0.1:
            log_m_chi.append(0)
        else:
            log_m_chi.append(math.log10(value))
    colorbarlabel = "log$_{10}k_{sn}$"

    this_cmap = cubehelix.cmap(rot=1, reverse=True,start=3,gamma=1.0,sat=2.0)
    sc = ax.scatter(easting,Ncoord,s=0.5, c=log_m_chi,cmap=this_cmap,edgecolors='none')

    # set the colour limits
    sc.set_clim(0, np.nanmax(log_m_chi))


    # This is the axis for the colorbar
    ax2 = fig.add_subplot(gs[10:15,15:70])
    plt.colorbar(sc,cmap=this_cmap,spacing='uniform', orientation='horizontal',cax=ax2)
    ax2.set_xlabel(colorbarlabel, fontname='Arial',labelpad=l_pad)

    # This affects all axes because we set share_all = True.
    ax.set_xlim(x_min,x_max)
    ax.set_ylim(y_max,y_min)

    ax.set_xticks(xlocs)
    ax.set_yticks(ylocs)

    print("The figure format is: " + FigFormat)
    if FigFormat == 'show':
        plt.show()
    elif FigFormat == 'return':
        return fig
    else:
        plt.savefig(FigFileName,format=FigFormat,dpi=dpi_save)
        fig.clf()



##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## This function plots the chi slope on a shaded relief map
## It uses a cubehelix colourmap over the log 10 of the channel steepness
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def BasicChiCoordinatePlot(FileName, DrapeName, csvfile, thiscmap='gray',drape_cmap='gray',
                            colorbarlabel='$\chi (m)$',clim_val = (0,0),
                            basin_order_list = [], basin_point_data = "None", basin_raster_name = "None",
                            drape_alpha = 0.6,FigFileName = 'Image.pdf',FigFormat = 'show',
                            size_format = "ESURF"):

    """This plots the chi coordinate, mimicking Sean Willet et al's plots

    Args:
        FileName (str): The name (with full path and extension) of the DEM.
        DrapenName (str): The name (with full path and extension) of the drape file (usually a hillshade, but could be anything)
        thisPointData (LSDMap_PointData): The point data object with the basic chi points
        thiscmap (colormap): The colourmap for the elevation raster
        drape_cmap (colormap):  The colourmap for the drape raster
        colorbarlabel (str): the text label on the colourbar.
        clim_val  (float,float): The colour limits for the drape file. If (0,0) it uses the minimum and maximum values of the drape file. Users can assign numbers to get consistent colourmaps between plots.
        basin_order_list (list of int): The basin indices to be selected
        basin_point_data (LSDM_PointData): The mapping between junctions and indices
        basin_raster_name (str): If a basin raster name is supplied the chi raster will be masked
        drape_alpha (float): The alpha value of the drape
        FigFileName (str): The name of the figure file
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).

    Returns:
        Prints a plot to file.

    Author:
        Simon M. Mudd

    """

    label_size = 10

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # get the data
    raster = LSDMap_IO.ReadRasterArrayBlocks(FileName)
    raster_drape = LSDMap_IO.ReadRasterArrayBlocks(DrapeName)

    raster_drape = LSDMap_BM.NanBelowThreshold(raster_drape,0)

    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(FileName)
    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # make a figure,
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    gs = plt.GridSpec(100,100,bottom=0.25,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])

    # now get the tick marks
    n_target_tics = 5
    xlocs,ylocs,new_x_labels,new_y_labels = LSDMap_BP.GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics)

    # Lets do the ticks in km
    n_hacked_digits = 3
    new_x_labels = LSDMap_BP.TickLabelShortenizer(new_x_labels,n_hacked_digits)
    new_y_labels = LSDMap_BP.TickLabelShortenizer(new_y_labels,n_hacked_digits)

    print("This cmap for base raster is: "+thiscmap)
    im1 = ax.imshow(raster[::-1], thiscmap, extent = extent_raster, interpolation="nearest")

    # set the colour limits
    print("Setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1]))
    if (clim_val == (0,0)):
        print("Im setting colour limits based on minimum and maximum values")
        im1.set_clim(0, np.nanmax(raster))
    else:
        print("Now setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1]))
        im1.set_clim(clim_val[0],clim_val[1])

    plt.hold(True)

    # Now for the drape: it is in grayscale
    #print "drape_cmap is: "+drape_cmap
    im3 = ax.imshow(raster_drape[::-1], drape_cmap, extent = extent_raster, alpha = drape_alpha, interpolation="nearest")

    # Set the colour limits of the drape
    im3.set_clim(0,np.nanmax(raster_drape))

    # Set up axes and ticks
    ax.spines['top'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)
    ax.set_xticklabels(new_x_labels,rotation=60)
    ax.set_yticklabels(new_y_labels)
    ax.set_xlabel("Easting (km)")
    ax.set_ylabel("Northing (km)")

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
    ax.tick_params(axis='both', width=1, pad = 2)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(2)

    # Now we get the chi points
    EPSG_string = LSDMap_IO.GetUTMEPSG(FileName)
    print("EPSG string is: " + EPSG_string)

    print("The file is: "+ csvfile)
    thisPointData = LSDMap_PD.LSDMap_PointData(csvfile)


    # Check if there are basins, and if we are to mask them
    if not len(basin_order_list) == 0:
        # check if there is a point data object
        if basin_point_data != "None":
            basin_junction_list = ConvertBasinIndexToJunction(basin_point_data,basin_order_list)
            print("I am thinning to the following basins: ")
            print(basin_junction_list)
            thisPointData.ThinDataSelection('basin_junction',basin_junction_list)

            if basin_raster_name != "None":
                basin_raster = LSDMap_IO.ReadRasterArrayBlocks(basin_raster_name)
                LSDMap_BM.MaskByCategory(raster_drape,basin_raster,basin_junction_list)


    # convert to easting and northing
    [easting,northing] = thisPointData.GetUTMEastingNorthing(EPSG_string)

    # The image is inverted so we have to invert the northing coordinate
    Ncoord = np.asarray(northing)
    Ncoord = np.subtract(extent_raster[3],Ncoord)
    Ncoord = np.add(Ncoord,extent_raster[2])

    chi = thisPointData.QueryData('chi')
    chi = [float(x) for x in chi]


    #this_cmap = 'brg_r'
    this_cmap = 'CMRmap_r'
    sc = ax.scatter(easting,Ncoord,s=0.5, c=chi,cmap=this_cmap,edgecolors='none')

    # set the colour limits
    sc.set_clim(0, np.nanmax(chi))


    # This is the axis for the colorbar
    ax2 = fig.add_subplot(gs[10:15,15:70])
    plt.colorbar(sc,cmap=this_cmap,spacing='uniform', orientation='horizontal',cax=ax2)
    ax2.set_xlabel(colorbarlabel, fontname='Arial',labelpad=l_pad)

    # This affects all axes because we set share_all = True.
    ax.set_xlim(x_min,x_max)
    ax.set_ylim(y_max,y_min)

    ax.set_xticks(xlocs)
    ax.set_yticks(ylocs)

    print("The figure format is: " + FigFormat)
    if FigFormat == 'show':
        plt.show()
    elif FigFormat == 'return':
        return fig
    else:
        plt.savefig(FigFileName,format=FigFormat,dpi=750)
        fig.clf()



##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## This function plots channels, color coded
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def BasicChannelPlotGridPlotCategories(FileName, DrapeName, chi_csv_fname, thiscmap='gray',drape_cmap='gray',
                            colorbarlabel='Elevation in meters',clim_val = (0,0),
                            drape_alpha = 0.6,FigFileName = 'Image.pdf',FigFormat = 'show',
                            elevation_threshold = 0, data_name = 'source_key',
                            source_thinning_threshold = 0,
                            size_format = "ESURF"):
    """This plots the channels over a draped plot, colour coded by source

    Args:
        FileName (str): The name (with full path and extension) of the DEM.
        DrapenName (str): The name (with full path and extension) of the drape file (usually a hillshade, but could be anything)
        chi_csv_fname (str): The name (with full path and extension) of the cdv file with chi, chi slope, etc information. This file is produced by the chi_mapping_tool.
        thiscmap (colormap): The colourmap for the elevation raster
        drape_cmap (colormap):  The colourmap for the drape raster
        colorbarlabel (str): the text label on the colourbar.
        clim_val  (float,float): The colour limits for the drape file. If (0,0) it uses the minimum and maximum values of the drape file. Users can assign numbers to get consistent colourmaps between plots.
        drape_alpha (float): The alpha value of the drape
        FigFileName (str): The name of the figure file
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        elevation_threshold (float): elevation_threshold chi points below this elevation are removed from plotting.
        data_name (str) = The name of the sources csv
        source_thinning_threshold (float) = Minimum chi length of a source segment. No thinning if 0.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).

    Returns:
        Prints a plot to file.

    Author:
        Simon M. Mudd

    """
    from matplotlib import colors
    label_size = 10

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size
    #plt.rc('text', usetex=True)

    # get the data
    raster = LSDMap_IO.ReadRasterArrayBlocks(FileName)
    raster_drape = LSDMap_IO.ReadRasterArrayBlocks(DrapeName)

    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(FileName)

    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # make a figure,
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,5))
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))

    gs = plt.GridSpec(100,100,bottom=0.25,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])

    # now get the tick marks
    n_target_tics = 5
    xlocs,ylocs,new_x_labels,new_y_labels = LSDMap_BP.GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics)

    im1 = ax.imshow(raster[::-1], thiscmap, extent = extent_raster, interpolation="nearest")

    # set the colour limits
    print("Setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1]))
    if (clim_val == (0,0)):
        print("Im setting colour limits based on minimum and maximum values")
        im1.set_clim(0, np.nanmax(raster))
    else:
        print("Now setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1]))
        im1.set_clim(clim_val[0],clim_val[1])

    plt.hold(True)

    # Now for the drape: it is in grayscale
    #print "drape_cmap is: "+drape_cmap
    im3 = ax.imshow(raster_drape[::-1], drape_cmap, extent = extent_raster, alpha = drape_alpha, interpolation="nearest")

    # Set the colour limits of the drape
    im3.set_clim(0,np.nanmax(raster_drape))


    ax.spines['top'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)

    ax.set_xticklabels(new_x_labels,rotation=60)
    ax.set_yticklabels(new_y_labels)

    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
    ax.tick_params(axis='both', width=1, pad = 2)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(2)

    # Now we get the chi points
    EPSG_string = LSDMap_IO.GetUTMEPSG(FileName)
    print("EPSG string is: " + EPSG_string)

    thisPointData = LSDMap_PD.LSDMap_PointData(chi_csv_fname)
    thisPointData.ThinData('elevation',elevation_threshold)

    # Logic for thinning the sources
    if source_thinning_threshold > 0:
        print("I am going to thin some sources out for you")
        source_info = FindSourceInformation(thisPointData)
        remaining_sources = FindShortSourceChannels(source_info,source_thinning_threshold)
        thisPointData.ThinDataSelection("source_key",remaining_sources)
    else:
        print("I am not thinning by source")

    # convert to easting and northing
    print("Converting to easting and northing")
    [easting,northing] = thisPointData.GetUTMEastingNorthing(EPSG_string)

    # The image is inverted so we have to invert the northing coordinate
    Ncoord = np.asarray(northing)
    Ncoord = np.subtract(extent_raster[3],Ncoord)
    Ncoord = np.add(Ncoord,extent_raster[2])

    these_data = thisPointData.QueryData(data_name)
    #print M_chi
    these_data = [int(x) for x in these_data]

    # make a color map of fixed colors
    NUM_COLORS = 15

    this_cmap = plt.cm.Set1
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    plt.cm.ScalarMappable(norm=cNorm, cmap=this_cmap)
    channel_data = [x % NUM_COLORS for x in these_data]

    print("Printing the scatter map")
    ax.scatter(easting,Ncoord,s=0.5, c=channel_data,norm=cNorm,cmap=this_cmap,edgecolors='none')

    # This affects all axes because we set share_all = True.
    ax.set_xlim(x_min,x_max)
    ax.set_ylim(y_max,y_min)

    ax.set_xticks(xlocs)
    ax.set_yticks(ylocs)
    ax.set_title('Channels colored by source number')

    print("The figure format is: " + FigFormat)
    if FigFormat == 'show':
        plt.show()
    elif FigFormat == 'return':
        return fig
    else:
        plt.savefig(FigFileName,format=FigFormat,dpi=500)
        fig.clf()

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## This function plots channels, color coded
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def BasicChannelPlotByBasin(FileName, DrapeName, chi_csv_fname, thiscmap='gray',drape_cmap='gray',
                            colorbarlabel='Elevation in meters',clim_val = (0,0),
                            drape_alpha = 0.6,FigFileName = 'Image.pdf',FigFormat = 'show',
                            elevation_threshold = 0, basin_key = 0, data_name = 'source_key',
                            source_thinning_threshold = 0,
                            size_format = "ESURF"):
    """This plots the channels over a draped plot, colour coded by source. It masks the data so that the channels
    are only plotted for a specific basin of interest, specified by the basin key from the chi csv file.

    Args:
        FileName (str): The name (with full path and extension) of the DEM.
        DrapenName (str): The name (with full path and extension) of the drape file (usually a hillshade, but could be anything)
        chi_csv_fname (str): The name (with full path and extension) of the cdv file with chi, chi slope, etc information. This file is produced by the chi_mapping_tool.
        thiscmap (colormap): The colourmap for the elevation raster
        drape_cmap (colormap):  The colourmap for the drape raster
        colorbarlabel (str): the text label on the colourbar.
        clim_val  (float,float): The colour limits for the drape file. If (0,0) it uses the minimum and maximum values of the drape file. Users can assign numbers to get consistent colourmaps between plots.
        drape_alpha (float): The alpha value of the drape
        FigFileName (str): The name of the figure file
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command. If 'return' then it returns the figure.
        elevation_threshold (float): elevation_threshold chi points below this elevation are removed from plotting.
        basin_key (int): the ID of the basin from the chi csv file
        data_name (str) = The name of the sources csv
        source_thinning_threshold (float) = Minimum chi length of a source segment. No thinning if 0.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).

    Returns:
        Prints a plot to file.

    Author:
        FJC (modified from SMM)

    """
    from matplotlib import colors
    label_size = 10

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size
    #plt.rc('text', usetex=True)

    # get the data
    raster = LSDMap_IO.ReadRasterArrayBlocks(FileName)
    raster_drape = LSDMap_IO.ReadRasterArrayBlocks(DrapeName)

    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(FileName)

    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # make a figure,
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,5))
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))

    gs = plt.GridSpec(100,100,bottom=0.25,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])

    # now get the tick marks
    n_target_tics = 5
    xlocs,ylocs,new_x_labels,new_y_labels = LSDMap_BP.GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics)

    im1 = ax.imshow(raster[::-1], thiscmap, extent = extent_raster, interpolation="nearest")

    # set the colour limits
    print("Setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1]))
    if (clim_val == (0,0)):
        print("Im setting colour limits based on minimum and maximum values")
        im1.set_clim(0, np.nanmax(raster))
    else:
        print("Now setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1]))
        im1.set_clim(clim_val[0],clim_val[1])

    plt.hold(True)

    # Now for the drape: it is in grayscale
    #print "drape_cmap is: "+drape_cmap
    im3 = ax.imshow(raster_drape[::-1], drape_cmap, extent = extent_raster, alpha = drape_alpha, interpolation="nearest")

    # Set the colour limits of the drape
    im3.set_clim(0,np.nanmax(raster_drape))


    ax.spines['top'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)

    ax.set_xticklabels(new_x_labels,rotation=60)
    ax.set_yticklabels(new_y_labels)

    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
    ax.tick_params(axis='both', width=1, pad = 2)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(2)

    # Now we get the chi points
    EPSG_string = LSDMap_IO.GetUTMEPSG(FileName)
    print("EPSG string is: " + EPSG_string)

    thisPointData = LSDMap_PD.LSDMap_PointData(chi_csv_fname)

    # mask for the basin
    thisPointData.ThinDataSelection("basin_key",basin_key)

    # thin elevation points below the threshold
    thisPointData.ThinData('elevation',elevation_threshold)

    # Logic for thinning the sources
    if source_thinning_threshold > 0:
        print("I am going to thin some sources out for you")
        source_info = FindSourceInformation(thisPointData)
        remaining_sources = FindShortSourceChannels(source_info,source_thinning_threshold)
        thisPointData.ThinDataSelection("source_key",remaining_sources)

    # convert to easting and northing
    [easting,northing] = thisPointData.GetUTMEastingNorthing(EPSG_string)

    # The image is inverted so we have to invert the northing coordinate
    Ncoord = np.asarray(northing)
    Ncoord = np.subtract(extent_raster[3],Ncoord)
    Ncoord = np.add(Ncoord,extent_raster[2])

    these_data = thisPointData.QueryData(data_name)
    #print M_chi
    these_data = [int(x) for x in these_data]

    # make a color map of fixed colors
    NUM_COLORS = 15

    this_cmap = plt.cm.Set1
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    plt.cm.ScalarMappable(norm=cNorm, cmap=this_cmap)
    channel_data = [x % NUM_COLORS for x in these_data]

    ax.scatter(easting,Ncoord,s=0.5, c=channel_data,norm=cNorm,cmap=this_cmap,edgecolors='none')

    # This affects all axes because we set share_all = True.
    ax.set_xlim(x_min,x_max)
    ax.set_ylim(y_max,y_min)

    ax.set_xticks(xlocs)
    ax.set_yticks(ylocs)
    ax.set_title('Channels colored by source number')

    print("The figure format is: " + FigFormat)
    if FigFormat == 'show':
        plt.show()
    elif FigFormat == 'return':
        return fig
    else:
        plt.savefig(FigFileName,format=FigFormat,dpi=500)
        fig.clf()

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## This function plots channels, color coded
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def ChiProfiles(chi_csv_fname, FigFileName = 'Image.pdf',FigFormat = 'show',
                basin_order_list = [],basin_rename_list = [],
                label_sources = False,
                elevation_threshold = 0,
                source_thinning_threshold = 0, plot_M_chi = False,
                size_format = "ESURF",
                plot_segments = False):
    """This function plots the chi vs elevation: lumps everything onto the same axis. This tends to make a mess.

    Args:
        chi_csv_fname (str): The name (with full path and extension) of the cdv file with chi, chi slope, etc information. This file is produced by the chi_mapping_tool.
        FigFileName (str): The name of the figure file
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        basin_order_list (int list): The basins to plot
        basin_rename_list (int list): A list for naming substitutions
        label_sources (bool): If true, label the sources.
        elevation_threshold (float): elevation_threshold chi points below this elevation are removed from plotting.
        source_thinning_threshold (float) = Minimum chi length of a source segment. No thinning if 0
        plot_MChi (bool): If true, plots chi against MChi
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).

    Returns:
         Does not return anything but makes a plot.

    Author: SMM
    """

    from matplotlib import colors
    from .adjust_text import adjust_text

    label_size = 10

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    if plot_M_chi:
        print("I will plot chi vs M_chi")
    else:
        print("I will plot ci vs elevation")

    # make a figure,
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,5))
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))

    gs = plt.GridSpec(100,100,bottom=0.25,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])

    thisPointData = LSDMap_PD.LSDMap_PointData(chi_csv_fname)
    thisPointData.ThinData('elevation',elevation_threshold)

    # Logic for thinning the sources
    if source_thinning_threshold > 0:
        print("I am going to thin some sources out for you")
        source_info = FindSourceInformation(thisPointData)
        remaining_sources = FindShortSourceChannels(source_info,source_thinning_threshold)
        thisPointData.ThinDataSelection("source_key",remaining_sources)

    # Logic for stacked labels. You need to run this after source thinning to
    # get an updated source dict
    if label_sources:
        source_info = FindSourceInformation(thisPointData)



    # Get the chi, m_chi, basin number, and source ID code
    chi = thisPointData.QueryData('chi')
    chi = [float(x) for x in chi]
    elevation = thisPointData.QueryData('elevation')
    elevation = [float(x) for x in elevation]
    fdist = thisPointData.QueryData('flow distance')
    fdist = [float(x) for x in fdist]
    m_chi = thisPointData.QueryData('m_chi')
    m_chi = [float(x) for x in m_chi]
    basin = thisPointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    source = thisPointData.QueryData('source_key')
    source = [int(x) for x in source]

    segments = thisPointData.QueryData('segment_number')
    segments = [int(x) for x in segments]
    segmented_elevation = thisPointData.QueryData('segmented_elevation')
    segmented_elevation = [float(x) for x in segmented_elevation]

    # Some booleans that tell if there are segments and segmented elevation
    have_segments = False
    if len(segments) == len(chi):
        have_segments = True
        print("I've got the segments")
    have_segmented_elevation = False
    if len(segmented_elevation) == len(chi):
        have_segmented_elevation = True
        print("I've got segmented elevation")
    else:
        print("I don't have the segmented elevation")

    print("The number of data points are: " +str(len(chi)))

    # need to convert everything into arrays so we can mask different basins
    Chi = np.asarray(chi)
    Elevation = np.asarray(elevation)
    #Fdist = np.asarray(fdist)
    M_chi = np.asarray(m_chi)
    Basin = np.asarray(basin)
    Source = np.asarray(source)

    Segments = np.asarray(segments)
    Segmented_elevation = np.asarray(segmented_elevation)

    #max_basin = np.amax(Basin)
    max_chi = np.amax(Chi)
    max_Elevation = np.amax(Elevation)
    max_M_chi = np.amax(M_chi)
    min_Elevation = np.amin(Elevation)

    if plot_M_chi:
        z_axis_min = 0
        z_axis_max = int(max_M_chi/10)*10+10
    else:
        z_axis_min = int(min_Elevation/10)*10
        z_axis_max = int(max_Elevation/10)*10+10

    chi_axis_max = int(max_chi/5)*5+5

    # make a color map of fixed colors
    NUM_COLORS = 2

    this_cmap = plt.cm.Set1
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    Basin_colors = [x % NUM_COLORS for x in Basin]


    dot_pos = FigFileName.rindex('.')
    newFilename = FigFileName[:dot_pos]+FigFileName[dot_pos:]
    print("newFilename: "+newFilename)

    if len(basin_order_list) == 0:
        print("No basins in this list!")
        print("I am defaulting to look at basin 0")
        basin_order_list.append(0)

    texts = []
    bbox_props = dict(boxstyle="circle,pad=0.1", fc="w", ec="k", lw=0.5,alpha = 0.25)
    for basin_number in basin_order_list:

        print(("This basin is: " +str(basin_number)))

        m = np.ma.masked_where(Basin!=basin_number, Basin)
        maskX = np.ma.masked_where(np.ma.getmask(m), Chi)
        if plot_M_chi:
            maskElevation = np.ma.masked_where(np.ma.getmask(m), M_chi)
        else:
            maskElevation = np.ma.masked_where(np.ma.getmask(m), Elevation)

        maskBasin = np.ma.masked_where(np.ma.getmask(m), Basin_colors)
        maskSource = np.ma.masked_where(np.ma.getmask(m), Source)

        if(have_segmented_elevation):
            # We need to loop through the sources.
            # We can do that by converting the sources to a set and then looping through the set
            sources_list = maskSource.tolist()
            myset = set(sources_list)
            print("The sources are: ")
            print(myset)
            sources_list = list(myset)
            for source in sources_list:
                m_mask = np.ma.masked_where(Source!=source, Source)
                mask_maskX = np.ma.masked_where(np.ma.getmask(m_mask), Chi)
                maskSegmentedElevation = np.ma.masked_where(np.ma.getmask(m_mask), Segmented_elevation)
                a_line, = ax.plot(mask_maskX,maskSegmentedElevation,'b',alpha = 0.6)
                a_line.set_dashes([3,1])


        # logic for source labeling
        if label_sources:

            # Convert the masked data to a list and then that list to a set and
            # back to a list (phew!)
            list_source = maskSource.tolist()
            set_source = set(list_source)
            list_source = list(set_source)

            # Now we have to get rid of stupid non values
            list_source = [x for x in list_source if x is not None]

            print("these sources are: ")
            print(list_source)

            for this_source in list_source:
                source_Chi= source_info[this_source]["Chi"]

                if plot_M_chi:
                    source_Elevation = source_info[this_source]["M_chi"]
                else:
                    source_Elevation = source_info[this_source]["Elevation"]
                #print("Source is: "+str(this_source))
                #print("Chi is: "+str(source_info[this_source]["Chi"]))
                #print("FlowDistance is is: "+str(source_info[this_source]["FlowDistance"]))
                #print("Elevation is: "+str(source_info[this_source]["Elevation"]))
                texts.append(ax.text(source_Chi, source_Elevation, str(this_source), style='italic',
                        verticalalignment='bottom', horizontalalignment='left',fontsize=8,bbox=bbox_props))


        ax.scatter(maskX,maskElevation,s=2.0, c=maskBasin,norm=cNorm,cmap=this_cmap,edgecolors='none')



    ax.spines['top'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)

    ax.set_xlabel("$\chi$")
    ax.set_ylabel("Elevation (m)")

    # This affects all axes because we set share_all = True.
    ax.set_ylim(z_axis_min,z_axis_max)
    ax.set_xlim(0,chi_axis_max)


    #print("Number of text elements is: "+str(len(texts)))
    adjust_text(texts)
    #print("Now the number of text elements is: "+str(len(texts)))


    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
    ax.tick_params(axis='both', width=1, pad = 2)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(2)

    print("The figure format is: " + FigFormat)
    if FigFormat == 'show':
        plt.show()
    elif FigFormat == 'return':
        return fig
    else:
        plt.savefig(newFilename,format=FigFormat,dpi=500)
        fig.clf()

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## This function plots channels, color coded
## Only plot the source colouring, not the chi gradient!!
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def StackedChiProfiles(chi_csv_fname, FigFileName = 'Image.pdf',
                       FigFormat = 'show',elevation_threshold = 0,
                       first_basin = 0, last_basin = 0,
                       basin_order_list = [],basin_rename_list = [],
                       X_offset = 5,label_sources = False,
                       source_thinning_threshold = 0,
                       size_format = "ESURF"):
    """This function plots the chi vs elevation: It stacks profiles (so the basins are spaced out) and colours them by the source number.

    Args:
        chi_csv_fname (str): The name (with full path and extension) of the cdv file with chi, chi slope, etc information. This file is produced by the chi_mapping_tool.
        FigFileName (str): The name of the figure file
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        elevation_threshold (float): elevation_threshold chi points below this elevation are removed from plotting.
        first_basin (int): The basin to start with (but overridden by the basin list)
        last_basin (int): The basin to end with (but overridden by the basin list)
        basin_order_list (int list): The basins to plot
        basin_rename_list (int list): A list for naming substitutions. Useful because LSDTopoTools might number basins in a way a human wouldn't, so a user can intervene in the names.
        X_offset (float): The offest in chi between the basins along the x-axis. Used to space out the profiles so you can see each of them.
        label_sources (bool): If true, label the sources.
        source_thinning_threshold (float) = Minimum chi length of a source segment. No thinning if 0.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).

    Returns:
         Does not return anything but makes a plot.

    Author: SMM
    """

    from .adjust_text import adjust_text
    from matplotlib import colors
    import matplotlib.patches as patches

    label_size = 10

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size


    # make a figure,
    if size_format == "geomorphology":
        print("I am plotting for geomorphology")
        fig = plt.figure(1, facecolor='white',figsize=(6.25,4))
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))

    gs = plt.GridSpec(100,100,bottom=0.25,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])

    thisPointData = LSDMap_PD.LSDMap_PointData(chi_csv_fname)
    thisPointData.ThinData('elevation',elevation_threshold)
    thisPointData.ThinData('chi',0)

    # Thin the sources.
    if source_thinning_threshold > 0:
        print("I am going to thin some sources out for you")
        source_info = FindSourceInformation(thisPointData)
        remaining_sources = FindShortSourceChannels(source_info,source_thinning_threshold)
        print("The remaining number of sources are: "+str(len(remaining_sources)))
        thisPointData.ThinDataSelection("source_key",remaining_sources)

    # Get the chi, m_chi, basin number, and source ID code
    chi = thisPointData.QueryData('chi')
    chi = [float(x) for x in chi]
    elevation = thisPointData.QueryData('elevation')
    elevation = [float(x) for x in elevation]
    fdist = thisPointData.QueryData('flow distance')
    fdist = [float(x) for x in fdist]
    m_chi = thisPointData.QueryData('m_chi')
    m_chi = [float(x) for x in m_chi]
    basin = thisPointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    source = thisPointData.QueryData('source_key')
    source = [int(x) for x in source]

    # need to convert everything into arrays so we can mask different basins
    Chi = np.asarray(chi)
    Elevation = np.asarray(elevation)
    #Fdist = np.asarray(fdist)
    #M_chi = np.asarray(m_chi)
    Basin = np.asarray(basin)
    Source = np.asarray(source)

    max_basin = np.amax(Basin)
    max_chi = np.amax(Chi)
    max_Elevation = np.amax(Elevation)
    #max_M_chi = np.amax(M_chi)
    min_Elevation = np.amin(Elevation)

    # determine the maximum and minimum elevations
    z_axis_min = int(min_Elevation/10)*10
    z_axis_max = int(max_Elevation/10)*10+10
    X_axis_max = int(max_chi/5)*5+5

    elevation_range = z_axis_max-z_axis_min
    z_axis_min = z_axis_min - 0.075*elevation_range

    # Now calculate the spacing of the stacks
    this_X_offset = 0
    if basin_order_list:
        basins_list = basin_order_list

        n_stacks = len(basins_list)
        added_X = X_offset*n_stacks
        X_axis_max = X_axis_max+added_X
    else:
        # now loop through a number of basins
        if last_basin >= max_basin:
            last_basin = max_basin-1

        if first_basin > last_basin:
            first_basin = last_basin
            print("Your first basin was larger than last basin. I won't plot anything")
        basins_list = list(range(first_basin,last_basin+1))

        n_stacks = last_basin-first_basin+1
        added_X = X_offset*n_stacks
        print(("The number of stacks is: "+str(n_stacks)+" the old max: "+str(X_axis_max)))
        X_axis_max = X_axis_max+added_X
        print(("The nex max is: "+str(X_axis_max)))


    # make a color map of fixed colors
    NUM_COLORS = 15

    this_cmap = plt.cm.Set1
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    #scalarMap = plt.cm.ScalarMappable(norm=cNorm, cmap=this_cmap)
    Source_colors = [x % NUM_COLORS for x in Source]
    plt.hold(True)

    # Logic for stacked labels. You need to run this after source thinning to
    # get an updated source dict
    if label_sources:
        source_info = FindSourceInformation(thisPointData)

    dot_pos = FigFileName.rindex('.')
    newFilename = FigFileName[:dot_pos]+'_Stack'+str(first_basin)+FigFileName[dot_pos:]

    texts = []
    # Format the bounding box of source labels
    bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="b", lw=0.5,alpha = 0.5)

    for basin_number in basins_list:

        print(("This basin is: " +str(basin_number)))

        m = np.ma.masked_where(Basin!=basin_number, Basin)
        maskX = np.ma.masked_where(np.ma.getmask(m), Chi)
        maskElevation = np.ma.masked_where(np.ma.getmask(m), Elevation)
        maskSource = np.ma.masked_where(np.ma.getmask(m), Source_colors)

        print(("adding an offset of: "+str(this_X_offset)))

        maskX = np.add(maskX,this_X_offset)

        this_min_x = np.nanmin(maskX)
        this_max_x =np.nanmax(maskX)
        width_box = this_max_x-this_min_x

        print(("Min: "+str(this_min_x)+" Max: "+str(this_max_x)))
        ax.add_patch(patches.Rectangle((this_min_x,z_axis_min), width_box, z_axis_max-z_axis_min,alpha = 0.01,facecolor='r',zorder=-10))

        if basin_number == basins_list[-1]:
            print(("last basin, geting maximum value,basin is: "+str(basin_number)))
            this_max = np.amax(maskX)
            this_max = int(this_max/5)*5+5
            print(("The rounded maximum is: "+str(this_max)))
            chi_axis_max = this_max

        #Source_colors = [x % NUM_COLORS for x in maskSource]

        # some logic for the basin rename
        if basin_rename_list:
            if len(basin_rename_list) == max_basin+1:
                this_basin_text = "Basin "+str(basin_rename_list[basin_number])
        else:
            this_basin_text = "Basin "+str(basin_number)


        ax.text(this_min_x+0.1*width_box, z_axis_min+0.025*elevation_range, this_basin_text, style='italic',
                verticalalignment='bottom', horizontalalignment='left',fontsize=8)

        # logic for source labeling
        if label_sources:

            # Convert the masked data to a list and then that list to a set and
            # back to a list (phew!)
            list_source = maskSource.tolist()
            set_source = set(list_source)
            list_source = list(set_source)

            # Now we have to get rid of stupid non values
            list_source = [x for x in list_source if x is not None]

            print("these sources are: ")
            print(list_source)

            for this_source in list_source:
                source_Chi= source_info[this_source]["Chi"]
                source_Elevation = source_info[this_source]["Elevation"]
                print(("Source is: "+str(this_source)))
                #print("Chi is: "+str(source_info[this_source]["Chi"]))
                #print("FlowDistance is is: "+str(source_info[this_source]["FlowDistance"]))
                #print("Elevation is: "+str(source_info[this_source]["Elevation"]))
                texts.append(ax.text(source_Chi+this_X_offset, source_Elevation, str(this_source), style='italic',
                        verticalalignment='bottom', horizontalalignment='left',fontsize=8,bbox=bbox_props))


        ax.scatter(maskX,maskElevation,s=2.0, c=maskSource,norm=cNorm,cmap=this_cmap,edgecolors='none')
        this_X_offset = this_X_offset+X_offset


    ax.spines['top'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)

    ax.set_xlabel("$\chi$ (m)")
    ax.set_ylabel("Elevation (m)")

    # This affects all axes because we set share_all = True.
    ax.set_ylim(z_axis_min,z_axis_max)
    ax.set_xlim(0,chi_axis_max)

    # adjust the text
    adjust_text(texts)

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
    ax.tick_params(axis='both', width=1, pad = 2)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(2)

    print("The figure format is: " + FigFormat)
    if FigFormat == 'show':
        plt.show()
    elif FigFormat == 'return':
        return fig
    else:
        plt.savefig(newFilename,format=FigFormat,dpi=500)
        fig.clf()

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## This function plots channels, color coded in chi space with a gradient
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def StackedProfilesGradient(chi_csv_fname, FigFileName = 'Image.pdf',
                       FigFormat = 'png',elevation_threshold = 0,
                       first_basin = 0, last_basin = 0, basin_order_list = [],
                       basin_rename_dict = {},
                       this_cmap = "viridis",axis_data_name = 'chi', colour_data_name = "m_chi", colorbarlabel = "Colourbar", cbar_loc = "bottom",
                       discrete_colours = False, NColours = 15, X_offset = 5,
                       plotting_data_format = 'log',
                       label_sources = False, source_thinning_threshold = 0,
                       size_format = "ESURF", aspect_ratio = 2, dpi = 500,
                       stack_patches = False, rotate_labels = False):
    """This function plots the chi vs elevation or flow distance vs elevation.

    It stacks profiles (so the basins are spaced out).
    It colours the plots by the chi steepness (which is equal to the normalised channel steepness if A_0 is set to 1).

    Args:
        chi_csv_fname (str): The name (with full path and extension) of the cdv file with chi, chi slope, etc information. This file is produced by the chi_mapping_tool.
        FigFileName (str): The name of the figure file
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        elevation_threshold (float): elevation_threshold chi points below this elevation are removed from plotting.
        first_basin (int): The basin to start with (but overridden by the basin list)
        last_basin (int): The basin to end with (but overridden by the basin list)
        basin_order_list (int list): The basins to plot
        basin_rename_list (int list): A list for naming substitutions. Useful because LSDTopoTools might number basins in a way a human wouldn't, so a user can intervene in the names.
        this_cmap (colormap): NOT USED! We now use a default colourmap but this may change.
        data_name (str): 'chi' or 'flow_distance' What to plot along the x-axis.
        X_offset (float): The offest in chi between the basins along the x-axis. Used to space out the profiles so you can see each of them.
        plotting_data_format: NOT USED previously if 'log' use logarithm scale, but we now automatically do this. Might change later.
        label_sources (bool): If true, label the sources.
        source_thinning_threshold (float) = Minimum chi length of a source segment. No thinning if 0.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        aspect_ratio (flt): the aspect ratio of the figure
        dpi (int): dots per inch of figure
        stack_patches (bool): if true, places a rectangular patch element behind each profile.
        rotate_labels (bool): if true, places the labels rotated at the top of each chi plot.

    Returns:
         Does not return anything but makes a plot.

    Author: SMM
    """

    import math
    import matplotlib.patches as patches
    #from adjust_text import adjust_text

    label_size = 10

    print("STARTING stacks. Cmap is: "+this_cmap+ " and the offset is: " + str(X_offset))
    print("The rename dict is: ")
    print(basin_rename_dict)

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure,
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        fig_size_inches = 6.25
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        fig_size_inches = 16
        l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))
        fig_size_inches = 4.92126
        l_pad = -35

    # Note all the below parameters are overwritten by the figure sizer routine
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])

    print("Getting data from the file: "+chi_csv_fname)
    thisPointData = LSDMap_PD.LSDMap_PointData(chi_csv_fname)

    print("I am going to thin your data to elevation and chi thresholds for you.")
    thisPointData.GetLongitude(PrintToScreen = True)
    thisPointData.selectValue('elevation',value = elevation_threshold,operator = ">")
    thisPointData.selectValue('chi',value = 0,operator = ">")



    #SK = thisPointData.QueryData("source_key")
    #print("Source keys are")
    #print(SK)

    # Thin the sources. Do this after the colouring so that thinned source colours
    # will be the same as unthinned source colours.
    if source_thinning_threshold > 0:
        print("I am going to thin some sources out for you")
        source_info = FindSourceInformation(thisPointData)
        remaining_sources = FindShortSourceChannels(source_info,source_thinning_threshold)
        print("The remaining number of sources are: "+str(len(remaining_sources)))
        print("The remaining sources are: ")
        print(remaining_sources)
        thisPointData.ThinDataSelection("source_key",remaining_sources)

    # Get the chi, m_chi, basin number, and source ID code
    if axis_data_name  == 'chi':
        x_data = thisPointData.QueryData('chi').values
        x_data = [float(x) for x in x_data]
    elif axis_data_name == 'flow_distance':
        x_data = thisPointData.QueryData('flow_distance').values
        x_data = [float(x) for x in x_data]
    else:
        print("I did not understand the data name. Choices are chi and flow distance. Defaulting to chi.")
        x_data = thisPointData.QueryData('chi')
        x_data = [float(x) for x in x_data]

    elevation = thisPointData.QueryData('elevation').values
    elevation = [float(x) for x in elevation]
    m_chi = thisPointData.QueryData(colour_data_name).values
    m_chi = [float(x) for x in m_chi]
    basin = thisPointData.QueryData('basin_key').values
    basin = [int(x) for x in basin]
    source = thisPointData.QueryData('source_key').values
    source = [int(x) for x in source]

    # make a color map of fixed colors
    if discrete_colours:

        this_cmap = plt.cm.Set1
        #cNorm  = colors.Normalize(vmin=0, vmax=NColours-1)
        m_chi = [x % NColours for x in m_chi]
    else:

        if (plotting_data_format == 'log'):
            log_m_chi = []
            for value in m_chi:
                if value < 0.1:
                    log_m_chi.append(0)
                else:
                    log_m_chi.append(math.log10(value))
            m_chi = log_m_chi


    # need to convert everything into arrays so we can mask different basins
    Xdata = np.asarray(x_data)
    Elevation = np.asarray(elevation)
    M_chi = np.asarray(m_chi)
    Basin = np.asarray(basin)
    Source = np.asarray(source)

    max_basin = np.amax(Basin)
    max_X = np.amax(Xdata)
    max_Elevation = np.amax(Elevation)
    max_M_chi = np.amax(M_chi)
    min_Elevation = np.amin(Elevation)
    min_M_chi = np.amin(M_chi)

    print(("Max M_chi is: "+str(max_M_chi)))

    if rotate_labels == False:
        z_axis_min = int(min_Elevation/10)*10
        z_axis_max = int(max_Elevation/10)*10+10
        X_axis_max = int(max_X/5)*5+5
        M_chi_axis_max = max_M_chi
        M_chi_axis_min = min_M_chi

        elevation_range = z_axis_max-z_axis_min
        z_axis_min = z_axis_min - 0.075*elevation_range
    else:
        z_axis_min = int(min_Elevation-500)
        z_axis_max = int(max_Elevation/10)*10+10
        X_axis_max = int(max_X/5)*5+5
        M_chi_axis_max = max_M_chi
        if min_M_chi < 0:
            M_chi_axis_min = 0
        else:
            M_chi_axis_min = min_M_chi

        elevation_range = z_axis_max-z_axis_min
        z_axis_min = z_axis_min - 0.075*elevation_range

    plt.hold(True)


    # Now calculate the spacing of the stacks
    this_X_offset = 0

    print("The chi offset is: "+ str(X_offset))
    if basin_order_list:
        print("You have supplied a basin order list so I am igoring the minimum and maximum basin number. ")
        basins_list = basin_order_list

        n_stacks = len(basins_list)
        added_X = X_offset*n_stacks
        X_axis_max = X_axis_max+added_X
    else:
        print("You did not supply a list of basins. I am using a first and last basin to plot.")
        # now loop through a number of basins
        if last_basin >= max_basin:
            last_basin = max_basin-1

        if first_basin > last_basin:
            first_basin = last_basin
            print("Your first basin was larger than last basin. I won't plot anything")

        # Make a list containing the basins
        basins_list = list(range(first_basin,last_basin+1))

        # This simply gets the maximum chi coordinate as the number of plots times the chi offset.
        # This isn't ideal since it could cut off the last stack
        n_stacks = last_basin-first_basin+1
        added_X = X_offset*n_stacks
        print(("The number of stacks is: "+str(n_stacks)+" the old max: "+str(X_axis_max)))
        X_axis_max = X_axis_max+added_X
        print(("The new max is: "+str(X_axis_max)))


    # Logic for stacked labels. You need to run this after source thinning to
    # get an updated source dict
    if label_sources:
        source_info = FindSourceInformation(thisPointData)

    # Now start looping through the basins
    #dot_pos = FigFileName.rindex('.')
    #newFilename = FigFileName[:dot_pos]+'_GradientStack'+str(first_basin)+FigFileName[dot_pos:]


    texts = []
    # Format the bounding box of source labels
    bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="b", lw=0.5,alpha = 0.5)

    print("I am going to loop through your basins.")
    for basin_number in basins_list:

        print(("This basin is: " +str(basin_number)))

        m = np.ma.masked_where(Basin!=basin_number, Basin)
        maskX = np.ma.masked_where(np.ma.getmask(m), Xdata)
        maskElevation = np.ma.masked_where(np.ma.getmask(m), Elevation)
        maskMChi = np.ma.masked_where(np.ma.getmask(m), M_chi)
        maskSource = np.ma.masked_where(np.ma.getmask(m), Source)

        print("adding an offset of: "+str(this_X_offset))

        # Get the minimum and maximum
        this_min_x = np.nanmin(maskX)
        if this_min_x < 0:
            this_min_x = 0
        this_max_x =np.nanmax(maskX)
        width_box = this_max_x-this_min_x

        # First subtract the minimum from the X data
        maskX = np.subtract(maskX,this_min_x)

        # Now add the offset to the minimum and maximum
        this_min_x = this_X_offset
        this_max_x = width_box

        # Now add the offset to the data
        maskX = np.add(maskX,this_X_offset)

        # This adds a rectangular coloured patch behind the stack. Helps to highlight which profiles are which.
        # Usually turned off since it looks a bit weird when they overlap
        if stack_patches:
            print("Min: "+str(this_min_x)+" Max: "+str(this_max_x))
            ax.add_patch(patches.Rectangle((this_min_x,z_axis_min), width_box, z_axis_max-z_axis_min,alpha = 0.01,facecolor='r',zorder=-10))

        # some logic for the basin rename. We use a dictionary for this
        if basin_rename_dict:
            if basin_number in basin_rename_dict:
                print("I found this basin in the rename list. I am renaming.")
                this_basin_text = str(basin_rename_dict[basin_number])
            else:
                this_basin_text = "Basin "+str(basin_number)
        else:
            this_basin_text = "Basin "+str(basin_number)

        # Add the basin text
        if rotate_labels == False:
            ax.text(this_min_x+0.1*width_box, z_axis_min+0.025*elevation_range, this_basin_text, style='italic',
                    verticalalignment='bottom', horizontalalignment='left',fontsize=8)
        else:
            ax.text(this_min_x+0.05*width_box, z_axis_min+0.055*elevation_range, this_basin_text, style='italic',
                    verticalalignment='bottom', horizontalalignment='left',fontsize=12, rotation=90)

        # Here is some logic for getting the maximum value for the last basin
        if basin_number == basins_list[-1]:
            print(("last basin, geting maximum value,basin is: "+str(basin_number)))
            value_max = np.amax(maskX)
            this_max = int(value_max/5)*5+5
            print(("The rounded maximum is: "+str(this_max)))
            if this_max-value_max <0.5:
                X_axis_max = this_max+1
            else:
                X_axis_max = this_max

        # logic for source labeling
        if label_sources:

            # Convert the masked data to a list and then that list to a set and
            # back to a list (phew!)
            list_source = maskSource.tolist()
            set_source = set(list_source)
            list_source = list(set_source)

            # Now we have to get rid of stupid non values
            list_source = [x for x in list_source if x is not None]

            #print("these sources are: ")
            #print list_source

            #print("the source info is: ")
            #print source_info

            for this_source in list_source:

                if data_name == 'chi':
                    source_X = source_info[this_source]["Chi"]
                elif data_name == 'flow_distance':
                    source_X = source_info[this_source]["FlowDistance"]
                else:
                    source_X = source_info[this_source]["Chi"]

                source_Elevation = source_info[this_source]["Elevation"]
                #print("Source is: "+str(this_source))
                #print("Chi is: "+str(source_info[this_source]["Chi"]))
                #print("FlowDistance is is: "+str(source_info[this_source]["FlowDistance"]))
                #print("Elevation is: "+str(source_info[this_source]["Elevation"]))
                texts.append(ax.text(source_X+this_X_offset, source_Elevation, str(this_source), style='italic',
                        verticalalignment='bottom', horizontalalignment='left',fontsize=8,bbox=bbox_props))

        # Now plot the scatter for this stack. The colour limits are for all plots
        cnorm = colors.Normalize( M_chi_axis_min, M_chi_axis_max)
        sc = ax.scatter(maskX,maskElevation,s=2.0, c=maskMChi,cmap=this_cmap,edgecolors='none',norm = cnorm, vmin = M_chi_axis_min, vmax = M_chi_axis_max)

        # increment the offset
        this_X_offset = this_X_offset+X_offset
        print("The new chi offset is: "+str(this_X_offset))


    # This is the axis for the colorbar
    if cbar_loc != "None":

        cbar_orient = "horizontal"
        if cbar_loc == "right" or cbar_loc == "left":
            cbar_orient = "vertical"

        ax2 = fig.add_axes([0.1,0.8,0.2,0.5])
        cbar = mpl.colorbar.ColorbarBase(ax2, cmap=this_cmap,
                                norm=cnorm,
                                orientation=cbar_orient)

        #Will's changes:
        # Changed rotation of colourbar text to 90 and the labelpad to -75 for "left"

        if cbar_loc == 'top':
            ax2.set_xlabel(colorbarlabel, fontname='Arial',labelpad=5)
        elif cbar_loc == 'bottom':
            ax2.set_xlabel(colorbarlabel, fontname='Arial',labelpad=5)
        elif cbar_loc == 'left':
            ax2.set_ylabel(colorbarlabel, fontname='Arial',labelpad=-75,rotation=90)
        elif cbar_loc == 'right':
            ax2.set_ylabel(colorbarlabel, fontname='Arial',labelpad=10,rotation=270)

        #ax2 = fig.add_axes([0.1,0.8,0.05,0.2])
        #cbar = plt.colorbar(sc,cmap=this_cmap,spacing='uniform', orientation='horizontal',cax=ax2)
        #cbar.set_label(colorbarlabel, fontsize=10)
        #ax2.set_xlabel(colorbarlabel, fontname='Arial',labelpad=5)


    ax.spines['top'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)

    ax.set_ylabel("Elevation (m)")

    # we need special formatting for the fow distance, since we want locations in kilometres
    if axis_data_name == 'flow_distance':
        # now get the tick marks
        n_target_tics = 5
        X_axis_min = 0
        xlocs,new_x_labels = LSDMap_BP.TickConverter(X_axis_min,X_axis_max,n_target_tics)

        ax.set_xticks(xlocs)

        ax.set_xticklabels(new_x_labels,rotation=60)

        ax.set_xlabel("Flow distance (km)")
    else:
        ax.set_xlabel("$\chi$ (m)")

    # This affects all axes because we set share_all = True.
    ax.set_ylim(z_axis_min,z_axis_max)
    ax.set_xlim(0,X_axis_max)

    # adjust the text
    #adjust_text(texts)

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
    ax.tick_params(axis='both', width=1, pad = 2)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(2)

    if rotate_labels:
        # remove the first y label. This is not working at the moment :( Fiona 30/03/18
        yticks = ax.yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)


    # Lets try to size the figure
    #cbar_L = "bottom"
    [fig_size_inches,map_axes,cbar_axes] = Helper.MapFigureSizer(fig_size_inches,aspect_ratio, cbar_loc = cbar_loc, title = "None")
    #print("Trying to use MapFigureSizer")
    #print(fig_size_inches)
    #print(map_axes)
    #print(cbar_axes)

    fig.set_size_inches(fig_size_inches[0], fig_size_inches[1])
    ax.set_position(map_axes)

    if cbar_loc != "None":
        ax2.set_position(cbar_axes)

    print("The figure format is: " + FigFormat)
    if FigFormat == 'show':
        plt.show()
    elif FigFormat == 'return':
        return fig
    else:
        plt.savefig(FigFileName,format=FigFormat,dpi=dpi)
        fig.clf()

def ChannelProfilePlot(DataDirectory, fname_prefix, FigFormat='png', size_format='ESURF',basin_key=[0], source_key=[0]):
    """
    This function makes a simple river long profile plot from the chi data map.

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the name of the DEM without extension
        FigFormat(str): format of the figure, e.g. png, svg
        size_format (str): size of the figure, can be either 'geomorphology', 'big', or 'ESURF'
        basin_key (list): basin keys to analyse
        source_key (list): source keys of the channels you want to plot.

    Returns:
        long profile plot

    Author: FJC
    """
    df = Helper.ReadChiDataMapCSV(DataDirectory,fname_prefix)

    # mask for the basin and the channel
    df = df[df['basin_key'].isin(basin_key)]
    df = df[df['source_key'].isin(source_key)]

    # set up the figure
    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        #l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        #l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
        #l_pad = -35

    gs = plt.GridSpec(100,100,bottom=0.2,left=0.05,right=0.95,top=0.95)
    ax = fig.add_subplot(gs[5:100,10:95])

    elevation = df['elevation'].tolist()
    flow_distance = df['flow_distance'].tolist()

    ax.plot(flow_distance, elevation, c='b')
    ax.set_xlabel('Distance upstream from outlet (m)')
    ax.set_ylabel('Elevation (m)')

    newFilename = DataDirectory+fname_prefix+"_profiles."+FigFormat
    plt.savefig(newFilename,format=FigFormat,dpi=300)
    ax.cla()
    plt.close(fig)


def map_Mchi_standard(DataDirectory, fname_prefix, size_format='ESURF', FigFormat='png', basin_list = [],outlier_detection_method = "None", log = False, colmanscal = [], bkbg = False, knickpoint = False):

    """
    This creates a basic knickpoint map

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_list (list): List of the basin ID you want.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        mancut (float): manual cutoff for the data selection.
        outlier_detection_method (str): determine the outlier detection method to select the right knickpoints
        colmanscal = manual color scale for the m_chi [min,max]
        bkbg (bool): turn to True to plot the data with a black background

    Returns:
        Shaded relief plot with the basins outlines and the knickpoints sized by intensity

    Author: BG - 05/10/2017
    """


    # check if a directory exists for the chi plots. If not then make it.
    raster_directory = DataDirectory+'raster_plots/'
    if not os.path.isdir(raster_directory):
        os.makedirs(raster_directory)

    # Set up fonts for plots
    basls = basin_list
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # set figure sizes based on format
    if size_format == "geomorphology":
        fig_width_inches = 6.25
    elif size_format == "big":
        fig_width_inches = 16
    else:
        fig_width_inches = 4.92126



    # get the rasters
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext
    BasinsName = fname_prefix+'_AllBasins'+raster_ext


    # create the map figure
    MF = MapFigure(HillshadeName, DataDirectory,coord_type="UTM_km", alpha = 0.7)
    if bkbg:
        MF.add_drape_image(HillshadeName,DataDirectory,colourmap = "gray",alpha=1,colour_min_max = [10000,10001],modify_raster_values=False,old_values=[], new_values=[],NFF_opti = True)
        # plot the basin outlines
        Basins = LSDP.GetBasinOutlines(DataDirectory, BasinsName)
        MF.plot_polygon_outlines(Basins, linewidth=0.5, colour = "white")
    else:
        # plot the basin outlines
        Basins = LSDP.GetBasinOutlines(DataDirectory, BasinsName)
        MF.plot_polygon_outlines(Basins, linewidth=0.5)

    # add the channel network without color
    if (knickpoint):
        knickpoint = "knickpoint"
    else:
        knickpoint = "normal"
    ChannelDF = Helper.ReadMChiSegCSV(DataDirectory,fname_prefix, type = knickpoint)
    ChannelPoints = LSDP.LSDMap_PointData(ChannelDF, data_type = "pandas", PANDEX = True)
    MF.add_point_data(ChannelPoints,this_colourmap = "RdBu_r", column_for_plotting = "m_chi",show_colourbar = True, colour_manual_scale = colmanscal, scale_points=True,scaled_data_in_log= True, column_for_scaling='drainage_area',alpha=0.4,max_point_size = 4,min_point_size = 1,zorder=100)

    # add the knickpoints plots
    if bkbg:
        ImageName = raster_directory+fname_prefix+"_Mchi_BK."+FigFormat
    else:
        ImageName = raster_directory+fname_prefix+"_Mchi_HS."+FigFormat
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 300) # Save the figure
