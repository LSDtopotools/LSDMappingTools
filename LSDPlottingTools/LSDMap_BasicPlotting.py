## LSDMap_BasicPlotting.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## These functions are tools to deal with rasters
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## SMM
## 26/07/2014
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from matplotlib import rcParams
from .adjust_text import adjust_text
import LSDPlottingTools.LSDMap_GDALIO as LSDMap_IO
import LSDPlottingTools.LSDMap_BasicManipulation as LSDMap_BM
import LSDPlottingTools.LSDMap_OSystemTools as LSDOst
from scipy import signal
import matplotlib.pyplot as plt
from LSDPlottingTools import colours


def TickSpineFormatter(ax, sizeformat = "esurf"):
    """This formats the line weights on the bounding box and ticks.

    Args:
        ax1 (axis object): the matplotlib axis object
        size_format (str): The size format. Can be geomorhpology, esurf or big

    returns:
        The axis object

    Author: SMM
    """

    import matplotlib.lines as mpllines

    # some formatting to make some of the ticks point outward
    for line in ax.get_xticklines():
        line.set_marker(mpllines.TICKDOWN)

    for line in ax.get_yticklines():
        line.set_marker(mpllines.TICKLEFT)

    if sizeformat == "esurf":
        lw = 1.0
        pd = 8
    elif sizeformat == "geomorphology":
        lw = 1.5
        pd = 10
    elif sizeformat == "big":
        lw = 2
        pd = 12
    else:
        lw = 1.0
        pd = 8

    ax.spines['top'].set_linewidth(lw)
    ax.spines['left'].set_linewidth(lw)
    ax.spines['right'].set_linewidth(lw)
    ax.spines['bottom'].set_linewidth(lw)

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
    ax.tick_params(axis='both', width=lw, pad = pd)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(pd)

    return ax



#==============================================================================
# This formats ticks if you want to convert metres to km
#==============================================================================
def TickLabelShortenizer(labels,n_hacked_digits):
    """This takes a list of labels and hacks off digits so you can go, say, from metres to km.

    Args:
        labels (str): A list of labels
        n_hacked_digits (int): The number of digits to hack off

    Return:
        The new labels

    Author: SMM
    """
    new_labels = []
    for label in labels:
        new_labels.append(label[0:-n_hacked_digits])

    return new_labels



def TickConverter(x_min,x_max,n_target_tics):
    """This function is used to convert ticks in metres to ticks in kilometres.

    Args:
        x_min (float): The minimum value on the axis (in metres).
        x_max (float): The maximum value on the axis (in metres)
        n_target_ticks (int): The number of ticks you want on the axis (this is optimised so you may not get exactly this number)

    Returns:
        new_xlocs (float list): List of locations of the ticks in metres.
        new_x_labels (str list): List of strings for ticks, will be location in kilometres.

    Author:
        Simon M Mudd
    """


    dx_fig = x_max-x_min
    dx_spacing = dx_fig/n_target_tics
    #print("spacing: "+str(dx_spacing))

    # This extracts the digits before the full stop
    str_dx = str(dx_spacing)
    str_dx = str_dx.split('.')[0]
    n_digits = str_dx.__len__()
    nd = int(n_digits)
    # We are left with the number of digits in the spacing. This will be used
    # to round tick locations

    first_digit = float(str_dx[0])

    dx_spacing_rounded = first_digit*pow(10,(nd-1))
    #print("dx spacing is: " + str(dx_spacing_rounded))

    str_xmin = str(x_min)
    #print("before split str_xmin: "+ str_xmin)
    str_xmin = str_xmin.split('.')[0]
    #print("after split str_xmin: "+ str_xmin)
    x_min = float(str_xmin)
    #print("x_min: "+ str(x_min))

    n_digx = str_xmin.__len__()

    if (n_digx-nd+1) >= 1:
        front_x = str_xmin[:(n_digx-nd+1)]
    else:
        front_x = str_xmin

    round_xmin = float(front_x)*pow(10,nd-1)
    #print("round xmin is: " + str(round_xmin))
    if round_xmin <0:
        round_xmin = 0

    # now we need to figure out where the xllocs and ylocs are
    xlocs = np.zeros(2*n_target_tics)
    xlocs_km = np.zeros(2*n_target_tics)
    new_x_labels = []

    for i in range(0,2*n_target_tics):
        xlocs[i] = round_xmin+(i)*dx_spacing_rounded
        xlocs_km[i] = xlocs[i]/1000.0
        new_x_labels.append( str(xlocs_km[i]).split(".")[0] )

    #print xlocs
    #print new_x_labels
    #print new_y_labels

    new_xlocs = []
    new_xlocs_km = []
    x_labels = []

    # Now loop through these to get rid of those not in range
    for index,xloc in enumerate(xlocs):
        #print xloc
        if (xloc <= x_max and xloc >= x_min):
            new_xlocs.append(xloc)
            new_xlocs_km.append(xlocs_km[index])
            x_labels.append(new_x_labels[index])


    #print "======================================="
    #print "I am getting the tick marks now"
    #print "X extent: " + str(x_min)+ " " +str(x_max)
    #print "x ticks: "
    #print new_xlocs
    #print x_labels

    #return xlocs,ylocs,new_x_labels,new_y_labels
    return new_xlocs,x_labels


#==============================================================================
# Formats ticks for an imshow plot in UTM
# Filename is the name of the file with full path
# x_max, x_min, y_max, y_min are the extent of the plotting area (NOT the DEM)
# n_target ticks are the number of ticks for plotting
#------------------------------------------------------------------------------
def GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics):
    """This fuction is used to set tick locations for UTM maps. It tries to optimise the spacing of these ticks.

    Args:
        x_min (float): The minimum value on the x axis (in metres).
        x_max (float): The maximum value on the x axis (in metres).
        y_min (float): The minimum value on the y axis (in metres).
        y_max (float): The maximum value on the y axis (in metres).
        n_target_ticks (int): The number of ticks you want on the axis (this is optimised so you may not get exactly this number)

    Returns:
        new_xlocs (float list): List of locations of the ticks in metres.
        new_x_labels (str list): List of strings for ticks, will be location in kilometres.
        new_ylocs (float list): List of locations of the ticks in metres.
        new_y_labels (str list): List of strings for ticks, will be location in kilometres.

    Author: SMM
    """

    CellSize,XMin,XMax,YMin,YMax = LSDMap_IO.GetUTMMaxMin(FileName)
    NDV, xsize, ysize, GeoT, Projection, DataType = LSDMap_IO.GetGeoInfo(FileName)

    #print("Getting ticks. YMin: "+str(YMin)+" and YMax: "+str(YMax))

    xmax_UTM = XMax
    xmin_UTM = XMin

    ymax_UTM = YMax
    ymin_UTM = YMin

    dy_fig = ymax_UTM-ymin_UTM
    dx_fig = xmax_UTM-xmin_UTM

    dx_spacing = dx_fig/n_target_tics
    dy_spacing = dy_fig/n_target_tics

    if (dx_spacing>dy_spacing):
        dy_spacing = dx_spacing

    str_dy = str(dy_spacing)
    str_dy = str_dy.split('.')[0]
    n_digits = str_dy.__len__()
    nd = int(n_digits)

    first_digit = float(str_dy[0])

    dy_spacing_rounded = first_digit*pow(10,(nd-1))

    str_xmin = str(xmin_UTM)
    str_ymin = str(ymin_UTM)
    str_xmin = str_xmin.split('.')[0]
    str_ymin = str_ymin.split('.')[0]
    xmin_UTM = float(str_xmin)
    ymin_UTM = float(str_ymin)



    n_digx = str_xmin.__len__()
    n_digy = str_ymin.__len__()


    if (n_digx-nd+1) >= 1:
        front_x = str_xmin[:(n_digx-nd+1)]
    else:
        front_x = str_xmin

    if (n_digy-nd+1) >= 1:
        front_y = str_ymin[:(n_digy-nd+1)]
    else:
        front_y = str_ymin


    round_xmin = float(front_x)*pow(10,nd-1)
    round_ymin = float(front_y)*pow(10,nd-1)

    #print("UTM y in: "+str(ymin_UTM)+" and rounded min: "+str(round_ymin))


    # now we need to figure out where the xllocs and ylocs are
    xUTMlocs = np.zeros(2*n_target_tics)
    yUTMlocs = np.zeros(2*n_target_tics)
    xlocs = np.zeros(2*n_target_tics)
    ylocs = np.zeros(2*n_target_tics)

    new_x_labels = []
    new_y_labels = []

    for i in range(0,2*n_target_tics):

        # Note we use dy spacing here in both x and y directions since we want
        # Ticks spaced the same in each direction
        xUTMlocs[i] = round_xmin+(i)*dy_spacing_rounded
        yUTMlocs[i] = round_ymin+(i)*dy_spacing_rounded
        xlocs[i] = xUTMlocs[i]

        # need to account for the rows starting at the upper boundary
        ylocs[i] = YMax-(yUTMlocs[i]-YMin)

        new_x_labels.append( str(xUTMlocs[i]).split(".")[0] )
        new_y_labels.append( str(yUTMlocs[i]).split(".")[0] )


    new_xlocs = []
    new_xUTMlocs = []
    x_labels = []

    # Now loop through these to get rid of those not in range
    for index,xloc in enumerate(xlocs):
        if (xloc < XMax and xloc > XMin):
            new_xlocs.append(xloc)
            new_xUTMlocs.append(xUTMlocs[index])
            x_labels.append(new_x_labels[index])

    new_ylocs = []
    new_yUTMlocs = []
    y_labels = []

    # Now loop through these to get rid of those not in range
    #you have to reverse the order of the lists.

    #print("before reverse: ")
    #print(ylocs)
    ylocs = ylocs[::-1]
    #print("After")
    #print(ylocs)
    yUTMlocs = yUTMlocs[::-1]
    new_y_labels = new_y_labels[::-1]

    #print("UTM y: ")
    #print(yUTMlocs)
    #print("and y locs:")
    #print(ylocs)
    #print("And labels:")
    #print(new_y_labels)


    for index,yloc in enumerate(ylocs):
        UTMloc = yUTMlocs[index]
        if (UTMloc < YMax and UTMloc > YMin):
            #print("UTMloc: "+str(UTMloc)+" yloc: "+str(yloc)+" label: "+new_y_labels[index])
            new_ylocs.append(yloc)
            new_yUTMlocs.append(yUTMlocs[index])
            y_labels.append(new_y_labels[index])


    #return xlocs,ylocs,new_x_labels,new_y_labels
    return new_xlocs,new_ylocs,x_labels,y_labels
#==============================================================================

#==============================================================================
# Formats ticks for an imshow plot in UTM
# Filename is the name of the file with full path
# x_max, x_min, y_max, y_min are the extent of the plotting area (NOT the DEM)
# n_target ticks are the number of ticks for plotting
#------------------------------------------------------------------------------
def GetTicksForUTMNoInversion(FileName,x_max,x_min,y_max,y_min,n_target_tics,minimum_tick_spacing=0):
    """This fuction is used to set tick locations for UTM maps. It tries to optimise the spacing of these ticks.

    Args:
        x_min (float): The minimum value on the x axis (in metres).
        x_max (float): The maximum value on the x axis (in metres).
        y_min (float): The minimum value on the y axis (in metres).
        y_max (float): The maximum value on the y axis (in metres).
        n_target_ticks (int): The number of ticks you want on the axis (this is optimised so you may not get exactly this number)
        minimum_tick_spacing (int): The minimum spacing between ticks to impose (e.g. 1000 to get ticks every km)

    Returns:
        new_xlocs (float list): List of locations of the ticks in metres.
        new_x_labels (str list): List of strings for ticks, will be location in kilometres.
        new_ylocs (float list): List of locations of the ticks in metres.
        new_y_labels (str list): List of strings for ticks, will be location in kilometres.

    Author: SMM
    """

    CellSize,XMin,XMax,YMin,YMax = LSDMap_IO.GetUTMMaxMin(FileName)
    NDV, xsize, ysize, GeoT, Projection, DataType = LSDMap_IO.GetGeoInfo(FileName)
    
    # take min and max specified as input arguments (MDH addition, these arguments were not otherwise used).
    Extents = [x_min,x_max,y_min,y_max]
    if None not in Extents:
      XMin = x_min
      XMax = x_max
      YMin = y_min
      YMax = y_max
      
    #print("Getting ticks. YMin: "+str(YMin)+" and YMax: "+str(YMax))

    xmax_UTM = XMax
    xmin_UTM = XMin

    ymax_UTM = YMax
    ymin_UTM = YMin

    dy_fig = ymax_UTM-ymin_UTM
    dx_fig = xmax_UTM-xmin_UTM

    dx_spacing = dx_fig/n_target_tics
    dy_spacing = dy_fig/n_target_tics

    if (dx_spacing>dy_spacing):
        dy_spacing = dx_spacing

    str_dy = str(dy_spacing)
    str_dy = str_dy.split('.')[0]
    n_digits = str_dy.__len__()
    nd = int(n_digits)

    first_digit = float(str_dy[0])

    dy_spacing_rounded = first_digit*pow(10,(nd-1))

    str_xmin = str(xmin_UTM)
    str_ymin = str(ymin_UTM)
    str_xmin = str_xmin.split('.')[0]
    str_ymin = str_ymin.split('.')[0]
    xmin_UTM = float(str_xmin)
    ymin_UTM = float(str_ymin)

    print("minimum values are x: "+str(xmin_UTM)+ " and y: "+str(ymin_UTM))

    n_digx = str_xmin.__len__()
    n_digy = str_ymin.__len__()


    if (n_digx-nd+1) >= 1:
        front_x = str_xmin[:(n_digx-nd+1)]
    else:
        front_x = str_xmin

    if (n_digy-nd+1) >= 1:
        front_y = str_ymin[:(n_digy-nd+1)]
    else:
        front_y = str_ymin
    ###### Sort a strange bug, but will create anormal values for scale
    if(front_x == "-"):
        front_x = 0
    ###### I'll try to sort it when I'll have time but contact me in case you need it before (Boris)
    round_xmin = float(front_x)*pow(10,nd-1)
    round_ymin = float(front_y)*pow(10,nd-1)

    #print("UTM y in: "+str(ymin_UTM)+" and rounded min: "+str(round_ymin))


    # now we need to figure out where the xllocs and ylocs are
    xUTMlocs = np.zeros(2*n_target_tics)
    yUTMlocs = np.zeros(2*n_target_tics)
    xlocs = np.zeros(2*n_target_tics)
    ylocs = np.zeros(2*n_target_tics)

    new_x_labels = []
    new_y_labels = []

    for i in range(0,2*n_target_tics):

        # Note we use dy spacing here in both x and y directions since we want
        # Ticks spaced the same in each direction
        xUTMlocs[i] = round_xmin+(i)*dy_spacing_rounded
        yUTMlocs[i] = round_ymin+(i)*dy_spacing_rounded
        xlocs[i] = xUTMlocs[i]
        ylocs[i] = yUTMlocs[i]

        new_x_labels.append( str(xUTMlocs[i]).split(".")[0] )
        new_y_labels.append( str(yUTMlocs[i]).split(".")[0] )


    new_xlocs = []
    new_xUTMlocs = []
    x_labels = []

    # Now loop through these to get rid of those not in range
    for index,xloc in enumerate(xlocs):
        if (xloc < XMax and xloc > XMin):
            new_xlocs.append(xloc)
            new_xUTMlocs.append(xUTMlocs[index])
            x_labels.append(new_x_labels[index])

    new_ylocs = []
    new_yUTMlocs = []
    y_labels = []



    for index,yloc in enumerate(ylocs):
        if (yloc < YMax and yloc > YMin):
            new_ylocs.append(yloc)
            new_yUTMlocs.append(yUTMlocs[index])
            y_labels.append(new_y_labels[index])


    #return xlocs,ylocs,new_x_labels,new_y_labels
    return new_xlocs,new_ylocs,x_labels,y_labels
#==============================================================================

#==============================================================================
def BasicDensityPlot(FileName, thiscmap='gray',colorbarlabel='Elevation in meters',
                             clim_val = (0,0),FigFileName = 'Image.pdf', FigFormat = 'show',
                             size_format = "esurf", is_log = False):
    """This creates a plot of a raster. The most basic plotting function. It uses AxisGrid to ensure proper placment of the raster.

    Args:
        FileName (str): The name of the raster (with full path and extension).
        thiscmap (colormap): The colourmap to be used.
        colorbarlabel (str): The label of the colourbar
        clim_val (float,float): The colour limits. If (0,0) then the min and max raster values are used.
        FigFilename (str): The name of the figure (with extension)
        FigFormat (str): the format of the figure (e.g., jpg, png, pdf). If "show" then the figure is plotted to screen.
        size_format (str): the size of the figure
        islog (bool): True if you want the figure to have a log density not the density

    Returns:
        A density plot of the raster

    Author:
        Simon M Mudd
    """

    import matplotlib.pyplot as plt

    # Set up fonts for plots
    label_size = 20

    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # get the data
    raster = LSDMap_IO.ReadRasterArrayBlocks(FileName)

    if is_log:
        # get the log of the raster
        raster = np.log10(raster)

    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(FileName)

    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # make a figure,
    # make a figure,
    print("The size format is: " + size_format)
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    gs = plt.GridSpec(100,100,bottom=0.1,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[15:100,15:95])

    # This is the axis for the colorbar
    ax2 = fig.add_subplot(gs[10:15,15:80])


    # now get the tick marks
    n_target_tics = 5
    xlocs,ylocs,new_x_labels,new_y_labels = GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics)

    im = ax.imshow(raster[::-1], thiscmap, extent = extent_raster, interpolation="nearest")


    plt.colorbar(im,cmap=thiscmap,spacing='uniform', orientation='horizontal',cax=ax2)
    ax2.text(0,1,colorbarlabel,horizontalalignment='left',
        verticalalignment='bottom',transform=ax2.transAxes)
    #ax2.set_xlabel(colorbarlabel, fontname='Arial',labelpad=0)
    #cbar = plt.colorbar(im)
    #cbar.set_label(colorbarlabel)

    # set the colour limits
    #print "Setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1])
    if (clim_val == (0,0)):
        #print "I don't think I should be here"
        im.set_clim(0, np.nanmax(raster))
    else:
        print("Now setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1]))
        im.set_clim(clim_val[0],clim_val[1])


    # This affects all axes because we set share_all = True.
    ax.set_xlim(x_min,x_max)
    ax.set_ylim(y_max,y_min)

    ax.set_xticks(xlocs)
    ax.set_yticks(ylocs)

    # convert to km
    n_hacked_digits = 3
    new_x_labels = TickLabelShortenizer(new_x_labels,n_hacked_digits)
    new_y_labels = TickLabelShortenizer(new_y_labels,n_hacked_digits)

    ax.set_xticklabels(new_x_labels,rotation=60)
    ax.set_yticklabels(new_y_labels)

    ax.set_xlabel("Easting (km)")
    ax.set_ylabel("Northing (km)")

    ax = TickSpineFormatter(ax,size_format)

#    # go through the ticks
#    ax.spines['top'].set_linewidth(2.5)
#    ax.spines['left'].set_linewidth(2.5)
#    ax.spines['right'].set_linewidth(2.5)
#    ax.spines['bottom'].set_linewidth(2.5)
#
#    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
#    ax.tick_params(axis='both', width=2.5, pad = 10)
#    for tick in ax.xaxis.get_major_ticks():
#        tick.set_pad(10)

    print("The figure format is: " + FigFormat)
    if FigFormat == 'show':
        plt.show()
    elif FigFormat == 'return':
        return fig
    else:
        plt.savefig(FigFileName,format=FigFormat)
        fig.clf()


#==============================================================================

#==============================================================================
def BasicDrapedPlotGridPlot(FileName, DrapeName, thiscmap='gray',drape_cmap='gray',
                            colorbarlabel='Elevation in meters',clim_val = (0,0),
                            drape_alpha = 0.6,FigFileName = 'Image.pdf',FigFormat = 'show', dpi_save = 250):
    """This creates a draped plot of a raster. It uses AxisGrid to ensure proper placment of the raster.

    Args:
        FileName (str): The name of the raster (with full path and extension).
        DrapeName (str): The name of the drape raster (with full path and extension). If the DrapeName is "None" it will calculate the hillshade.
        thiscmap (colormap): The colourmap to be used.
        drape_cmap (colormap): The colourmap to be used for the drape.
        colorbarlabel (str): The label of the colourbar
        clim_val (float,float): The colour limits. If (0,0) then the min and max raster values are used.
        drape_alpha (float): The alpha value (transparency) of the drape
        FigFilename (str): The name of the figure (with extension)
        FigFormat (str): the format of the figure (e.g., jpg, png, pdf). If "show" then the figure is plotted to screen.

    Returns:
        A density plot of the draped raster

    Author:
        Simon M Mudd
    """


    import matplotlib.pyplot as plt

    # Set up fonts for plots
    label_size = 20
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # get the data
    raster = LSDMap_IO.ReadRasterArrayBlocks(FileName)

    from scipy import ndimage
    if DrapeName == "None":
        #filtered = ndimage.filters.gaussian_filter(raster, 3)
        #filtered = signal.wiener(raster)
        raster_drape = Hillshade(raster)
    else:
        raster_drape = LSDMap_IO.ReadRasterArrayBlocks(DrapeName)

    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(FileName)

    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # make a figure, sized for a ppt slide
    fig = plt.figure(1, facecolor='white',figsize=(10,7.5))

    gs = plt.GridSpec(100,75,bottom=0.1,left=0.1,right=0.9,top=1.0)
    ax = fig.add_subplot(gs[10:100,10:75])

    # now get the tick marks
    n_target_tics = 5
    xlocs,ylocs,new_x_labels,new_y_labels = GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics)

    im1 = ax.imshow(raster[::-1], thiscmap, extent = extent_raster, interpolation="nearest")

    cbar = plt.colorbar(im1)
    cbar.set_label(colorbarlabel)

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
    im3 = ax.imshow(raster_drape[::-1], drape_cmap, extent = extent_raster, alpha = drape_alpha, interpolation="nearest")

    # Set the colour limits of the drape
    im3.set_clim(0,np.nanmax(raster_drape))


    ax.spines['top'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)

    # This affects all axes because we set share_all = True.
    ax.set_xlim(x_min,x_max)
    ax.set_ylim(y_max,y_min)

    ax.set_xticks(xlocs)
    ax.set_yticks(ylocs)

    ax.set_xticklabels(new_x_labels,rotation=60)
    ax.set_yticklabels(new_y_labels)

    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap
    ax.tick_params(axis='both', width=1.5, pad = 10)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(10)

    print("The figure format is: " + FigFormat)
    if FigFormat == 'show':
        plt.show()
    elif FigFormat == 'return':
        return fig
    else:
        plt.savefig(FigFileName,format=FigFormat,dpi=dpi_save)
        fig.clf()

#==============================================================================


#==============================================================================
def DrapedOverHillshade(FileName, DrapeName, thiscmap='gray',drape_cmap='gray',
                            colorbarlabel='Elevation in meters',clim_val = (0,0),
                            drape_alpha = 0.6, ShowColorbar = False,
                            ShowDrapeColorbar=False, drape_cbarlabel=None):
    """This creates a draped plot of a raster.

    It uses AxisGrid to ensure proper placment of the raster.
    It also includes a hillshde to make the figure look nicer (so there are three raster layers).

    Note:
        Remember, this has THREE layers: a base layer, a hillshade and a drape.

    Args:
        FileName (str): The name of the raster (with full path and extension).
        DrapeName (str): The name of the drape raster (with full path and extension).
        thiscmap (colormap): The colourmap to be used.
        drape_cmap (colormap): The colourmap to be used for the drape.
        colorbarlabel (str): The label of the colourbar
        clim_val (float,float): The colour limits. If (0,0) then the min and max raster values are used.
        drape_alpha (float): The alpha value (transparency) of the drape
        ShowColorbar (bool): Whether you want to show the colorbar
        drape_cbarlabel (str): The label of the drape colourbar

    Returns:
        A density plot of the draped raster

    Author:
        SMM and DAV

    """


    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import AxesGrid

    label_size = 20


    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    hillshade = Hillshade(FileName)

    # DAV - option to supply array directly (after masking for example, rather
    # than reading directly from a file. Should not break anyone's code)
    # (You can't overload functions in Python...)
    if isinstance(DrapeName, str):
      raster_drape = LSDMap_IO.ReadRasterArrayBlocks(DrapeName)
    elif isinstance(DrapeName, np.ndarray):
      raster_drape = DrapeName
    else:
      print("DrapeName supplied is of type: ", type(DrapeName))
      raise ValueError('DrapeName must either be a string to a filename, \
      or a numpy ndarray type. Please try again.')

    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(FileName)

    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # make a figure, sized for a ppt slide
    fig = plt.figure(1, facecolor='white',figsize=(10,7.5))

    if ShowColorbar:
        grid = AxesGrid(fig, 111,
                        nrows_ncols=(1, 1),
                        axes_pad=(0.45, 0.15),
                        label_mode="1",
                        share_all=True,
                        cbar_location="right",
                        cbar_mode="each",
                        cbar_size="7%",
                        cbar_pad="2%",
                        )

    else:
        grid = AxesGrid(fig, 111,
                        nrows_ncols=(1, 1),
                        axes_pad=(0.45, 0.15),
                        label_mode="1",
                        share_all=True,
                        )

    # now get the tick marks
    n_target_tics = 5
    xlocs,ylocs,new_x_labels,new_y_labels = GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics)

    im = grid[0].imshow(hillshade[::-1], thiscmap, extent = extent_raster, interpolation="nearest")
    #im = grid[0].imshow(raster, thiscmap, interpolation="nearest")
    if ShowColorbar:
        cbar = grid.cbar_axes[0].colorbar(im)
        cbar.set_label_text(colorbarlabel)

    # set the colour limits
    print("Setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1]))
    if (clim_val == (0,0)):
        print("I don't think I should be here")
        im.set_clim(0, np.nanmax(hillshade))
    else:
        print("Now setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1]))
        im.set_clim(clim_val[0],clim_val[1])

    # Now for the drape: it is in grayscape
    im2 = grid[0].imshow(raster_drape[::-1], drape_cmap, extent = extent_raster, alpha = drape_alpha, interpolation="none")

    if ShowDrapeColorbar:
      cbar2 = grid.cbar_axes[0].colorbar(im2)
      cbar2.set_label_text(drape_cbarlabel)


    # This affects all axes because we set share_all = True.
    grid.axes_llc.set_xlim(x_min,x_max)
    grid.axes_llc.set_ylim(y_max,y_min)

    grid.axes_llc.set_xticks(xlocs)
    grid.axes_llc.set_yticks(ylocs)

    grid.axes_llc.set_xticklabels(new_x_labels,rotation=60)
    grid.axes_llc.set_yticklabels(new_y_labels)

    grid.axes_llc.set_xlabel("Easting (m)")
    grid.axes_llc.set_ylabel("Northing (m)")

    plt.show()

# =============================================================================
def DrapedOverHillshade_Categories(FileName, DrapeName, nCategories,
                                   category_min_max, thiscmap='gray',
                                   drape_cmap='gray', clim_val=(0, 0),
                                   drape_alpha=0.6, ShowDrapeColorbar=False,
                                   drape_cbarlabel=None, category_labels=None):
    """This creates a draped plot of a categorical raster.

    It uses AxisGrid to ensure proper placment of the raster.
    It also includes a hillshade to make the figure look nicer
    (so there are three raster layers).

    Note:
        Remember, this has THREE layers: a base layer, a hillshade and a drape.

    Args:
        FileName (str): The name of the raster (with full path and extension).
        DrapeName (str): The name of the drape raster (with full path and extension).
        nCategories (int): The number of categories in the dataset
        category_min_max (float, float): The min and max values in the categorical raster.
        thiscmap (colormap): The colourmap to be used for the hillshade.
        drape_cmap (colormap): The colourmap to be used for the drape.
        clim_val (float,float): The colour limits. If (0,0) then the min and max raster values are used.
        drape_alpha (float): The alpha value (transparency) of the drape.
        ShowDrapeColorbar (bool): Toggles the display of a categorical colorbar.
        drape_cbarlabel (str): The label of the drape colourbar.
        category_labels (list): List of strings used as category labels.

    Returns:
        A density plot of the draped categorical raster

    Author:
        SWDG, after SMM and DAV

    """

    from mpl_toolkits.axes_grid1 import AxesGrid

    label_size = 20

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    hillshade = Hillshade(FileName)

    # DAV - option to supply array directly (after masking for example, rather
    # than reading directly from a file. Should not break anyone's code)
    # (You can't overload functions in Python...)
    if isinstance(DrapeName, str):
        raster_drape = LSDMap_IO.ReadRasterArrayBlocks(DrapeName)
    elif isinstance(DrapeName, np.ndarray):
        raster_drape = DrapeName
    else:
        print("DrapeName supplied is of type: ", type(DrapeName))
        raise ValueError('DrapeName must either be a string to a filename,'
                         ' or a numpy ndarray type. Please try again.')

    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(FileName)

    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # make a figure, sized for a ppt slide
    fig = plt.figure(1, facecolor='white', figsize=(10, 7.5))

    if ShowDrapeColorbar:
        grid = AxesGrid(fig, 111,
                        nrows_ncols=(1, 1),
                        axes_pad=(0.45, 0.15),
                        label_mode="1",
                        share_all=True,
                        cbar_location="right",
                        cbar_mode="each",
                        cbar_size="7%",
                        cbar_pad="2%",
                        )

    else:
        grid = AxesGrid(fig, 111,
                        nrows_ncols=(1, 1),
                        axes_pad=(0.45, 0.15),
                        label_mode="1",
                        share_all=True,
                        )

    # now get the tick marks
    n_target_tics = 5
    xlocs, ylocs, new_x_labels, new_y_labels = GetTicksForUTM(FileName,
                                                              x_max,
                                                              x_min,
                                                              y_max,
                                                              y_min,
                                                              n_target_tics)

    im = grid[0].imshow(hillshade[::-1], thiscmap, extent=extent_raster,
                        interpolation="nearest")

    # set the colour limits
    if (clim_val == (0, 0)):
        im.set_clim(0, np.nanmax(hillshade))
    else:
        im.set_clim(clim_val[0], clim_val[1])

    # Now for the drape: it is in grayscape
    im2 = grid[0].imshow(raster_drape[::-1], drape_cmap, extent=extent_raster,
                         alpha=drape_alpha, interpolation="none")

    if ShowDrapeColorbar:
        cbar = colours.colorbar_index(plt.gcf(), grid.cbar_axes[0], nCategories,
                                      drape_cmap, category_min_max[0],
                                      category_min_max[1])
        cbar.set_label(drape_cbarlabel)
        if category_labels:
            cbar.set_ticklabels(category_labels)

    # This affects all axes because we set share_all = True.
    grid.axes_llc.set_xlim(x_min, x_max)
    grid.axes_llc.set_ylim(y_max, y_min)

    grid.axes_llc.set_xticks(xlocs)
    grid.axes_llc.set_yticks(ylocs)

    grid.axes_llc.set_xticklabels(new_x_labels, rotation=60)
    grid.axes_llc.set_yticklabels(new_y_labels)

    grid.axes_llc.set_xlabel("Easting (m)")
    grid.axes_llc.set_ylabel("Northing (m)")

    plt.show()


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def DrapedOverFancyHillshade(FileName, HSName, DrapeName, thiscmap='gray',drape_cmap='gray',
                            colorbarlabel='Basin number',clim_val = (0,0),
                            drape_alpha = 0.6,FigFileName = 'Image.pdf',FigFormat = 'show',
                            elevation_threshold = 0):
    """This creates a draped plot of a raster. It uses AxisGrid to ensure proper placment of the raster. It also includes a hillshde to make the figure look nicer (so there are three raster layers). In this case you need to tell it the name of the hillshade raster.

    Args:
        FileName (str): The name of the raster (with full path and extension).
        HSName (str): The name of the hillshade raster (with full path and extension).
        DrapeName (str): The name of the drape raster (with full path and extension).
        thiscmap (colormap): The colourmap to be used.
        drape_cmap (colormap): The colourmap to be used for the drape.
        colorbarlabel (str): The label of the colourbar
        clim_val (float,float): The colour limits. If (0,0) then the min and max raster values are used.
        drape_alpha (float): The alpha value (transparency) of the drape
        FigFilename (str): The name of the figure (with extension)
        FigFormat (str): the format of the figure (e.g., jpg, png, pdf). If "show" then the figure is plotted to screen.
        elevation_threshold (float): If raster values are less than this threshold they become nodata.

    Returns:
        A density plot of the draped raster

    Author:
        SMM

    """
    label_size = 10

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # get the data
    raster = LSDMap_IO.ReadRasterArrayBlocks(FileName)
    raster_HS = LSDMap_IO.ReadRasterArrayBlocks(HSName)
    raster_drape = LSDMap_IO.ReadRasterArrayBlocks(DrapeName)

    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(FileName)

    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # make a figure, sized for a ppt slide
    fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))

    gs = plt.GridSpec(100,100,bottom=0.25,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])

    # This is the axis for the colorbar
    ax2 = fig.add_subplot(gs[10:15,15:70])

    # now get the tick marks
    n_target_tics = 5
    xlocs,ylocs,new_x_labels,new_y_labels = GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics)
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
    HS_cmap = 'gray'
    im2 = ax.imshow(raster_HS[::-1], HS_cmap, extent = extent_raster, alpha = 0.4, interpolation="nearest")

    # Set the colour limits of the drape
    im2.set_clim(0,np.nanmax(raster_HS))

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

    # This affects all axes because we set share_all = True.
    ax.set_xlim(x_min,x_max)
    ax.set_ylim(y_max,y_min)

    ax.set_xticks(xlocs)
    ax.set_yticks(ylocs)

    cbar = plt.colorbar(im3,cmap=drape_cmap,spacing='uniform', orientation='horizontal',cax=ax2)
    cbar.set_label(colorbarlabel, fontsize=10)
    ax2.set_xlabel(colorbarlabel, fontname='Arial',labelpad=-35)

    print("The figure format is: " + FigFormat)
    if FigFormat == 'show':
        plt.show()
    elif FigFormat == 'return':
        return fig
    else:
        plt.savefig(FigFileName,format=FigFormat,dpi=500)
        fig.clf()

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## This function plots the chi slope on a shaded relief map
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def BasinsOverFancyHillshade(FileName, HSName, BasinName, Basin_csv_name, basin_point_data, thiscmap='gray',drape_cmap='gray',
                             clim_val = (0,0), drape_alpha = 0.6,FigFileName = 'Image.pdf',
                             FigFormat = 'show',elevation_threshold = 0,
                             grouped_basin_list = [], basin_rename_list = [],spread  = 20,
                             chanPointData = "None",
                             label_sources = False, source_chi_threshold = 10,
                             size_format = "esurf"):
    """This creates a plot with a hillshade draped over elevation (or any other raster) with the basins on them.

    Args:
        FileName (str): The name of the raster (with full path and extension).
        HSName (str): The name of the hillshade raster (with full path and extension).
        BasinName (str): The name of the basin raster (with full path and extension).
        Basin_csv_name (str): The name of the csv file where basin info is stored
        basin_point_data (LSDMap_PointData): The basin point data
        thiscmap (colormap): The colourmap to be used.
        drape_cmap (colormap): The colourmap to be used for the drape.
        clim_val (float,float): The colour limits. If (0,0) then the min and max raster values are used.
        drape_alpha (float): The alpha value (transparency) of the drape
        FigFilename (str): The name of the figure (with extension)
        FigFormat (str): the format of the figure (e.g., jpg, png, pdf). If "show" then the figure is plotted to screen.
        elevation_threshold (float): If raster values are less than this threshold they become nodata.
        grouped_basin_list (int list): A list of lists with basins to be grouped.
        basin_rename_list (int list): A list of updated names for the basins. So if you wanted basin 4 to be renamed basin 6 the fourth element in this list would be 6.
        spread (float): Basins get a different number each, this is the spread between groups that controls how different the grouped basins are.
        chanPointData (str ir LSDMap_PointData): Either "none" or a point data object with the channel network
        label_sources (bool): Wheteher or not to label the sources
        source_chi_threshold (float): Sources with length less than this will be plotted.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).

    Returns:
        A density plot of the draped raster

    Author:
        SMM

    """
    from . import LSDMap_ChiPlotting as LSDMap_CP
    label_size = 10

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # get the data
    raster = LSDMap_IO.ReadRasterArrayBlocks(FileName)
    raster_HS = LSDMap_IO.ReadRasterArrayBlocks(HSName)
    raster_drape = LSDMap_IO.ReadRasterArrayBlocks(BasinName)

    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(FileName)

    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    #x_range = x_max-x_min
    #y_range = y_max-y_min

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
    xlocs,ylocs,new_x_labels,new_y_labels = GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics)
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
    HS_cmap = 'gray'
    im2 = ax.imshow(raster_HS[::-1], HS_cmap, extent = extent_raster, alpha = 0.4, interpolation="nearest")

    # Set the colour limits of the drape
    im2.set_clim(0,np.nanmax(raster_HS))

    # now we need to update the basin colours
    # this groups basins
    if grouped_basin_list:
        key_groups = LSDMap_BM.BasinKeyToJunction(grouped_basin_list,basin_point_data)
        raster_drape = LSDMap_BM.RedefineIntRaster(raster_drape,key_groups,spread)


    # Now for the drape: it is in grayscale
    #print "drape_cmap is: "+drape_cmap
    im3 = ax.imshow(raster_drape[::-1], drape_cmap, extent = extent_raster, alpha = drape_alpha, interpolation="nearest")

    # Set the colour limits of the drape
    im3.set_clim(0,np.nanmax(raster_drape))

    #========================================================
    # now we need to label the basins
    # Now we get the chi points
    #print("Getting the EPSG string from: "+FileName)
    #print("Type of fname: "+str(type(FileName)))
    FileName= str(FileName)
    #print("Now the type is: "+str(type(FileName)))
    EPSG_string = LSDMap_IO.GetUTMEPSG(FileName)
    #print("EPSG string is: " + EPSG_string)

    # Now plot channel data
    if chanPointData != "None":

        # convert to easting and northing
        [easting_c,northing_c] = chanPointData.GetUTMEastingNorthing(EPSG_string)

        # The image is inverted so we have to invert the northing coordinate
        Ncoord_c = np.asarray(northing_c)
        Ncoord_c = np.subtract(extent_raster[3],Ncoord_c)
        Ncoord_c = np.add(Ncoord_c,extent_raster[2])

        # plot the rivers
        ax.scatter(easting_c,Ncoord_c,s=0.2, marker='.',color= 'b',alpha=0.2)


    thisPointData = basin_point_data

    # convert to easting and northing
    [easting,northing] = thisPointData.GetUTMEastingNorthingFromQuery(EPSG_string,"outlet_latitude","outlet_longitude")

    # The image is inverted so we have to invert the northing coordinate
    Ncoord = np.asarray(northing)
    Ncoord = np.subtract(extent_raster[3],Ncoord)
    Ncoord = np.add(Ncoord,extent_raster[2])

    these_data = thisPointData.QueryData("outlet_junction")
    #print M_chi
    these_data = [int(x) for x in these_data]

    # plot the centroids
    ax.scatter(easting,Ncoord,s=1, marker='o',color= 'r',alpha=0.7)

    # add a load of xmin, xmax points
    buffered_east = []
    buffered_north = []

    for loc in easting:
        buffered_east.append(loc)
        buffered_east.append(loc)
        buffered_east.append(loc)
        buffered_east.append(extent_raster[0])
        buffered_east.append(extent_raster[1])

    for loc in Ncoord:
        buffered_north.append(loc)
        buffered_north.append(extent_raster[2])
        buffered_north.append(extent_raster[3])
        buffered_north.append(loc)
        buffered_north.append(loc)


    # add text for basins
    texts = []
    bbox_props = dict(boxstyle="circle,pad=0.1", fc="w", ec="k", lw=0.5,alpha = 0.5)
    for index, datum in enumerate(these_data):
        this_easting = easting[index]
        this_northing = Ncoord[index]

        # Check to see if basins rename list works
        if basin_rename_list:
            if len(basin_rename_list) == len(these_data):
                texts.append(ax.text(this_easting,this_northing, str(basin_rename_list[index]),fontsize = 8, color= "r",alpha=0.7,bbox=bbox_props))
        else:
            texts.append(ax.text(this_easting,this_northing, str(index),fontsize = 8, color= "r",alpha=0.7,bbox=bbox_props))



    # Now for sources labelling
    texts2 = []
    if label_sources:

        # Format the bounding box
        bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="b", lw=0.5,alpha = 0.5)

        # First get the source node information
        source_nodes = LSDMap_CP.FindSourceInformation(chanPointData)

        # Find all the channels with chi greater than a threshold
        selected_sources = LSDMap_CP.FindShortSourceChannels(source_nodes,source_chi_threshold)

        # Now loop through these sources and get their eaasting and northing
        for source in selected_sources:

            # Get the latitude and longitude
            latitude = source_nodes[source]["Latitude"]
            longitude = source_nodes[source]["Longitude"]

            # Convert to UTM
            [this_easting,this_northing] = LSDMap_BM.GetUTMEastingNorthing(EPSG_string,latitude,longitude)

            # convert the northing for imshow
            this_NCoord = LSDMap_BM.ConvertNorthingForImshow(FileName,this_northing)

            # now append the text
            texts2.append(ax.text(this_easting,this_NCoord, str(source),fontsize = 6, color= "b",alpha=0.7,bbox=bbox_props))


    # Adjust the basin text
    adjust_text(texts,x=buffered_east,y=buffered_north,autoalign='xy',ax=ax)

    # Adjust the sources texts
    adjust_text(texts2, arrowprops=dict(arrowstyle="->", color='r', lw=0.5),expand_points=(3, 3),
            force_points=1)


    # Now to fix up the axes
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

    # This affects all axes because we set share_all = True.
    ax.set_xlim(x_min,x_max)
    ax.set_ylim(y_max,y_min)

    # Adjust the text
    #adjust_text(texts,x=buffered_east,y=buffered_north,autoalign='xy',ax=ax)
    #adjust_text(texts)


    ax.set_xticks(xlocs)
    ax.set_yticks(ylocs)

    print("The figure format is: " + FigFormat)
    if FigFormat == 'show':
        plt.show()
    elif FigFormat == 'return':
        return fig
    else:
        plt.savefig(FigFileName,format=FigFormat,dpi=500)
        fig.clf()



#==============================================================================
# Make a simple hillshade plot
def Hillshade(raster_file, azimuth = 315, angle_altitude = 45, NoDataValue = -9999,z_factor = 1):
    """Creates a hillshade raster

    Args:
        raster_file (str): The name of the raster file with path and extension.
        azimuth (float): Azimuth of sunlight
        angle_altitude (float): Angle altitude of sun
        NoDataValue (float): The nodata value of the raster

    Returns:
        HSArray (numpy.array): The hillshade array

    Author:
        DAV and SWDG
    """

    #print("The raster file is: "+raster_file)

    # You have passed a filepath to be read in as a raster
    if isinstance(raster_file, str):
      array = LSDMap_IO.ReadRasterArrayBlocks(raster_file,raster_band=1)

    # You already have an array and just want the hill shade
    elif isinstance(raster_file, np.ndarray):
      array = raster_file
    else:
        print("raster_file must be either a filepath (string) or a numpy array. Try again.")

    # DAV attempting mask nodata vals
    nodata_mask = array == NoDataValue
    array[nodata_mask] = np.nan

    x, y = np.gradient(array)
    slope = np.pi/2. - np.arctan(np.multiply(z_factor,np.sqrt(x*x + y*y)))
    aspect = np.arctan2(-x, y)
    azimuthrad = azimuth*np.pi / 180.
    altituderad = angle_altitude*np.pi / 180.


    shaded = np.sin(altituderad) * np.sin(slope)\
     + np.cos(altituderad) * np.cos(slope)\
     * np.cos(azimuthrad - aspect)



    #this_array = 255*(shaded + 1)/2
    return 255*(shaded + 1)/2
#==============================================================================


#==============================================================================
def SwathPlot(path, filename, axis):
    """A function that creates a swath in either the x or y direction only.
       Averages across entire DEM. Exceedingly basic.

    Args:
        path (str): the path to the raster
        filename (str): the name of the file
        axis (int): if 0, swath along x-axis, if not swath along y-axis

    Returns:
        A plot of the swath

    Author:
        SMM
    """

    # get the path to the raster file
    NewPath = LSDOst.AppendSepToDirectoryPath(path)
    FileName = NewPath+filename

    # get the data vectors
    means,medians,std_deviations,twentyfifth_percentile,seventyfifth_percentile = LSDMap_BM.SimpleSwath(path, filename, axis)

    print("Means shape is: ")
    print(means.shape)

    x_vec,y_vec = LSDMap_IO.GetLocationVectors(FileName)


    print("X shape is: ")
    print(x_vec.shape)

    print("Y shape is: ")
    print(y_vec.shape)

    import matplotlib.pyplot as plt

    # Set up fonts for plots
    label_size = 20
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure, sized for a ppt slide
    fig = plt.figure(1, facecolor='white',figsize=(10,7.5))

    gs = plt.GridSpec(100,75,bottom=0.1,left=0.1,right=0.9,top=1.0)
    ax = fig.add_subplot(gs[10:100,10:75])

    if axis == 0:
        dir_vec = x_vec
    else:
        dir_vec = y_vec

    min_sd = np.subtract(means,std_deviations)
    plus_sd = np.add(means,std_deviations)

    ax.plot(dir_vec,means, linewidth = 2, color = "red")
    #ax.fill_between(dir_vec, twentyfifth_percentile, seventyfifth_percentile, facecolor='green', alpha = 0.7, interpolate=True)
    ax.fill_between(dir_vec, min_sd, plus_sd, facecolor='blue', alpha = 0.5, interpolate=True)

    ax.set_xlim(dir_vec[0],dir_vec[-1])

    plt.show()
#==============================================================================


def LongitudinalSwathAnalysisPlot(full_file_path, ax):
    """Longitudinal channel swath profiles from the swath analysis driver
        output.

    Author:
        DAV & DTM
    """

    # FileList = dir +
    f = open(full_file_path, 'r')
    lines = f.readlines()
    N_Segments = len(lines)-1

    # Initialise a bunch of vectors
    distance = np.zeros(N_Segments)
    mean = np.zeros(N_Segments)
    sd = np.zeros(N_Segments)
    minimum = np.zeros(N_Segments)
    LQ = np.zeros(N_Segments)
    median = np.zeros(N_Segments)
    UQ = np.zeros(N_Segments)
    maximum = np.zeros(N_Segments)
    for i in range (0, N_Segments):

       line = lines[i+1].strip().split(" ")
       distance[i] = float(line[0])
       mean[i] = float(line[1])
       sd[i] = float(line[2])
       minimum[i] = float(line[3])
       LQ[i] = float(line[4])
       median[i] = float(line[5])
       UQ[i] = float(line[6])
       maximum[i] = float(line[7])

       # if there is nodata (-9999) replace with the numpy nodata entry
       if (mean[i] == -9999):
         mean[i] = np.nan
         sd[i] = np.nan
         minimum[i] = np.nan
         LQ[i] = np.nan
         median[i] = np.nan
         UQ[i] = np.nan
         maximum[i] = np.nan

    f.close()

    #######################
    #                     #
    #  Plot data          #
    #                     #
    #######################
    OutputFigureFormat = 'pdf'
    # SET UP COLOURMAPS
    # PLOT 1 -> transverse profile
    #plt.figure(1, facecolor='White',figsize=[10,6])

    #plt.plot(distance,median,'-')
    #plt.plot(distance,LQ,'--',color='green')
    #plt.plot(distance,UQ,'--',color='green')
    #plt.plot(distance,minimum,'-.',color='black')
    #plt.plot(distance,maximum,'-.',color='black')

    ax.plot(distance,mean,'-', alpha=0.6)

    # Configure final plot
    ax.spines['top'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)
    ax.tick_params(axis='both', width=1)
    plt.ylabel('Average elevation difference (m)')
    plt.xlabel('Distance along channel longitudinal profile (m)')
    plt.subplots_adjust(bottom=0.15,left=0.18)

def MultiLongitudinalSwathAnalysisPlot(data_dir, wildcard_fname, maximum=0):
    """For multiple overlaid channel swath profiles.

    Arguments:
        data_dir: Full path to the directory containing the swath profile text
                    files.
        wildcard_fname: Any string that can be expanded by python's `glob` to
                        give a list of files, e.g. "some_common_prefix_*.txt"
        maximum (optional): Hacky solution, but give this a float value and it
                            will plot the 'zero' line on your swath profile.
                            C.f. maximum length of channel)

    Author:
        DAV
    """
    import glob

    fig, ax = plt.subplots()

    for f in sorted(glob.glob(data_dir + wildcard_fname)):
        print(f)
        LongitudinalSwathAnalysisPlot(f, ax)

    # Plot the zero line on the graph, if length supplied.
    x, y = function_sketcher((lambda x: x*0), np.linspace(0, maximum, 100))
    ax.plot( x, y, linestyle="--", color="k" )

    fig.canvas.draw()
    return fig, ax

def function_sketcher(formula, x_range):
    """Creates x and y values over specified range for a given lambda function.

    Arguments:
        formula: a lambda function describing a mathematical function
            example: "lambda x: x**3+2*x-4" (ommit the quotations when passing)
                    would represent the line $y = x^3 + 2*x - 4$.
        x_range: A range of numbers to calculate y from = range(0,10) or np.linspace

    Returns:
        x: the x coordinates of the function graph as a numpy array (same as x_range)
        y: the y coordinates of the function graph as a numpy array

    Author:
        DAV
    """
    x = np.array(x_range)
    y = formula(x)
    return x, y

#==============================================================================
def round_to_n(x, n):
    """A rounding function

    Args:
        x (float): The number to be rounded
        n (int): The number of digits to be rounded to.

    Returns:
        Rounded (float): The rounded number

    Author:
        SMM
    """

    if n < 1:
        raise ValueError("number of significant digits must be >= 1")
    # Use %e format to get the n most significant digits, as a string.
    format = "%." + str(n-1) + "e"
    as_string = format % x
    return float(as_string)
#==============================================================================

def init_plotting_DV():
    """Initial plotting parameters.

    Author:
        DAV
    """
    plt.rcParams['figure.figsize'] = (8, 8)
    plt.rcParams['font.size'] = 17
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['axes.labelsize'] = 1.2*plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = 1.2*plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']
#    plt.rcParams['savefig.dpi'] = 2*plt.rcParams['savefig.dpi']
    plt.rcParams['xtick.major.size'] = 3
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.major.size'] = 3
    plt.rcParams['ytick.minor.size'] = 3
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['ytick.minor.width'] = 1
    plt.rcParams['legend.frameon'] = True
    plt.rcParams['legend.loc'] = 'center left'
    plt.rcParams['axes.linewidth'] = 1
    plt.rcParams['xtick.minor.visible'] = False
    plt.rcParams['ytick.minor.visible'] = False
    plt.rcParams['lines.linewidth'] = 2



















#
