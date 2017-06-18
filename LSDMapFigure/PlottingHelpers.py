# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 10:29:18 2017

@author: smudd
"""

#==============================================================================
def MapFigureSizer(figure_width_inches,aspect_ratio, cbar_loc = "None",
                   cbar_width = 0.2, 
                   cbar_text_width = 0.4,
                   cbar_padding = 0.1,
                   cbar_fraction = 1,
                   whitespace_padding = 0.1,
                   map_text_width = 0.65,
                   map_text_height = 0.45):
    """This function takes a string size argument and calculates the size of the
    various components of a plot based on a map image making up the centre of the
    figure.

    We use inches because bloody yanks wrote matplotlib and figures in matplotlib use inches.
    Luckily we do not have to calculate bushels or furlongs of pixels somewhere,
    and the inches stupidity is somewhat limited.

    Args:
        figure_width_inches (flt): The figure width in inches
        aspect_ratio (flt): The width to height ratio of the data
        cbar_loc (string): the location of the colourbar, either "left", "right", "top",
        "bottom", or "none"
        cbar_width (flt): the width of the colorbar
        text_padding (list): the padding around the map from the whitespace and the axis/tick labels


    """


    # The text padding list holds the height and width, in inches of
    # [0] = left/right tick marks+tick labels
    # [1] = left/right text (i.e., and exis label)
    # [2] = top/bottom tick marks+tick labels
    # [3] = top/bottom text (i.e., and exis label)
    # [4] = whitespace at the edge of the figure + between map and colourbar
    # [5] = colorbar text (e.g. tick labels+axis label)


    # This gets returned, we use it to make the figure.

    #======================================================

    # Now we need to figure out where the axis are. Sadly this requires
    # a load of tedious conditional statments about the location of the axes

    #Will's changes:
    #1/ Make space for the colourbar label on the left
    #2/ Double the cbar_padding to leave space for the colourbar values on the right.
    # NB: this should later be a function of font size
    #3/ Changed rotation of colourbar text to 90 and the labelpad to -75 in PlottingRaster.py

    if cbar_loc == "left":
        cbar_left_inches = whitespace_padding + cbar_text_width
        #cbar_left_inches = whitespace_padding
        map_left_inches = cbar_left_inches+cbar_width+map_text_width + 2*cbar_padding
        #map_left_inches = cbar_left_inches+cbar_width+cbar_text_width+map_text_width+cbar_padding
        map_width_inches = figure_width_inches-map_left_inches-whitespace_padding
        map_height_inches = map_width_inches/aspect_ratio



        map_bottom_inches = whitespace_padding+map_text_height
        cbar_bottom_inches = map_bottom_inches
        figure_height_inches = map_bottom_inches+map_height_inches+whitespace_padding

        fig_size_inches = [figure_width_inches,figure_height_inches]

        print("cbar_left: "+str(cbar_left_inches)+" map left: "+str(map_left_inches))
        print("cbar_bottom: "+str(cbar_bottom_inches)+ " map bottom: "+str(map_bottom_inches))

        map_axes = [map_left_inches/figure_width_inches,
                    map_bottom_inches/figure_height_inches,
                    map_width_inches/figure_width_inches,
                    map_height_inches/figure_height_inches]
        cbar_axes = [cbar_left_inches/figure_width_inches,
                    map_bottom_inches/figure_height_inches,
                    cbar_width/figure_width_inches,
                    cbar_fraction*(map_height_inches/figure_height_inches)]

    elif cbar_loc == "right":

        map_left_inches = whitespace_padding+map_text_width
        cbar_left_inches= figure_width_inches-whitespace_padding-cbar_width-cbar_text_width
        map_right_inches = cbar_left_inches-cbar_padding

        map_width_inches = map_right_inches-map_left_inches
        map_height_inches = map_width_inches/aspect_ratio

        map_bottom_inches = whitespace_padding+map_text_height
        cbar_bottom_inches = map_bottom_inches
        figure_height_inches = map_bottom_inches+map_height_inches+whitespace_padding

        fig_size_inches = [figure_width_inches,figure_height_inches]

        print("cbar_left: "+str(cbar_left_inches)+" map left: "+str(map_left_inches))
        print("cbar_bottom: "+str(cbar_bottom_inches)+ " map bottom: "+str(map_bottom_inches))


        map_axes = [map_left_inches/figure_width_inches,
                    map_bottom_inches/figure_height_inches,
                    map_width_inches/figure_width_inches,
                    map_height_inches/figure_height_inches]
        cbar_axes = [cbar_left_inches/figure_width_inches,
                    map_bottom_inches/figure_height_inches,
                    cbar_width/figure_width_inches,
                    cbar_fraction*(map_height_inches/figure_height_inches)]

    elif cbar_loc == "top":
        print("I am placing the colourbar on the top")

        map_left_inches = whitespace_padding+map_text_width
        map_right_inches = figure_width_inches-whitespace_padding
        map_width_inches = map_right_inches-map_left_inches
        map_height_inches = map_width_inches/aspect_ratio

        cbar_left_inches= map_left_inches

        map_bottom_inches = whitespace_padding+map_text_height
        cbar_bottom_inches = map_bottom_inches+map_height_inches+cbar_padding+cbar_text_width

        figure_height_inches = cbar_bottom_inches+cbar_width+whitespace_padding

        fig_size_inches = [figure_width_inches,figure_height_inches]

        print("cbar_left: "+str(cbar_left_inches)+" map left: "+str(map_left_inches))
        print("cbar_bottom: "+str(cbar_bottom_inches)+ " map bottom: "+str(map_bottom_inches))


        map_axes = [map_left_inches/figure_width_inches,
                    map_bottom_inches/figure_height_inches,
                    map_width_inches/figure_width_inches,
                    map_height_inches/figure_height_inches]
        cbar_axes = [cbar_left_inches/figure_width_inches,
                    cbar_bottom_inches/figure_height_inches,
                    cbar_fraction*(map_width_inches/figure_width_inches),
                    cbar_width/figure_height_inches]

    elif cbar_loc == "bottom":
        print("I am placing the colourbar on the bottom")

        map_left_inches = whitespace_padding+map_text_width
        map_right_inches = figure_width_inches-whitespace_padding
        map_width_inches = map_right_inches-map_left_inches
        map_height_inches = map_width_inches/aspect_ratio

        cbar_left_inches= map_left_inches

        cbar_bottom_inches = whitespace_padding+cbar_text_width
        map_bottom_inches = cbar_bottom_inches+cbar_width+cbar_padding+map_text_height

        whitespace_padding+map_text_height

        figure_height_inches = map_bottom_inches+map_height_inches+whitespace_padding

        fig_size_inches = [figure_width_inches,figure_height_inches]

        print("cbar_left: "+str(cbar_left_inches)+" map left: "+str(map_left_inches))
        print("cbar_bottom: "+str(cbar_bottom_inches)+ " map bottom: "+str(map_bottom_inches))


        map_axes = [map_left_inches/figure_width_inches,
                    map_bottom_inches/figure_height_inches,
                    map_width_inches/figure_width_inches,
                    map_height_inches/figure_height_inches]
        cbar_axes = [cbar_left_inches/figure_width_inches,
                    cbar_bottom_inches/figure_height_inches,
                    cbar_fraction*(map_width_inches/figure_width_inches),
                    cbar_width/figure_height_inches]

    else:
        print("No colourbar")

        map_left_inches = whitespace_padding+map_text_width
        map_right_inches = figure_width_inches-whitespace_padding
        map_width_inches = map_right_inches-map_left_inches
        map_height_inches = map_width_inches/aspect_ratio

        map_bottom_inches = whitespace_padding+map_text_height

        figure_height_inches = map_bottom_inches+map_height_inches+whitespace_padding

        fig_size_inches = [figure_width_inches,figure_height_inches]

        map_axes = [map_left_inches/figure_width_inches,
                    map_bottom_inches/figure_height_inches,
                    map_width_inches/figure_width_inches,
                    map_height_inches/figure_height_inches]
        cbar_axes = None

    print("The figure size is: ")
    print fig_size_inches
    print("Map axes are:")
    print(map_axes)
    print("cbar_axes are:")
    print(cbar_axes)
    return fig_size_inches, map_axes, cbar_axes


#==============================================================================
# Formats ticks for an imshow plot in UTM
# Filename is the name of the file with full path
# x_max, x_min, y_max, y_min are the extent of the plotting area (NOT the DEM)
# n_target ticks are the number of ticks for plotting
#------------------------------------------------------------------------------
def GetTicksForUTMNoInversion(FileName,x_max,x_min,y_max,y_min,n_target_tics):
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
