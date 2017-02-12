# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 10:29:18 2017

@author: smudd
"""

#==============================================================================
def MapFigureSizer(figure_width_inches,aspect_ratio, cbar_loc = "None"):
    """This function takes a string size argument and calculates the size of the
    various components of a plot based on a map image making up the centre of the
    figure. 
    
    We use inches because bloody yanks wrote matplotlib and figures in matplotlib use inches.
    Luckily we do not have to calculate bushels or furlongs of pixels somewhere, 
    and the inches stupidity is somewhat limited. 
    
    Args:
        figure_width_inches (flt): The figure width in inches
        aspect_ratio (flt): The width to height ratio of the data
    
    """
    
    # By default the axes labels with tick label take up ~0.9 inches.
    # later we will adjust this with the size of the font but for now give a default
    cumulative_label_width = 0.9    
    cumulative_label_height = 0.9    
    
    # first check if cbar is on the left or right
    if cbar_loc == "left" or cbar_loc == "right":
        #The colourbar takes up ~1 inch of space
        cumulative_label_width = cumulative_label_width+1.5
    elif cbar_loc == "top" or cbar_loc == "bottom":
        cumulative_label_height = cumulative_label_height+1.5
        
    # Now get the width of the map
    map_width_inches = figure_width_inches-cumulative_label_width
    print("Map width is: "+str(map_width_inches))
    
    # now calculate the height of map and then the figure
    map_height_inches = map_width_inches/aspect_ratio
    print("Map height is: "+str(map_height_inches))
    
    # now get the full height of the figure
    figure_height_inches = cumulative_label_height+map_height_inches
    
    # This gets returned, we use it to make the figure. 
    fig_size_inches = [figure_width_inches,figure_height_inches]

    # Now we need to figure out where the axis are. Sadly this requires
    # a load of tedious conditional statments about the location of the axes
    if cbar_loc == "left":
        cbar_left_inches = 0.9
        cbar_bottom_inches = 0.9
        map_left_inches = 2.4
        mab_bottom_inches = 0.9
        
        map_axes = [map_left_inches/figure_width_inches,
                    mab_bottom_inches/figure_height_inches,
                    map_width_inches/figure_width_inches,
                    map_height_inches/figure_height_inches]        
        cbar_axes = [cbar_left_inches/figure_width_inches,
                    cbar_bottom_inches/figure_height_inches,
                    0.6/figure_width_inches,
                    map_height_inches/map_height_inches]         
        
    elif cbar_loc == "right":
        cbar_left_inches = 0.9
        cbar_bottom_inches = 0.9
        map_left_inches = 1.5
        mab_bottom_inches = 0.9        

        map_axes = [map_left_inches/figure_width_inches,
                    mab_bottom_inches/figure_height_inches,
                    map_width_inches/figure_width_inches,
                    map_height_inches/figure_height_inches]        
        cbar_axes = [cbar_left_inches/figure_width_inches,
                    cbar_bottom_inches/figure_height_inches,
                    0.6/figure_width_inches,
                    map_height_inches/map_height_inches] 
    
    elif cbar_loc == "top":
        print("I am placing the colourbar on the top")        
        
        cbar_left_inches = 0.9
        mab_bottom_inches = 0.9
        map_left_inches = 0.9
        cbar_bottom_inches = 0.9+map_height_inches+0.9
        
        
        
        

        map_axes = [map_left_inches/figure_width_inches,
                    mab_bottom_inches/figure_height_inches,
                    map_width_inches/figure_width_inches,
                    map_height_inches/figure_height_inches]        
        cbar_axes = [cbar_left_inches/figure_width_inches,
                    cbar_bottom_inches/figure_height_inches,
                    map_width_inches/figure_width_inches,
                    0.4/map_height_inches] 
        
    elif cbar_loc == "bottom":
        cbar_left_inches = 0.9
        cbar_bottom_inches = 0.9
        map_left_inches = 0.9
        mab_bottom_inches = 0.9

        map_axes = [map_left_inches/figure_width_inches,
                    mab_bottom_inches/figure_height_inches,
                    map_width_inches/figure_width_inches,
                    map_height_inches/figure_height_inches]        
        cbar_axes = [cbar_left_inches/figure_width_inches,
                    cbar_bottom_inches/figure_height_inches,
                    0.6/figure_width_inches,
                    map_height_inches/map_height_inches] 
        
    else:
        map_left_inches = 0.9
        mab_bottom_inches = 0.9       

        map_axes = [map_left_inches/figure_width_inches,
                    mab_bottom_inches/figure_height_inches,
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