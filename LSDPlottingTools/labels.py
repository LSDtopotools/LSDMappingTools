# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 16:09:29 2016

A series of functions to provide extra functionality to matplotlib 
involving the creation of labels for plots.

    Author: DAV
@stackoverflow: http://stackoverflow.com/questions/16992038/inline-labels-in-matplotlib
"""
from __future__ import absolute_import, division, print_function, unicode_literals

from math import atan2 as _atan2, degrees as _degrees
import numpy as _np
import re as _re

#Label line with line2D label data
def labelLine(line,x,label=None,align=True,**kwargs):
    """Places a label on a line plot and orients it to run parallel with the line.

    Given a matplotlib Line instance and x-coordinate, places `label` at the x-coord
    on the given line and orientates it parallel to the line. 

    Author: http://stackoverflow.com/questions/16992038/inline-labels-in-matplotlib

    Arguments:
        line: Matplotlib Line instance
        x: x-coordinate on the line at which the label will be positioned.
        label (str): The label to be added to the line.
        align (bool): whether or not to align the label parallel to the line

    """

    ax = line.get_axes()
    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if (x < xdata[0]) or (x > xdata[-1]):
        print('x label location is outside data range!')
        return

    #Find corresponding y co-ordinate and angle of the
    ip = 1
    for i in range(len(xdata)):
        if x < xdata[i]:
            ip = i
            break

    y = ydata[ip-1] + (ydata[ip]-ydata[ip-1])*(x-xdata[ip-1])/(xdata[ip]-xdata[ip-1])

    if not label:
        label = line.get_label()

    if align:
        #Compute the slope
        dx = xdata[ip] - xdata[ip-1]
        dy = ydata[ip] - ydata[ip-1]
        ang = _degrees(_atan2(dy,dx))

        #Transform to screen co-ordinates
        pt = _np.array([x,y]).reshape((1,2))
        trans_angle = ax.transData.transform_angles(_np.array((ang,)),pt)[0]

    else:
        trans_angle = 0

    #Set a bunch of keyword arguments
    if 'color' not in kwargs:
        kwargs['color'] = line.get_color()

    if ('horizontalalignment' not in kwargs) and ('ha' not in kwargs):
        kwargs['ha'] = 'center'

    if ('verticalalignment' not in kwargs) and ('va' not in kwargs):
        kwargs['va'] = 'center'

    if 'backgroundcolor' not in kwargs:
        kwargs['backgroundcolor'] = ax.get_axis_bgcolor()

    if 'clip_on' not in kwargs:
        kwargs['clip_on'] = True

    if 'zorder' not in kwargs:
        kwargs['zorder'] = 2.5

    ax.text(x,y,label,rotation=trans_angle,**kwargs)

def labelLines(lines,align=True,xvals=None,**kwargs):
    """Version of labelLine that assigns labels for all lines
       in a plot.

    Similar to labelLine, except a list of lines is passed.

    Argumnets: 
        lines (list): A list of the lines to be labeled.
        xvals: A list of x-coordinates where the labels should be anchored.
    """

       

    ax = lines[0].get_axes()
    labLines = []
    labels = []

    #Take only the lines which have labels other than the default ones
    for line in lines:
        label = line.get_label()
        if "_line" not in label:
            labLines.append(line)
            labels.append(label)

    if xvals is None:
        xmin,xmax = ax.get_xlim()
        xvals = _np.linspace(xmin,xmax,len(labLines)+2)[1:-1]

    for line,x,label in zip(labLines,xvals,labels):
        labelLine(line,x,label,align,**kwargs)

def make_line_label(fname):
    """Makes a string (label) by splitting a file name. 

    Warning:
       A lot of this is hard coded to split according to certain filenaming
       conventions, separated by underscored. e.g. MyFile_part1_part2_part3.file
       So you should modify this to fit your own file naming
       convention. 

    Todo: Rewrite this as a more generic function.

    Arguments:
        fname (str): Filename to create labels from.

    Author: DAV

    """
    parts = []
    # Passing a list of delimiters to the re.split function
    try:
        part1 = _re.split("[_.]", fname)[0]
        print(part1)
        parts.append(part1)
    except:
        print("Unable to extract file name parts as strings by splitting \
              after first underscore")
        part1 = None
        
    try:
        part2 = _re.split("[_.]", fname)[1]
        print(part2)
        parts.append(part2)
    except:
        print("Unable to extract file name parts as strings by splitting \
              after second underscore")     
        part2 = None
        
    try:
        part3 = _re.split("[_.]", fname)[2]
        print(part3)
        parts.append(part3)
    except:
        print("Unable to extract file name parts as strings by splitting \
              after third underscore")
        part3 = None
    
    print(parts)  
    label = '_'.join(parts[1:])

    print("FIGURE LABEL IS:", label)
    return label
