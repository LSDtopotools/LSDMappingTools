# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 16:09:29 2016

@author: NauticalMile
@stackoverflow: http://stackoverflow.com/questions/16992038/inline-labels-in-matplotlib
"""

from math import atan2 as _atan2, degrees as _degrees
import numpy as _np
import re as _re

#Label line with line2D label data
def labelLine(line,x,label=None,align=True,**kwargs):

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
    """
    Makes a string by splitting a file name. 
    """
    # Passing a list of delimiters to the re.split function
    part1 = _re.split("[_.]", fname)[0]
    part2 = _re.split("[_.]", fname)[1]
    part3 = _re.split("[_.]", fname)[2]

    part = part2 + '_' + part3
    print part
    return part
