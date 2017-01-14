## LSDMatplotlibExtensions.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## These functions are wrappers and extensions around the matplotlib 
## library to provide additionaly plotting functionality without having
## to copy and paste code all the time.
## They are not necessarily specific to LSDTopoTools usage.
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## DAV
## 11/01/2017
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
"""
Created on Thu Jan 12 14:33:21 2017

@author: DAV
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import re

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=-1):
    """
    Truncates a standard matplotlib colourmap so
    that you can use part of the colourange in your plots.
    Handy when the colourmap you like has very light values at
    one end of the map.
    
    Usage example:
       minColor = 0.00
       maxColor = 0.85
       inferno_t = truncate_colormap(plt.get_cmap("inferno"), minColor, maxColor) 
    """
    cmap = plt.get_cmap(cmap)
    
    if n == -1:
        n = cmap.N
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
         'trunc({name},{a:.2f},{b:.2f})'.format(name=cmap.name, a=minval, b=maxval),
         cmap(np.linspace(minval, maxval, n)))
    return new_cmap
    
def discrete_colourmap(N, base_cmap=None):
    """
    Create an N-bin discrete colourmap from the specified input colormap.
    github.com/jakevdp
    
    DAV: modified so you can pass in the string name of a colourmap
    or a Colormap object.
    """

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:
    if isinstance(base_cmap, mcolors.Colormap):
        base = base_cmap
    elif isinstance(base_cmap, str):
        base = plt.cm.get_cmap(base_cmap)
    else:
        print "DrapeName supplied is of type: ", type(base_cmap)
        raise ValueError('DrapeName must either be a string name of a colormap, \
                         or a Colormap. Please try again.')
        
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

    
def make_line_label(fname):
    """
    Makes a string by splitting a file name. 
    """
    # Passing a list of delimiters to the re.split function
    part1 = re.split("[_.]", fname)[0]
    part2 = re.split("[_.]", fname)[1]
    part3 = re.split("[_.]", fname)[2]

    part = part2 + '_' + part3
    print part
    return part