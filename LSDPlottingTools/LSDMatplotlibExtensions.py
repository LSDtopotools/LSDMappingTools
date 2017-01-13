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