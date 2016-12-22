## LSDMap_Subplots.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## These functions are tools to deal with creating nice subplots from multiple 
## rasters
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## FJC
## 22/12/2016
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

from glob import glob
import numpy as np
import matplotlib.pyplot as pp
import string
import matplotlib.image as mpimg
import matplotlib.cm as cmx
from matplotlib import rcParams
import LSDMap_GDALIO as LSDMap_IO
import LSDMap_BasicManipulation as LSDMap_BM
import LSDOSystemTools as LSDOst
import LSDMap_BasicPlotting as LSDMap_BP

#==============================================================================
# Convert cm to inch for figure sizing
#------------------------------------------------------------------------------
def cm2inch(value):
    return value/2.54
    
#==============================================================================
# Function to create nice field sites figure from all the hillshades in a folder
# N_HSFiles = number of field sites
# Also reads in a series of map files to show locations of each site. At the moment
# you have to manually indicate on these maps where the field site is - would
# be nice to automate this when I have some time to mess around with it.
#------------------------------------------------------------------------------

def field_sites(DataDirectory, N_HSFiles, NRows, NCols, n_target_ticks):
    
    # Set up fonts
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = 6
    
    # Read in the files for each site
    HSFiles = sorted(glob(DataDirectory+'*_HS.bil'), key=str)
    MapFiles = sorted(glob(DataDirectory+'*_map.eps'), key=str)
    print MapFiles
      
    n_files = len(HSFiles)+len(MapFiles)
    print "Number of files = ", n_files
    
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
            print "The hillshade file name is: ", HSFiles[i]
             
            hillshade_raster = LSDMap_IO.ReadRasterArrayBlocks(HSFiles[i])
            hillshade_raster = np.ma.masked_where(hillshade_raster == -9999, hillshade_raster)
            
            # now get the extent
            extent_raster = LSDMap_IO.GetRasterExtent(HSFiles[i])
        
            # get DEM info
            CellSize,XMin,XMax,YMin,YMax = LSDMap_IO.GetUTMMaxMin(HSFiles[i])
            print YMin, YMax
            
            #plot the rasters
            ax[i].imshow(hillshade_raster, extent = extent_raster, cmap=cmx.gray)
            ax[i].text(0.05,0.97, alphabet[i], horizontalalignment='left', verticalalignment='top', bbox=dict(facecolor='white', edgecolor='k', pad=3), fontsize = 8, transform=ax[i].transAxes)
    
            #change ticks
            xlocs = ax[i].xaxis.get_ticklocs()
            ylocs = ax[i].yaxis.get_ticklocs()
        
            new_xlocs,new_ylocs,new_x_labels,new_y_labels = LSDMap_BP.GetTicksForUTM(HSFiles[i],xlocs.max(),xlocs.min(),ylocs.max(),ylocs.min(),n_target_ticks)
        
            # change the location of the ticks depending on subplot placement
            if i < 2:
                ax[i].xaxis.tick_top()
            if i % 2 != 0:
                ax[i].yaxis.tick_right()
            
            ax[i].set_xticklabels(new_x_labels, rotation=30)
            ax[i].set_yticklabels(new_y_labels)
            ax[i].tick_params(axis='x', pad=7)
            ax[i].tick_params(axis='y', pad=7)
        
        if i >= N_HSFiles:
            print "The map file name is: ", MapFiles[file_counter]
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
        
    pp.xlabel('Easting (m)', fontsize=8, labelpad=-122)
    pp.ylabel('Northing (m)', fontsize=8, labelpad=10, position=(0.0,0.67))
    pp.tight_layout(pad=0.1, w_pad = 0.1, h_pad = 0.2)
    OutputFigureName = "field_sites"
    OutputFigureFormat = 'pdf'
    pp.savefig(DataDirectory+OutputFigureName + '.' + OutputFigureFormat, format=OutputFigureFormat, dpi=300)
    #pp.show()





