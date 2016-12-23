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
    print "Number of files = ", n_files
    
    # Now make the subplots
    fig, ax = pp.subplots(2,3, figsize=(cm2inch(15),cm2inch(11)))
    ax = ax.ravel()
       
    #use seaborn to get a nice color palette
    cmap_oranges = sns.light_palette("#ff8f66", input="hex", as_cmap=True, reverse=True)  
    cmap_ice = sns.light_palette("#00ffff", input="hex", as_cmap=True, reverse=True)
    
    alphabet = list(string.ascii_lowercase)
    
    for i in range (n_files):
        print "The floodplain file name is: ", FPFiles[i]
        
        # get the name of the field site
        fname = FPFiles[i].split('/')
        print fname
        split_fname = fname[-1].split('_')
        print split_fname
        HSFile = DataDirectory+split_fname[1]+'_'+split_fname[2]+"_HS.bil"
        print HSFile
         
        hillshade_raster = LSDMap_IO.ReadRasterArrayBlocks(HSFile)
        FP_raster = LSDMap_IO.ReadRasterArrayBlocks(FPFiles[i])
        FP_raster = np.ma.masked_where(FP_raster <= 0, FP_raster)
        
        # now get the extent
        extent_raster = LSDMap_IO.GetRasterExtent(HSFile)
    
        # get DEM info
        CellSize,XMin,XMax,YMin,YMax = LSDMap_IO.GetUTMMaxMin(HSFile)
        print YMin, YMax
        
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
    
def flood_maps_with_shapefile(DataDirectory):
    """
    Make subplots showing the difference between the mapped and predicted
    floodplain initiation points.
    Uses Fiona (yaaay) to read in the shapefile
    
    """
    from fiona import collection
    
    # Set up fonts
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = 8
    
    # Read in the files for each site
    HSFiles = sorted(glob(DataDirectory+'*_HS.bil'), key=str)
    FPFiles = sorted(glob(DataDirectory+'*_FP*.bil'), key=str)
    SHPFiles = sorted(glob(DataDirectory+'*_FIPs.shp'), key=str)
      
    n_files = len(FPFiles)
    print "Number of files = ", n_files
    
    # Now make the subplots
    fig, ax = pp.subplots(1,2, figsize=(cm2inch(12),cm2inch(7)))
    ax = ax.ravel()
    
    #get a list with the figure letterings
    figure_letter = ["a", "b"]
    
    for i in range (n_files):
        
        print "The hillshade file name is: ", HSFiles[i]
        print "The floodplain file name is: ", FPFiles[i]
        print "The shapefile name is: ", SHPFiles[i]
        
        # get the name of the field site
        split_fname = FPFiles[i].split('_FP')
        print split_fname[0]
         
        hillshade_raster = LSDMap_IO.ReadRasterArrayBlocks(HSFiles[i])
        FP_raster = LSDMap_IO.ReadRasterArrayBlocks(FPFiles[i])
        FP_raster = np.ma.masked_where(FP_raster <= 0, FP_raster)
        
        # now get the extent
        extent_raster = LSDMap_IO.GetRasterExtent(HSFiles[i])

        # get DEM info
        CellSize,XMin,XMax,YMin,YMax = LSDMap_IO.GetUTMMaxMin(HSFiles[i])
        print YMin, YMax
        
        #plot the rasters
        ax[i].imshow(hillshade_raster, extent = extent_raster, cmap=cmx.gray)
        ax[i].imshow(FP_raster, extent = extent_raster, cmap=cmx.bwr, alpha=0.6)
        
        #plot the mapped points
        with collection(SHPFiles[i],'r') as input:
            for point in input:
                x = point['geometry']['coordinates'][0]
                y = point['geometry']['coordinates'][1]
                ax[i].scatter(x,y, c="red", s=15)
     
        ax[i].text(0.03,0.97, figure_letter[i], bbox=dict(facecolor='white', edgecolor='k', pad=3), horizontalalignment='left', verticalalignment='top', fontsize = 8, transform=ax[i].transAxes)
        
        #change ticks
        xlocs = ax[i].xaxis.get_ticklocs()
        ylocs = ax[i].yaxis.get_ticklocs()
    
        n_target_tics = 7
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
  
    # Save figure
    pp.tight_layout(pad=0.5)
    OutputFigureName = "Comparison_with_mapped_FIPs"
    OutputFigureFormat = 'pdf'
    pp.savefig(DataDirectory+OutputFigureName + '.' + OutputFigureFormat, format=OutputFigureFormat, dpi=300)





