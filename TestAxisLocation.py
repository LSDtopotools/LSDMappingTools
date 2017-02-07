# -*- coding: utf-8 -*-
"""
Created on Mon Feb 06 08:42:57 2017

@author: smudd
"""

import LSDPlottingTools as LSDP
import matplotlib.pyplot as plt
from matplotlib import rcParams

def TestAxisLocation():

    tick_label_size = 10
    text_size = 12

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = text_size 
    rcParams['xtick.labelsize'] = tick_label_size
    rcParams['ytick.labelsize'] = tick_label_size 
    
    #DataDirectory = "C:\\Vagrantboxes\\LSDTopoTools\\Topographic_projects\\Meghalaya\\Divides\\"    
    DataDirectory = "T:\\analysis_for_papers\\Meghalaya\\divide_migration\\"
    Base_file = "Mega_divide"
     
    bil = ".bil"
    
    #Filename = "Mega_clip.bil"
    #HSFilename = "Mega_clip_hs.bil"
    #BasinFilename = "Mega_clip_AllBasins.bil"
    
    DEMname = DataDirectory+Base_file+bil 
    FigFileName = DataDirectory+Base_file+"Picture.png"
    FigFormat = "png"

    # get the data
    #raster = LSDP.ReadRasterArrayBlocks(DEMname)
    

    # now get the extent
    extent_raster = LSDP.GetRasterExtent(DEMname)
    
    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]
    
    print(extent_raster)
    
    fig = plt.figure(1, facecolor='white',figsize=(10,10))
    
    # Add an axis. This will be used to check how high the text is. 
    ax2 = fig.add_axes([0.1,0.1,0.7,0.7],zorder = -1)
    
    # turn off the ticks
    ax2.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off

    ax2.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='off') # labels along the bottom edge are off
    
    ax1 = fig.add_axes([0.0,0.09,0.85,0.85])
    #ax1.set(alpha=0.25)    
    ax1.text(0.1,0.1,
             "x: Tick 10pt, Font 12pt need 0.9 inches.\n",
             transform=ax1.transAxes,)

    #im = ax1.imshow(raster[::-1], "jet", extent = extent_raster, interpolation="nearest")

    
    n_target_tics = 5
    xlocs,ylocs,new_x_labels,new_y_labels = LSDP.GetTicksForUTM(DEMname,x_max,x_min,y_max,y_min,n_target_tics)  

    ax1.set_xticks(xlocs)
    ax1.set_yticks(ylocs)      
    ax1.set_xticklabels(new_x_labels,rotation=60)
    ax1.set_yticklabels(new_y_labels) 
    
    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap        
    ax1.tick_params(axis='both', width=1, pad = 2)
    for tick in ax1.xaxis.get_major_ticks():
        tick.set_pad(2)  

    

    # This affects all axes because we set share_all = True.
    ax1.set_xlim(x_min,x_max)    
    ax1.set_ylim(y_max,y_min)
    ax1.set_xlabel("YumYumDonuts")
    ax1.set_ylabel("Bigmama")
    
    #ax3 = fig.add_axes([0.0,0.5,0.5,0.5])
    #ax1.text(0,0,"Hello",fontsize = 94)    
    plt.savefig(FigFileName,format=FigFormat,dpi=100)
    
    yo= 284/3
    print(yo)
    
    # So in every inch there are 94 points of text
    # The bottom 
    
    plt.show()
    
    
    

if __name__ == "__main__":
    TestAxisLocation()