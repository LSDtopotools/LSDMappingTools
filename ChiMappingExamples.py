# -*- coding: utf-8 -*-
"""
Created on Mon Jun 05 12:28:29 2017

@author: smudd
"""


import matplotlib

# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')


import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rcParams
from LSDMapFigure.PlottingRaster import MapFigure


if __name__ == "__main__":

    # Change these filenames and paths to suit your own files
    DataDirectory = 'C:\\VagrantBoxes\\LSDTopoTools\\Topographic_projects\\LSDTT_chi_examples\\'
    Base_file = 'Xian'


    # This is for the first example. Uncomment to get a hillshade image
    ExampleOne_SimpleHillshade(DataDirectory,Base_file)



def ExampleOne_SimpleHillshade(DataDirectory,Base_file):

    # specify the figure size and format
    size_format='ESURF'
    FigFormat = 'png'
    fig_size_inches = 12
    ax_style = "Normal"
    
    # Get the filenames you want    
    BackgroundRasterName = Base_file+"_hs.bil"    
    DrapeRasterName = Base_file+".bil"

    # clear the plot
    plt.clf() 
    
    # this is where we want the colourbar
    cbar_loc = "right"
    
    # set up the base image and the map
    MF = MapFigure(BackgroundRasterName, DataDirectory,coord_type="UTM_km",colourbar_location = cbar_loc)
    MF.add_drape_image(DrapeRasterName,DataDirectory,colourmap = "jet", alpha = 0.6)
    
    # Save the image
    ImageName = DataDirectory+"Xian_example1.png" 
    MF.save_fig(fig_width_inches = fig_size_inches, FigFileName = ImageName, axis_style = ax_style, Fig_dpi = 250)








