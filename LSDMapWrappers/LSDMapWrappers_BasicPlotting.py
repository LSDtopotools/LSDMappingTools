"""
    This contains wrapper functions that simplify plotting raster
    and vector data for publication-ready figures. 
    
    The documentation of the examples can be found here:
    https://lsdtopotools.github.io/LSDTopoTools_ChiMudd2014/

    Simon Mudd and Fiona Clubb, June 2017

    Released under GPL3


"""

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')



import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rcParams


"""
    IMPORTANT: You must call this function from a lower level driectory
    where both LSDPlottingTools and LSDMapFigure are in the python path!

    That is, it will not work if you call it from outside the directory structure.

"""
import LSDPlottingTools as LSDP
from LSDMapFigure.PlottingRaster import MapFigure
import LSDMapFigure.PlottingHelpers as PlotHelp
#import LSDPlottingTools.LSDMap_VectorTools as LSDMap_VT


def SimpleHillshade(DataDirectory,Base_file, cmap = "jet", cbar_loc = "right", size_format = "ESURF", fig_format = "png", dpi = 250):
    """
    This function makes a shaded relief plot of the DEM with the basins coloured
    by the basin ID.

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        Base_file (str): The prefix for the m/n csv files
        cmap (str or colourmap): The colourmap to use for the plot
        cbar_loc (str): where you want the colourbar. Options are none, left, right, top and botton. The colourbar will be of the elevation.
                        If you want only a hillshade set to none and the cmap to "gray"
        size_format (str): Either geomorphology or big. Anything else gets you a 4.9 inch wide figure (standard ESURF size)
        fig_format (str): An image format. png, pdf, eps, svg all valid
        dpi (int): The dots per inch of the figure

    Returns:
        Shaded relief plot. The elevation is also included in the plot. 

    Author: FJC, SMM
    """
    # specify the figure size and format
    # set figure sizes based on format
    if size_format == "geomorphology":
        fig_size_inches = 6.25
    elif size_format == "big":
        fig_size_inches = 16
    else:
        fig_size_inches = 4.92126
    ax_style = "Normal"

    # Get the filenames you want
    BackgroundRasterName = Base_file+"_hs.bil"
    DrapeRasterName = Base_file+".bil"

    # clear the plot
    plt.clf()

    # set up the base image and the map
    MF = MapFigure(BackgroundRasterName, DataDirectory,coord_type="UTM_km",colourbar_location = cbar_loc)
    MF.add_drape_image(DrapeRasterName,DataDirectory,colourmap = cmap, alpha = 0.6, colorbarlabel = "Elevation (m)")

    # Save the image
    ImageName = DataDirectory+Base_file+"_hillshade."+fig_format
    MF.save_fig(fig_width_inches = fig_size_inches, FigFileName = ImageName, axis_style = ax_style, FigFormat=fig_format, Fig_dpi = dpi)






def PrintBasins(DataDirectory,fname_prefix, add_basin_labels = True, cmap = "jet", cbar_loc = "right", size_format = "ESURF", fig_format = "png", dpi = 250):
    """
    This function makes a shaded relief plot of the DEM with the basins coloured
    by the basin ID.

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        add_basin_labels (bool): If true, label the basins with text. Otherwise use a colourbar. 
        cmap (str or colourmap): The colourmap to use for the plot
        cbar_lox (str): where you want the colourbar. Options are none, left, right, top and botton. The colourbar will be of the elevation.
                        If you want only a hillshade set to none and the cmap to "gray"
        size_format (str): Either geomorphology or big. Anything else gets you a 4.9 inch wide figure (standard ESURF size)
        fig_format (str): An image format. png, pdf, eps, svg all valid
        dpi (int): The dots per inch of the figure
        
        
    Returns:
        Shaded relief plot with the basins coloured by basin ID. Uses a colourbar to show each basin

    Author: FJC, SMM
    """
    #import modules
    from LSDMapFigure.PlottingRaster import MapFigure

    # set figure sizes based on format
    if size_format == "geomorphology":
        fig_width_inches = 6.25
    elif size_format == "big":
        fig_width_inches = 16
    else:
        fig_width_inches = 4.92126

    # get the basin IDs to make a discrete colourmap for each ID
    BasinInfoDF = PlotHelp.ReadBasinInfoCSV(DataDirectory, fname_prefix)

    basin_keys = list(BasinInfoDF['basin_key'])
    basin_keys = [int(x) for x in basin_keys]

    basin_junctions = list(BasinInfoDF['outlet_junction'])
    basin_junctions = [float(x) for x in basin_junctions]

    print ('Basin keys are: ')
    print basin_keys

    # going to make the basin plots - need to have bil extensions.
    print("I'm going to make the basin plots. Your topographic data must be in ENVI bil format or I'll break!!")

    # get the rasters
    raster_ext = '.bil'
    #BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext
    BasinsName = fname_prefix+'_AllBasins'+raster_ext
    print (BasinsName)
    Basins = LSDP.GetBasinOutlines(DataDirectory, BasinsName)
    
    # If wanted, add the labels
    if add_basin_labels:
        print("I am going to add basin labels, there will be no colourbar.")
        MF = MapFigure(HillshadeName, DataDirectory,coord_type="UTM_km", colourbar_location="None")
        MF.plot_polygon_outlines(Basins, linewidth=0.8)
        MF.add_drape_image(BasinsName, DataDirectory, colourmap = cmap, alpha = 0.8, colorbarlabel='Basin ID', discrete_cmap=True, n_colours=len(basin_keys), show_colourbar = False, modify_raster_values=True, old_values=basin_junctions, new_values=basin_keys, cbar_type = int)

        # This is used to label the basins
        label_dict = dict(zip(basin_junctions,basin_keys))
        # this dict has the basin junction as the key and the basin_key as the value
    
        Points = LSDP.GetPointWithinBasins(DataDirectory, BasinsName)
        MF.add_text_annotation_from_shapely_points(Points, text_colour='k', label_dict=label_dict)
    else:
        print("I am showing the basins without text labels.")
        MF = MapFigure(HillshadeName, DataDirectory,coord_type="UTM_km", colourbar_location=cbar_loc)
        MF.plot_polygon_outlines(Basins, linewidth=0.8)
        MF.add_drape_image(BasinsName, DataDirectory, colourmap = cmap, alpha = 0.8, colorbarlabel='Basin ID', discrete_cmap=True, n_colours=len(basin_keys), show_colourbar = True, modify_raster_values=True, old_values=basin_junctions, new_values=basin_keys, cbar_type = int)

    ImageName = DataDirectory+fname_prefix+'_coloured_basins.'+fig_format
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=fig_format, Fig_dpi = dpi) # Save the figure




def MaskBasins(DataDirectory, fname_prefix, Basins_to_mask= []):
    """
    This function makes a shaded relief plot of the DEM with the basins coloured
    by the basin ID. It shows how to mask certain basins.
    Search the function for "Basins_to_mask".

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files

    Returns:
        Shaded relief plot with the basins coloured by basin ID

    Author: SMM
    """
    import numpy as np

    FigFormat = "png"
    size_format = "geomorphology"

    #import modules
    # from LSDMapFigure.PlottingRaster import MapFigure
    # from LSDMapFigure.PlottingRaster import BaseRaster
    # import LSDPlottingTools.LSDMap_VectorTools as LSDMap_VT
    # import LSDPlottingTools.LSDMap_PointTools as LSDMap_PT

    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # set figure sizes based on format
    if size_format == "geomorphology":
        fig_width_inches = 6.25
    elif size_format == "big":
        fig_width_inches = 16
    else:
        fig_width_inches = 4.92126



    # get the basin IDs to make a discrete colourmap for each ID
    BasinInfoDF = PlotHelp.ReadBasinInfoCSV(DataDirectory, fname_prefix)

    basin_keys = list(BasinInfoDF['basin_key'])
    basin_keys = [int(x) for x in basin_keys]

    basin_junctions = list(BasinInfoDF['outlet_junction'])
    basin_junctions = [float(x) for x in basin_junctions]
    
    
    # get the junctions to mask
    key_to_index_dict = dict(zip(basin_keys,basin_junctions))
    junctions_to_mask = []
    for basin in Basins_to_mask:
        junctions_to_mask.append( key_to_index_dict[basin])
    print("The junctions to mask are")
    print(junctions_to_mask)

    print ('Basin keys are: ')
    print basin_keys
    
    print("Let me mask those for you")
    new_keys = []
    for key in basin_keys:
        if key in Basins_to_mask:
            new_keys.append(np.nan)
        else:
            new_keys.append(key)
    print("The new keys are: ")
    print(new_keys)
    basin_keys = new_keys

    # get a discrete colormap
    cmap = plt.cm.jet

    # going to make the basin plots - need to have bil extensions.
    print("I'm going to make the basin plots. Your topographic data must be in ENVI bil format or I'll break!!")

    # get the rasters
    raster_ext = '.bil'
    #BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext
    BasinsName = fname_prefix+'_AllBasins'+raster_ext
    print (BasinsName)

    # create the map figure
    # We set colourbar location to none since we are labelling the figures
    MF = MapFigure(HillshadeName, DataDirectory,coord_type="UTM_km", colourbar_location='none')

    # add the basins drape
    MF.add_drape_image(BasinsName, DataDirectory, colourmap = cmap, alpha = 0.8, colorbarlabel='Basin ID', discrete_cmap=True, n_colours=len(basin_keys), show_colourbar = True, modify_raster_values=True, old_values=basin_junctions, new_values=basin_keys, cbar_type = int)

    # add the basin outlines
    Basins = LSDP.GetBasinOutlines(DataDirectory, BasinsName)
    
    # get rid of the basins that are being masked
    for junction in junctions_to_mask:
        del Basins[junction]
        
        
        
    
    # note that at this stage the Basins are keyed with the junction index
    MF.plot_polygon_outlines(Basins, linewidth=0.8)

    # add the basin labelling
    label_dict = dict(zip(basin_junctions,basin_keys))
    # this dict has the basin junction as the key and the basin_key as the value
    

    Points = LSDP.GetPointWithinBasins(DataDirectory, BasinsName)
    
    # get rid of points as well
    for junction in junctions_to_mask:
        del Points[junction]
        del label_dict[junction]
    
    MF.add_text_annotation_from_shapely_points(Points, text_colour='k', label_dict=label_dict)

    # Save the figure
    ImageName = DataDirectory+fname_prefix+'_labelled_basins.'+FigFormat
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 250)

def ExampleOne_PartFive_MaskBasinsAndRename(DataDirectory, fname_prefix):
    """
    This function makes a shaded relief plot of the DEM with the basins coloured
    by the basin ID.

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files

    Returns:
        Shaded relief plot with the basins coloured by basin ID

    Author: SMM
    """

    FigFormat = "png"
    size_format = "geomorphology"


    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # set figure sizes based on format
    if size_format == "geomorphology":
        fig_width_inches = 6.25
    elif size_format == "big":
        fig_width_inches = 16
    else:
        fig_width_inches = 4.92126

    # going to make the basin plots - need to have bil extensions.
    print("I'm going to make the basin plots. Your topographic data must be in ENVI bil format or I'll break!!")

    # get the rasters
    raster_ext = '.bil'
    #BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext
    BasinsName = fname_prefix+'_AllBasins'+raster_ext
    print (BasinsName)

    # create the map figure
    # We set colourbar location to none since we are labelling the figures
    MF = MapFigure(HillshadeName, DataDirectory,coord_type="UTM_km", colourbar_location='bottom',basemap_colourmap = "gray")

    # add the basins drape
    Remove_Basins = [4,8]
    Rename_Basins = { 12: 'chumbox', 14: 'zeppo'}
    Value_dict= { 1: 0.2, 2:0.3, 3:0.4, 5:0.9,6:0.7, 7:0.3, 9:0.5, 10:0.5}
    MF.add_basin_plot(BasinsName,fname_prefix,DataDirectory, mask_list = Remove_Basins, 
                      rename_dict = Rename_Basins, value_dict = Value_dict,
                      use_keys_not_junctions = True, show_colourbar = True, 
                      discrete_cmap=True, n_colours=8, colorbarlabel = "$m/n$",
                      colourmap = plt.cm.jet, adjust_text = False)
    
    # Save the figure
    ImageName = DataDirectory+fname_prefix+'_test_Coloured_basins.'+FigFormat
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 250)