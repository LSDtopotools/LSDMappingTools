"""
    This contains wrapper functions that simplify plotting raster
    and vector data for publication-ready figures.

    The documentation of the examples can be found here:
    https://lsdtopotools.github.io/LSDTopoTools_ChiMudd2014/

    Simon Mudd and Fiona Clubb, January 2018

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
import LSDPlottingTools.LSDMap_PointTools as LSDMap_PD
from LSDMapFigure.PlottingRaster import MapFigure
import LSDMapFigure.PlottingHelpers as PlotHelp
import LSDPlottingTools.LSDMap_ChiPlotting as LSDCP
#import LSDPlottingTools.adjust_text


def PrintChiChannels(DataDirectory,fname_prefix, ChannelFileName, add_basin_labels = True, cmap = "jet", cbar_loc = "right", size_format = "ESURF", fig_format = "png", dpi = 250,plotting_column = "source_key",discrete_colours = False, NColours = 10, out_fname_prefix = ""):
    """
    This function prints a channel map over a hillshade.

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        ChannelFileName (str): The name of the channel file (a csv) without path but with extension
        add_basin_labels (bool): If true, label the basins with text. Otherwise use a colourbar.
        cmap (str or colourmap): The colourmap to use for the plot
        cbar_lox (str): where you want the colourbar. Options are none, left, right, top and botton. The colourbar will be of the elevation.
                        If you want only a hillshade set to none and the cmap to "gray"
        size_format (str): Either geomorphology or big. Anything else gets you a 4.9 inch wide figure (standard ESURF size)
        fig_format (str): An image format. png, pdf, eps, svg all valid
        dpi (int): The dots per inch of the figure
        plotting_column (str): the name of the column to plot
        discrete_colours (bool): if true use a discrete colourmap
        NColours (int): the number of colours to cycle through when making the colourmap
        out_fname_prefix (str): The prefix of the image file. If blank uses the fname_prefix


    Returns:
        Shaded relief plot with the channels coloured by a plotting column designated by the plotting_column keyword. Uses a colourbar to show each basin

    Author: SMM
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
    BackgroundRasterName = fname_prefix+"_hs.bil"
    DrapeRasterName = fname_prefix+".bil"
    chi_csv_fname = DataDirectory+ChannelFileName

    thisPointData = LSDMap_PD.LSDMap_PointData(chi_csv_fname)


    # clear the plot
    plt.clf()

    # set up the base image and the map
    MF = MapFigure(BackgroundRasterName, DataDirectory, coord_type="UTM_km",colourbar_location = "None")
    MF.add_drape_image(DrapeRasterName,DataDirectory,colourmap = "gray", alpha = 0.6)
    MF.add_point_data(thisPointData,column_for_plotting = plotting_column,this_colourmap = cmap,
                       scale_points = True,column_for_scaling = "drainage_area",
                       scaled_data_in_log = True,
                       max_point_size = 5, min_point_size = 1,discrete_colours = discrete_colours, NColours = NColours)

    # Save the image
    if len(out_fname_prefix) == 0:
        ImageName = DataDirectory+fname_prefix+"_chi_channels."+fig_format
    else:
        ImageName = DataDirectory+out_fname_prefix+"_chi_channels."+fig_format

    MF.save_fig(fig_width_inches = fig_size_inches, FigFileName = ImageName, axis_style = ax_style, FigFormat=fig_format, Fig_dpi = dpi)



def PrintChiChannelsAndBasins(DataDirectory,fname_prefix, ChannelFileName, add_basin_labels = True, cmap = "jet", cbar_loc = "right", size_format = "ESURF", fig_format = "png", dpi = 250,plotting_column = "source_key",discrete_colours = False, NColours = 10, colour_log = True, colorbarlabel = "Colourbar", Basin_remove_list = [], Basin_rename_dict = {} , value_dict = {}, out_fname_prefix = "", show_basins = True, min_channel_point_size = 0.5, max_channel_point_size = 2):
    """
    This function prints a channel map over a hillshade.

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        add_basin_labels (bool): If true, label the basins with text. Otherwise use a colourbar.
        cmap (str or colourmap): The colourmap to use for the plot
        cbar_loc (str): where you want the colourbar. Options are none, left, right, top and botton. The colourbar will be of the elevation.
                        If you want only a hillshade set to none and the cmap to "gray"
        size_format (str): Either geomorphology or big. Anything else gets you a 4.9 inch wide figure (standard ESURF size)
        fig_format (str): An image format. png, pdf, eps, svg all valid
        dpi (int): The dots per inch of the figure
        plotting_column (str): the name of the column to plot
        discrete_colours (bool): if true use a discrete colourmap
        NColours (int): the number of colours to cycle through when making the colourmap
        colour_log (bool): If true the colours are in log scale
        Basin_remove_list (list): A lists containing either key or junction indices of basins you want to remove from plotting
        Basin_rename_dict (dict): A dict where the key is either basin key or junction index, and the value is a new name for the basin denoted by the key
        out_fname_prefix (str): The prefix of the image file. If blank uses the fname_prefix
        show_basins (bool): If true, plot the basins
        min_channel_point_size (float): The minimum size of a channel point in points
        max_channel_point_size (float): The maximum size of a channel point in points

    Returns:
        Shaded relief plot with the basins coloured by basin ID. Includes channels. These can be plotted by various metrics denoted but the plotting_column parameter.

    Author: SMM
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

    chi_csv_fname = DataDirectory+ChannelFileName
    chi_csv_fname = DataDirectory+ChannelFileName

    thisPointData = LSDMap_PD.LSDMap_PointData(chi_csv_fname)

    #thisPointData.ThinDataSelection("basin_key",[10])

    thisPointData.selectValue("basin_key",value = Basin_remove_list, operator = "!=")
    #print("The new point data is:")
    #print(thisPointData.GetLongitude())


    # clear the plot
    plt.clf()

    # set up the base image and the map
    print("I am showing the basins without text labels.")
    MF = MapFigure(HillshadeName, DataDirectory,coord_type="UTM_km", colourbar_location="None")

    # This adds the basins
    if show_basins:
        MF.add_basin_plot(BasinsName,fname_prefix,DataDirectory, mask_list = Basin_remove_list, rename_dict = Basin_rename_dict, value_dict = value_dict, label_basins = add_basin_labels, show_colourbar = False,
                      colourmap = "gray")

    if discrete_colours:
        print("I am printing discrete colours.")

    MF.add_point_data(thisPointData,column_for_plotting = plotting_column,
                       scale_points = True,column_for_scaling = "drainage_area", show_colourbar = True, colourbar_location = cbar_loc,
                       colorbarlabel = colorbarlabel, this_colourmap = cmap,
                       scaled_data_in_log = True,
                       max_point_size = max_channel_point_size, min_point_size = min_channel_point_size,zorder=10, colour_log = colour_log, discrete_colours = discrete_colours, NColours = NColours)

    # Save the image
    if len(out_fname_prefix) == 0:
        ImageName = DataDirectory+fname_prefix+"_chi_channels_and_basins."+fig_format
    else:
        ImageName = DataDirectory+out_fname_prefix+"_chi_channels_and_basins."+fig_format

    MF.save_fig(fig_width_inches = fig_size_inches, FigFileName = ImageName, axis_style = ax_style, FigFormat=fig_format, Fig_dpi = dpi)



def PrintChiCoordChannelsAndBasins(DataDirectory,fname_prefix, ChannelFileName, add_basin_labels = True, cmap = "cubehelix", cbar_loc = "right", size_format = "ESURF", fig_format = "png", dpi = 250,plotting_column = "source_key",discrete_colours = False, NColours = 10, colour_log = True, colorbarlabel = "Colourbar", Basin_remove_list = [], Basin_rename_dict = {} , value_dict = {}, plot_chi_raster = False, out_fname_prefix = "", show_basins = True, min_channel_point_size = 0.5, max_channel_point_size = 2):
    """
    This function prints a channel map over a hillshade.

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        add_basin_labels (bool): If true, label the basins with text. Otherwise use a colourbar.
        cmap (str or colourmap): The colourmap to use for the plot
        cbar_loc (str): where you want the colourbar. Options are none, left, right, top and botton. The colourbar will be of the elevation.
                        If you want only a hillshade set to none and the cmap to "gray"
        size_format (str): Either geomorphology or big. Anything else gets you a 4.9 inch wide figure (standard ESURF size)
        fig_format (str): An image format. png, pdf, eps, svg all valid
        dpi (int): The dots per inch of the figure
        plotting_column (str): the name of the column to plot
        discrete_colours (bool): if true use a discrete colourmap
        NColours (int): the number of colours to cycle through when making the colourmap
        colour_log (bool): If true the colours are in log scale
        Basin_remove_list (list): A lists containing either key or junction indices of basins you want to remove from plotting
        Basin_rename_dict (dict): A dict where the key is either basin key or junction index, and the value is a new name for the basin denoted by the key
        out_fname_prefix (str): The prefix of the image file. If blank uses the fname_prefix
        show_basins (bool): If true, plot the basins
        min_channel_point_size (float): The minimum size of a channel point in points
        max_channel_point_size (float): The maximum size of a channel point in points

    Returns:
        Shaded relief plot with the basins coloured by basin ID. Includes channels. These can be plotted by various metrics denoted but the plotting_column parameter.

    Author: SMM
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
    ChiCoordName = fname_prefix+'_Maskedchi'+raster_ext
    print (BasinsName)
    Basins = LSDP.GetBasinOutlines(DataDirectory, BasinsName)

    chi_csv_fname = DataDirectory+ChannelFileName
    chi_csv_fname = DataDirectory+ChannelFileName

    thisPointData = LSDMap_PD.LSDMap_PointData(chi_csv_fname)

    #thisPointData.ThinDataSelection("basin_key",[10])

    thisPointData.selectValue("basin_key",value = Basin_remove_list, operator = "!=")
    #print("The new point data is:")
    #print(thisPointData.GetLongitude())


    # clear the plot
    plt.clf()

    # set up the base image and the map
    print("I am showing the basins without text labels.")
    MF = MapFigure(HillshadeName, DataDirectory,coord_type="UTM_km", colourbar_location="None")

    # This adds the basins


    if plot_chi_raster:
        if show_basins:
            MF.add_basin_plot(BasinsName,fname_prefix,DataDirectory, mask_list = Basin_remove_list, rename_dict = Basin_rename_dict, value_dict = value_dict, label_basins = add_basin_labels, show_colourbar = False, colourmap = "gray", alpha = 1, outlines_only = True)

        MF.add_drape_image(ChiCoordName,DataDirectory,colourmap = "cubehelix",alpha=0.6,zorder = 0.5)

        MF.add_point_data(thisPointData,column_for_plotting = plotting_column,scale_points = True,column_for_scaling = "drainage_area", show_colourbar = True, colourbar_location = cbar_loc,colorbarlabel = colorbarlabel, this_colourmap = cmap,scaled_data_in_log = True,max_point_size = max_channel_point_size, min_point_size = min_channel_point_size,zorder=0.4, colour_log = colour_log, discrete_colours = discrete_colours, NColours = NColours)
    else:
        if show_basins:
            MF.add_basin_plot(BasinsName,fname_prefix,DataDirectory, mask_list = Basin_remove_list, rename_dict = Basin_rename_dict, value_dict = value_dict, label_basins = add_basin_labels, show_colourbar = False, colourmap = "gray", alpha = 0.7, outlines_only = False)

        MF.add_point_data(thisPointData,column_for_plotting = plotting_column,scale_points = True,column_for_scaling = "drainage_area", show_colourbar = True, colourbar_location = cbar_loc,colorbarlabel = colorbarlabel, this_colourmap = cmap,scaled_data_in_log = True,max_point_size = 2, min_point_size = 0.5,zorder=10, colour_log = colour_log, discrete_colours = discrete_colours, NColours = NColours)





    # Save the image
    if len(out_fname_prefix) == 0:
        ImageName = DataDirectory+fname_prefix+"_chicoord_and_basins."+fig_format
    else:
        ImageName = DataDirectory+out_fname_prefix+"_chicoord_and_basins."+fig_format

    MF.save_fig(fig_width_inches = fig_size_inches, FigFileName = ImageName, axis_style = ax_style, FigFormat=fig_format, Fig_dpi = dpi)



def PrintChiStacked(DataDirectory,fname_prefix, ChannelFileName, cmap = "jet", cbar_loc = "bottom", size_format = "ESURF", fig_format = "png", dpi = 250,plotting_column = "source_key",discrete_colours = False, NColours = 10,colorbarlabel = "Colourbar", axis_data_name = "chi", plot_data_name = "m_chi", plotting_data_format = 'log', Basin_select_list = [], Basin_rename_dict = {}, out_fname_prefix = "", first_basin = 0, last_basin = 0, figure_aspect_ratio = 2, X_offset = 5, rotate_labels=False):
    """
    This function prints a channel map over a hillshade.

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        add_basin_labels (bool): If true, label the basins with text. Otherwise use a colourbar.
        cmap (str or colourmap): The colourmap to use for the plot
        cbar_loc (str): where you want the colourbar. Options are none, left, right, top and botton. The colourbar will be of the elevation.
                        If you want only a hillshade set to none and the cmap to "gray"
        size_format (str): Either geomorphology or big. Anything else gets you a 4.9 inch wide figure (standard ESURF size)
        fig_format (str): An image format. png, pdf, eps, svg all valid
        dpi (int): The dots per inch of the figure
        plotting_column (str): the name of the column to plot
        discrete_colours (bool): if true use a discrete colourmap
        NColours (int): the number of colours to cycle through when making the colourmap
        Basin_remove_list (list): A lists containing either key or junction indices of basins you want to remove from plotting
        Basin_rename_dict (dict): A dict where the key is either basin key or junction index, and the value is a new name for the basin denoted by the key
        out_fname_prefix (str): The prefix of the image file. If blank uses the fname_prefix
        axis_data_name (str): the data used as the x axis
        plot_data_name (str): the data name used to colour the plot

    Returns:
        Shaded relief plot with the basins coloured by basin ID. Includes channels. These can be plotted by various metrics denoted but the plotting_column parameter.
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

    # get the basin IDs to make a discrete colourmap for each ID
    BasinInfoDF = PlotHelp.ReadBasinInfoCSV(DataDirectory, fname_prefix)

    basin_keys = list(BasinInfoDF['basin_key'])
    basin_keys = [int(x) for x in basin_keys]

    basin_junctions = list(BasinInfoDF['outlet_junction'])
    basin_junctions = [float(x) for x in basin_junctions]

    print ('Basin keys are: ')
    print basin_keys


    chi_csv_fname = DataDirectory+ChannelFileName

    # Save the image
    if len(out_fname_prefix) == 0:
        ImageName = DataDirectory+fname_prefix+"_stacked_chi."+fig_format
    else:
        ImageName = DataDirectory+out_fname_prefix+"_stacked_chi."+fig_format

    if axis_data_name == "flow_distance" and X_offset <= 10:
        print("WARNING! You have a weird flow distance offset. I think it is the chi offset. Check your offset.")
        x_offset = 50000
    else:
        x_offset = X_offset

    # print("The colourbar is located on the "+cbar_loc)
    # print("Cmap is: "+cmap)

    print("About to go into the stacks. My x_offset is: " +str(x_offset)+ ", and my rename dict is:" )
    print(Basin_rename_dict)
    LSDCP.StackedProfilesGradient(chi_csv_fname, FigFileName = ImageName,
                       FigFormat = fig_format,elevation_threshold = 0,
                       first_basin = first_basin, last_basin = last_basin, basin_order_list = Basin_select_list,
                       basin_rename_dict = Basin_rename_dict,
                       this_cmap = cmap,axis_data_name = axis_data_name, colour_data_name = plot_data_name,
                       discrete_colours = discrete_colours, NColours = NColours,
                       colorbarlabel = colorbarlabel, cbar_loc = cbar_loc, X_offset = x_offset,
                       plotting_data_format = plotting_data_format,
                       label_sources = False, source_thinning_threshold = 0,
                       size_format = size_format, aspect_ratio = figure_aspect_ratio, dpi = dpi, rotate_labels=rotate_labels)
