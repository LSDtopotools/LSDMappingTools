## LSDMap_ChiPlotting.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## These functions are tools to deal with chi maps
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## SMM
## 14/12/2016
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

import osgeo.gdal as gdal
import numpy as np
import numpy.ma as ma
from osgeo import osr
from os.path import exists
from osgeo.gdalconst import GA_ReadOnly
from numpy import uint8
import matplotlib.pyplot as plt
#from cycler import cycler
from matplotlib import rcParams
import LSDMap_GDALIO as LSDMap_IO
import LSDMap_BasicManipulation as LSDMap_BM
import LSDOSystemTools as LSDOst
import LSDMap_BasicPlotting as LSDMap_BP
import LSDMap_PointData as LSDMap_PD

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## This function plots the chi slope on a shaded relief map
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def BasicChiPlotGridPlot(FileName, DrapeName, chi_csv_fname, thiscmap='gray',drape_cmap='gray',
                            colorbarlabel='Elevation in meters',clim_val = (0,0),
                            drape_alpha = 0.6,FigFileName = 'Image.pdf',FigFormat = 'show',
                            elevation_threshold = 0):
    
    import matplotlib.lines as mpllines
    from mpl_toolkits.axes_grid1 import AxesGrid
    from matplotlib import colors

    label_size = 10
    #title_size = 30
    axis_size = 12

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size 
    #plt.rc('text', usetex=True)

    # get the data
    raster = LSDMap_IO.ReadRasterArrayBlocks(FileName)
    raster_drape = LSDMap_IO.ReadRasterArrayBlocks(DrapeName)

    # SMM attempting mask nodata vals
    NoDataValue = -9999    
    nodata_mask = raster == NoDataValue
    raster[nodata_mask] = np.nan

    nodata_mask2 = raster_drape == NoDataValue
    raster_drape[nodata_mask2] = np.nan 

    #print("Min, max raster is: "+ str(np.nanmax(raster))+","+str(np.min(raster)))
    #print("Min, max raster_drape is: "+ str(np.nanmax(raster_drape))+","+str(np.min(raster_drape)))


    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(FileName)
    
    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # make a figure, sized for a ppt slide
    fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))

    gs = plt.GridSpec(100,100,bottom=0.25,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])
    
    # This is the axis for the colorbar
    ax2 = fig.add_subplot(gs[10:15,15:70])

    #grid = AxesGrid(fig, 111, 
    #                nrows_ncols=(1, 1),
    #                axes_pad=(0.45, 0.15),
    #                label_mode="1",
    #                share_all=True,
    #                cbar_location="right",
    #                cbar_mode="each",
    #                cbar_size="7%",
    #                cbar_pad="2%",
    #                )


    # now get the tick marks    
    n_target_tics = 5
    xlocs,ylocs,new_x_labels,new_y_labels = LSDMap_BP.GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics)  


    print "xmax: " + str(x_max)
    print "xmin: " + str(x_min)
    print "ymax: " + str(y_max)
    print "ymin: " + str(y_min)
    
    #Z1 = np.array(([0, 1]*4 + [1, 0]*4)*4)
    #Z1.shape = (8, 8)  # chessboard
    #im2 = ax.imshow(Z1, cmap=plt.cm.gray, interpolation='nearest',
    #             extent=extent_raster)  
 
    #plt.hold(True)

    im1 = ax.imshow(raster[::-1], thiscmap, extent = extent_raster, interpolation="nearest")
    

    # set the colour limits
    print "Setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1])
    if (clim_val == (0,0)):
        print "Im setting colour limits based on minimum and maximum values"
        im1.set_clim(0, np.nanmax(raster))
    else:
        print "Now setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1])
        im1.set_clim(clim_val[0],clim_val[1])
   
    plt.hold(True)
   
    # Now for the drape: it is in grayscale
    #print "drape_cmap is: "+drape_cmap
    im3 = ax.imshow(raster_drape[::-1], drape_cmap, extent = extent_raster, alpha = drape_alpha, interpolation="nearest")

    # Set the colour limits of the drape
    im3.set_clim(0,np.nanmax(raster_drape))

    
    ax.spines['top'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1) 
     
    #ax.spines['bottom'].set_capstyle('projecting')

    #for spine in ax.spines.values():
    #    spine.set_capstyle('projecting')


    
    ax.set_xticklabels(new_x_labels,rotation=60)
    ax.set_yticklabels(new_y_labels)  
    
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")  

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap        
    ax.tick_params(axis='both', width=1, pad = 2)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(2)    

    # Now we get the chi points
    EPSG_string = LSDMap_IO.GetUTMEPSG(FileName)
    print "EPSG string is: " + EPSG_string
    
    thisPointData = LSDMap_PD.LSDMap_PointData(chi_csv_fname) 
    thisPointData.ThinData('elevation',elevation_threshold)
    
    # convert to easting and northing
    [easting,northing] = thisPointData.GetUTMEastingNorthing(EPSG_string)
    
    
    # The image is inverted so we have to invert the northing coordinate
    Ncoord = np.asarray(northing)
    Ncoord = np.subtract(extent_raster[3],Ncoord)
    Ncoord = np.add(Ncoord,extent_raster[2])
    
    M_chi = thisPointData.QueryData('m_chi')
    #print M_chi
    M_chi = [float(x) for x in M_chi]
    
    
    # make a color map of fixed colors
    this_cmap = colors.ListedColormap(['#2c7bb6','#abd9e9','#ffffbf','#fdae61','#d7191c'])
    bounds=[0,50,100,175,250,1205]
    norm = colors.BoundaryNorm(bounds, this_cmap.N)
    
    sc = ax.scatter(easting,Ncoord,s=0.5, c=M_chi,cmap=this_cmap,norm=norm,edgecolors='none')

    # This affects all axes because we set share_all = True.
    ax.set_xlim(x_min,x_max)    
    ax.set_ylim(y_max,y_min)     

    ax.set_xticks(xlocs)
    ax.set_yticks(ylocs)   
    
    cbar = plt.colorbar(sc,cmap=this_cmap,norm=norm,spacing='uniform', ticks=bounds, boundaries=bounds,orientation='horizontal',cax=ax2)
    cbar.set_label(colorbarlabel, fontsize=10)
    ax2.set_xlabel(colorbarlabel, fontname='Arial',labelpad=-35)    

    print "The figure format is: " + FigFormat
    if FigFormat == 'show':    
        plt.show()
    elif FigFormat == 'return':
        return fig 
    else:
        plt.savefig(FigFileName,format=FigFormat,dpi=500)
        fig.clf()
        
        
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## This function plots channels, color coded
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def BasicChannelPlotGridPlotCategories(FileName, DrapeName, chi_csv_fname, thiscmap='gray',drape_cmap='gray',
                            colorbarlabel='Elevation in meters',clim_val = (0,0),
                            drape_alpha = 0.6,FigFileName = 'Image.pdf',FigFormat = 'show',
                            elevation_threshold = 0, data_name = 'source_key'):
    
    import matplotlib.lines as mpllines
    from mpl_toolkits.axes_grid1 import AxesGrid
    from matplotlib import colors

    label_size = 10
    #title_size = 30
    axis_size = 12

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size 
    #plt.rc('text', usetex=True)

    # get the data
    raster = LSDMap_IO.ReadRasterArrayBlocks(FileName)
    raster_drape = LSDMap_IO.ReadRasterArrayBlocks(DrapeName)
 
    # SMM attempting mask nodata vals
    NoDataValue = -9999    
    nodata_mask = raster == NoDataValue
    raster[nodata_mask] = np.nan

    nodata_mask2 = raster_drape == NoDataValue
    raster_drape[nodata_mask2] = np.nan 

   
    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(FileName)
    
    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # make a figure, sized for a ppt slide
    fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))

    gs = plt.GridSpec(100,100,bottom=0.25,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])
    
    # This is the axis for the colorbar
    #ax2 = fig.add_subplot(gs[10:15,15:70])

    #grid = AxesGrid(fig, 111, 
    #                nrows_ncols=(1, 1),
    #                axes_pad=(0.45, 0.15),
    #                label_mode="1",
    #                share_all=True,
    #                cbar_location="right",
    #                cbar_mode="each",
    #                cbar_size="7%",
    #                cbar_pad="2%",
    #                )


    # now get the tick marks    
    n_target_tics = 5
    xlocs,ylocs,new_x_labels,new_y_labels = LSDMap_BP.GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics)  


    print "xmax: " + str(x_max)
    print "xmin: " + str(x_min)
    print "ymax: " + str(y_max)
    print "ymin: " + str(y_min)


    #Z1 = np.array(([0, 1]*4 + [1, 0]*4)*4)
    #Z1.shape = (8, 8)  # chessboard
    #im2 = ax.imshow(Z1, cmap=plt.cm.gray, interpolation='nearest',
    #             extent=extent_raster)  
 
    #plt.hold(True)

    im1 = ax.imshow(raster[::-1], thiscmap, extent = extent_raster, interpolation="nearest")
    

    # set the colour limits
    print "Setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1])
    if (clim_val == (0,0)):
        print "Im setting colour limits based on minimum and maximum values"
        im1.set_clim(0, np.nanmax(raster))
    else:
        print "Now setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1])
        im1.set_clim(clim_val[0],clim_val[1])
   
    plt.hold(True)
   
    # Now for the drape: it is in grayscale
    #print "drape_cmap is: "+drape_cmap
    im3 = ax.imshow(raster_drape[::-1], drape_cmap, extent = extent_raster, alpha = drape_alpha, interpolation="nearest")

    # Set the colour limits of the drape
    im3.set_clim(0,np.nanmax(raster_drape))
    
    
    ax.spines['top'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1) 
     
    #ax.spines['bottom'].set_capstyle('projecting')

    #for spine in ax.spines.values():
    #    spine.set_capstyle('projecting')


    
    ax.set_xticklabels(new_x_labels,rotation=60)
    ax.set_yticklabels(new_y_labels)  
    
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")  

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap        
    ax.tick_params(axis='both', width=1, pad = 2)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(2)    

    # Now we get the chi points
    EPSG_string = LSDMap_IO.GetUTMEPSG(FileName)
    print "EPSG string is: " + EPSG_string
    
    thisPointData = LSDMap_PD.LSDMap_PointData(chi_csv_fname) 
    thisPointData.ThinData('elevation',elevation_threshold)
    
    # convert to easting and northing
    [easting,northing] = thisPointData.GetUTMEastingNorthing(EPSG_string)
 
    # The image is inverted so we have to invert the northing coordinate
    Ncoord = np.asarray(northing)
    Ncoord = np.subtract(extent_raster[3],Ncoord)
    Ncoord = np.add(Ncoord,extent_raster[2])
    
    these_data = thisPointData.QueryData(data_name)
    #print M_chi
    these_data = [int(x) for x in these_data]

    # make a color map of fixed colors
    NUM_COLORS = 15

    this_cmap = plt.cm.Set1
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    scalarMap = plt.cm.ScalarMappable(norm=cNorm, cmap=this_cmap)
    channel_data = [x % NUM_COLORS for x in these_data]

    sc = ax.scatter(easting,Ncoord,s=0.5, c=channel_data,norm=cNorm,cmap=this_cmap,edgecolors='none')

    # This affects all axes because we set share_all = True.
    ax.set_xlim(x_min,x_max)    
    ax.set_ylim(y_max,y_min)     

    ax.set_xticks(xlocs)
    ax.set_yticks(ylocs)   
    ax.set_title('Channels colored by source number')
    #ax.text(1.1, 0.01, 'Channels colored by source number',
    #    verticalalignment='top', horizontalalignment='right',
    #    transform=ax.transAxes,
    #    color='green', fontsize=15)
    
    #cbar = plt.colorbar(sc,cmap=this_cmap,norm=cNorm,spacing='uniform', orientation='horizontal',cax=ax2)
    #cbar.set_label(colorbarlabel, fontsize=10)
    #ax2.set_xlabel(colorbarlabel, fontname='Arial',labelpad=-35)    

    print "The figure format is: " + FigFormat
    if FigFormat == 'show':    
        plt.show()
    elif FigFormat == 'return':
        return fig 
    else:
        plt.savefig(FigFileName,format=FigFormat,dpi=500)
        fig.clf() 
        
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## This function plots channels, color coded
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        
def ChiProfiles(chi_csv_fname, FigFileName = 'Image.pdf',FigFormat = 'show',
                elevation_threshold = 0):

    import matplotlib.lines as mpllines
    from mpl_toolkits.axes_grid1 import AxesGrid
    from matplotlib import colors

    label_size = 10
    #title_size = 30
    axis_size = 12

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size     
   

    # make a figure, sized for a ppt slide
    fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))

    gs = plt.GridSpec(100,100,bottom=0.25,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])

    # Now we get the chi points
    #EPSG_string = LSDMap_IO.GetUTMEPSG(FileName)
    #print "EPSG string is: " + EPSG_string
    
    thisPointData = LSDMap_PD.LSDMap_PointData(chi_csv_fname) 
    thisPointData.ThinData('elevation',elevation_threshold)
    
    # Get the chi, m_chi, basin number, and source ID code
    chi = thisPointData.QueryData('chi')
    chi = [float(x) for x in chi]
    elevation = thisPointData.QueryData('elevation')
    elevation = [float(x) for x in elevation]
    fdist = thisPointData.QueryData('flow distance')
    fdist = [float(x) for x in fdist]     
    m_chi = thisPointData.QueryData('m_chi')
    m_chi = [float(x) for x in m_chi]    
    basin = thisPointData.QueryData('basin_key')
    basin = [int(x) for x in basin] 
    source = thisPointData.QueryData('source_key')
    source = [int(x) for x in source]
    
    #print source
    
 
    # need to convert everything into arrays so we can mask different basins
    Chi = np.asarray(chi)
    Elevation = np.asarray(elevation)
    Fdist = np.asarray(fdist)
    M_chi = np.asarray(m_chi)
    Basin = np.asarray(basin)
    Source = np.asarray(source)
    
    max_basin = np.amax(Basin)
    max_chi = np.amax(Chi)
    max_Elevation = np.amax(Elevation)
    max_M_chi = np.amax(M_chi)
    min_Elevation = np.amin(Elevation)
    
    z_axis_min = int(min_Elevation/10)*10 
    z_axis_max = int(max_Elevation/10)*10+10
    chi_axis_max = int(max_chi/5)*5+5
    
    # make a color map of fixed colors
    NUM_COLORS = 15

    this_cmap = plt.cm.Set1
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    #scalarMap = plt.cm.ScalarMappable(norm=cNorm, cmap=this_cmap)      
    Source_colors = [x % NUM_COLORS for x in Source]
    
    #print "The source colours!!"
    #print Source_colours
    
    # create a mask
    basin_number = 10
    
    dot_pos = FigFileName.rindex('.')
    newFilename = FigFileName[:dot_pos]+'_'+str(basin_number)+FigFileName[dot_pos:]
    print "newFilename: "+newFilename
    
    
    m = np.ma.masked_where(Basin!=basin_number, Basin)
    maskChi = np.ma.masked_where(np.ma.getmask(m), Chi)
    maskElevation = np.ma.masked_where(np.ma.getmask(m), Elevation)
    maskSource = np.ma.masked_where(np.ma.getmask(m), Source_colors)
        
    source_colors = [x % NUM_COLORS for x in maskSource]

    sc = ax.scatter(maskChi,maskElevation,s=2.0, c=Source_colors,norm=cNorm,cmap=this_cmap,edgecolors='none')

    ax.spines['top'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1) 

    ax.set_xlabel("$\chi$")
    ax.set_ylabel("Elevation (m)") 
    
    # This affects all axes because we set share_all = True.
    ax.set_ylim(z_axis_min,z_axis_max)    
    ax.set_xlim(0,chi_axis_max)      

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap        
    ax.tick_params(axis='both', width=1, pad = 2)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(2)       

    print "The figure format is: " + FigFormat
    if FigFormat == 'show':    
        plt.show()
    elif FigFormat == 'return':
        return fig 
    else:
        plt.savefig(newFilename,format=FigFormat,dpi=500)
        fig.clf()     
        
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## This function plots channels, color coded
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        
def StackedChiProfiles(chi_csv_fname, FigFileName = 'Image.pdf',
                       FigFormat = 'show',elevation_threshold = 0, 
                       first_basin = 0, last_basin = 0):

    import matplotlib.lines as mpllines
    from mpl_toolkits.axes_grid1 import AxesGrid
    from matplotlib import colors

    label_size = 10
    #title_size = 30
    axis_size = 12

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size     
   

    # make a figure, sized for a ppt slide
    fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))

    gs = plt.GridSpec(100,100,bottom=0.25,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])

    thisPointData = LSDMap_PD.LSDMap_PointData(chi_csv_fname) 
    thisPointData.ThinData('elevation',elevation_threshold)
    
    # Get the chi, m_chi, basin number, and source ID code
    chi = thisPointData.QueryData('chi')
    chi = [float(x) for x in chi]
    elevation = thisPointData.QueryData('elevation')
    elevation = [float(x) for x in elevation]
    fdist = thisPointData.QueryData('flow distance')
    fdist = [float(x) for x in fdist]     
    m_chi = thisPointData.QueryData('m_chi')
    m_chi = [float(x) for x in m_chi]    
    basin = thisPointData.QueryData('basin_key')
    basin = [int(x) for x in basin] 
    source = thisPointData.QueryData('source_key')
    source = [int(x) for x in source]

    # need to convert everything into arrays so we can mask different basins
    Chi = np.asarray(chi)
    Elevation = np.asarray(elevation)
    Fdist = np.asarray(fdist)
    M_chi = np.asarray(m_chi)
    Basin = np.asarray(basin)
    Source = np.asarray(source)
    
    max_basin = np.amax(Basin)
    max_chi = np.amax(Chi)
    max_Elevation = np.amax(Elevation)
    max_M_chi = np.amax(M_chi)
    min_Elevation = np.amin(Elevation)
    
    z_axis_min = int(min_Elevation/10)*10 
    z_axis_max = int(max_Elevation/10)*10+10
    chi_axis_max = int(max_chi/5)*5+5
    
    # make a color map of fixed colors
    NUM_COLORS = 15

    this_cmap = plt.cm.Set1
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    #scalarMap = plt.cm.ScalarMappable(norm=cNorm, cmap=this_cmap)      
    Source_colors = [x % NUM_COLORS for x in Source]
    
    # now loop through a number of basins
    if last_basin >= max_basin:
        last_basin = max_basin-1
    
    if first_basin > last_basin:
        first_basin = last_basin
        print("Your first basin was larger than last basin. I won't plot anything")
     
    plt.hold(True) 
    chi_offset = 5
    this_chi_offset = 0
    
    n_stacks = last_basin-first_basin+1
    added_chi = chi_offset*n_stacks
    print("The number of stacks is: "+str(n_stacks)+" the old max: "+str(chi_axis_max))    
    chi_axis_max = chi_axis_max+added_chi
    print("The nex max is: "+str(chi_axis_max))
    
    
    dot_pos = FigFileName.rindex('.')
    newFilename = FigFileName[:dot_pos]+'_Stack'+str(first_basin)+FigFileName[dot_pos:]
    basins_list = list(range(first_basin,last_basin+1))
    
    for basin_number in basins_list:
        
        print ("This basin is: " +str(basin_number))

        m = np.ma.masked_where(Basin!=basin_number, Basin)
        maskChi = np.ma.masked_where(np.ma.getmask(m), Chi)
        maskElevation = np.ma.masked_where(np.ma.getmask(m), Elevation)
        maskSource = np.ma.masked_where(np.ma.getmask(m), Source_colors)
        
        print("adding an offset of: "+str(this_chi_offset))
        
        maskChi = np.add(maskChi,this_chi_offset)
        this_chi_offset = this_chi_offset+chi_offset
        
        if basin_number == basins_list[-1]:
            print("last basin, geting maximum value,basin is: "+str(basin_number))
            this_max = np.amax(maskChi)
            this_max = int(this_max/5)*5+5
            print("The rounded maximum is: "+str(this_max))
            chi_axis_max = this_max
        
        source_colors = [x % NUM_COLORS for x in maskSource]

        sc = ax.scatter(maskChi,maskElevation,s=2.0, c=Source_colors,norm=cNorm,cmap=this_cmap,edgecolors='none')
        
    ax.spines['top'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1) 

    ax.set_xlabel("$\chi$ (m)")
    ax.set_ylabel("Elevation (m)") 
    
    # This affects all axes because we set share_all = True.
    ax.set_ylim(z_axis_min,z_axis_max)    
    ax.set_xlim(0,chi_axis_max)      

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap        
    ax.tick_params(axis='both', width=1, pad = 2)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(2)       

    print "The figure format is: " + FigFormat
    if FigFormat == 'show':    
        plt.show()
    elif FigFormat == 'return':
        return fig 
    else:
        plt.savefig(newFilename,format=FigFormat,dpi=500)
        fig.clf()           
        
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## This function plots channels, color coded in chi space with a gradient
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        
def StackedProfilesGradient(chi_csv_fname, FigFileName = 'Image.pdf',
                       FigFormat = 'show',elevation_threshold = 0, 
                       first_basin = 0, last_basin = 0,
                       this_cmap = plt.cm.cubehelix,data_name = 'chi', X_offset = 5,
                       plotting_data_format = 'log'):

    import math
    import matplotlib.lines as mpllines
    from mpl_toolkits.axes_grid1 import AxesGrid
    from matplotlib import colors

    label_size = 10
    #title_size = 30
    axis_size = 12

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size     
   

    # make a figure, sized for a ppt slide
    fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))

    gs = plt.GridSpec(100,100,bottom=0.25,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])

    thisPointData = LSDMap_PD.LSDMap_PointData(chi_csv_fname) 
    thisPointData.ThinData('elevation',elevation_threshold)
    
    
    # Get the chi, m_chi, basin number, and source ID code
    if data_name  == 'chi':
        x_data = thisPointData.QueryData('chi')
        x_data = [float(x) for x in x_data]
    elif data_name == 'flow_distance':
        x_data = thisPointData.QueryData('flow distance')
        x_data = [float(x) for x in fdist]   
    else:
        print("I did not understand the data name. Choices are chi and elevation. Defaulting to chi.")
        x_data = thisPointData.QueryData('chi')
        x_data = [float(x) for x in x_data] 

    elevation = thisPointData.QueryData('elevation')
    elevation = [float(x) for x in elevation]        
    m_chi = thisPointData.QueryData('m_chi')
    m_chi = [float(x) for x in m_chi]    
    basin = thisPointData.QueryData('basin_key')
    basin = [int(x) for x in basin] 
    source = thisPointData.QueryData('source_key')
    source = [int(x) for x in source]

    colorbarlabel = "$k_{sn}$"
    if (plotting_data_format == 'log'):
        log_m_chi = []
        for value in m_chi:
            if value < 0.1:
                log_m_chi.append(0)
            else:
                log_m_chi.append(math.log10(value))
        m_chi = log_m_chi
        colorbarlabel = "log$_{10}k_{sn}$"
    
    # need to convert everything into arrays so we can mask different basins
    Xdata = np.asarray(x_data)
    Elevation = np.asarray(elevation)
    M_chi = np.asarray(m_chi)
    Basin = np.asarray(basin)
    Source = np.asarray(source)  
    
    max_basin = np.amax(Basin)
    max_X = np.amax(Xdata)
    max_Elevation = np.amax(Elevation)
    max_M_chi = np.amax(M_chi)
    min_Elevation = np.amin(Elevation)
    
    print("Max M_chi is: "+str(max_M_chi))
    
    z_axis_min = int(min_Elevation/10)*10 
    z_axis_max = int(max_Elevation/10)*10+10
    X_axis_max = int(max_X/5)*5+5
    M_chi_axis_max = max_M_chi
    
    # now loop through a number of basins
    if last_basin >= max_basin:
        last_basin = max_basin-1
    
    if first_basin > last_basin:
        first_basin = last_basin
        print("Your first basin was larger than last basin. I won't plot anything")
     
    plt.hold(True)
    
    # X offest is now given by user so they can tune to the channel profiles
    #X_offset = 5
    this_X_offset = 0
    
    n_stacks = last_basin-first_basin+1
    added_X = X_offset*n_stacks
    print("The number of stacks is: "+str(n_stacks)+" the old max: "+str(X_axis_max))    
    X_axis_max = X_axis_max+added_X
    print("The nex max is: "+str(X_axis_max))
    
    
    dot_pos = FigFileName.rindex('.')
    newFilename = FigFileName[:dot_pos]+'_GradientStack'+str(first_basin)+FigFileName[dot_pos:]
    basins_list = list(range(first_basin,last_basin+1))
    
    for basin_number in basins_list:
        
        print ("This basin is: " +str(basin_number))

        m = np.ma.masked_where(Basin!=basin_number, Basin)
        maskX = np.ma.masked_where(np.ma.getmask(m), Xdata)
        maskElevation = np.ma.masked_where(np.ma.getmask(m), Elevation)
        maskMChi = np.ma.masked_where(np.ma.getmask(m), M_chi)
        #maskSource = np.ma.masked_where(np.ma.getmask(m), Source_colors)
        
        print("adding an offset of: "+str(this_X_offset))
        
        maskX = np.add(maskX,this_X_offset)
        this_X_offset = this_X_offset+X_offset
        
        if basin_number == basins_list[-1]:
            print("last basin, geting maximum value,basin is: "+str(basin_number))
            this_max = np.amax(maskX)
            this_max = int(this_max/5)*5+5
            print("The rounded maximum is: "+str(this_max))
            X_axis_max = this_max
        
        sc = ax.scatter(maskX,maskElevation,s=2.0, c=maskMChi,cmap=this_cmap,edgecolors='none')
 
    # set the colour limits
    sc.set_clim(0, M_chi_axis_max)
    bounds = (0, M_chi_axis_max)

    # This is the axis for the colorbar
    
    ax2 = fig.add_subplot(gs[10:15,15:70])
    cbar = plt.colorbar(sc,cmap=this_cmap,spacing='uniform', orientation='horizontal',cax=ax2)   
    cbar.set_label(colorbarlabel, fontsize=10)
    ax2.set_xlabel(colorbarlabel, fontname='Arial',labelpad=-35)        

    ax.spines['top'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1) 

    ax.set_xlabel("$\chi$ (m)")
    ax.set_ylabel("Elevation (m)") 
    
    # This affects all axes because we set share_all = True.
    ax.set_ylim(z_axis_min,z_axis_max)    
    ax.set_xlim(0,X_axis_max)      

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap        
    ax.tick_params(axis='both', width=1, pad = 2)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(2)       

    print "The figure format is: " + FigFormat
    if FigFormat == 'show':    
        plt.show()
    elif FigFormat == 'return':
        return fig 
    else:
        plt.savefig(newFilename,format=FigFormat,dpi=500)
        fig.clf()                   