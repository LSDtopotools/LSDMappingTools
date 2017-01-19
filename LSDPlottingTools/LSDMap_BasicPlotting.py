## LSDMap_BasicPlotting.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## These functions are tools to deal with rasters
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## SMM
## 26/07/2014
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

import osgeo.gdal as gdal
import numpy as np
import numpy.ma as ma
from osgeo import osr
from os.path import exists
from osgeo.gdalconst import GA_ReadOnly
from numpy import uint8
from matplotlib import rcParams
from adjust_text import adjust_text
import LSDMap_GDALIO as LSDMap_IO
import LSDMap_BasicManipulation as LSDMap_BM
import LSDOSystemTools as LSDOst
import LSDMap_PointData as LSDMap_PD 
import matplotlib.pyplot as plt


#==============================================================================
# This formats ticks if you want to convert metres to km
#==============================================================================
def TickConverter(x_min,x_max,n_target_tics):
    dx_fig = x_max-x_min
    dx_spacing = dx_fig/n_target_tics
    #print("spacing: "+str(dx_spacing))

    # This extracts the digits before the full stop
    str_dx = str(dx_spacing)
    str_dx = str_dx.split('.')[0]
    n_digits = str_dx.__len__()
    nd = int(n_digits)
    # We are left with the number of digits in the spacing. This will be used
    # to round tick locations
        
    first_digit = float(str_dx[0])
       
    dx_spacing_rounded = first_digit*pow(10,(nd-1))
    #print("dx spacing is: " + str(dx_spacing_rounded))
 
    str_xmin = str(x_min)
    #print("before split str_xmin: "+ str_xmin)
    str_xmin = str_xmin.split('.')[0]
    #print("after split str_xmin: "+ str_xmin)
    x_min = float(str_xmin)
    #print("x_min: "+ str(x_min))
    
    n_digx = str_xmin.__len__() 
      
    if (n_digx-nd+1) >= 1:
        front_x = str_xmin[:(n_digx-nd+1)]
    else:
        front_x = str_xmin
        
    round_xmin = float(front_x)*pow(10,nd-1)
    #print("round xmin is: " + str(round_xmin))
    if round_xmin <0:
        round_xmin = 0

    # now we need to figure out where the xllocs and ylocs are
    xlocs = np.zeros(2*n_target_tics) 
    xlocs_km = np.zeros(2*n_target_tics) 
    new_x_labels = []

    for i in range(0,2*n_target_tics):
        xlocs[i] = round_xmin+(i)*dx_spacing_rounded       
        xlocs_km[i] = xlocs[i]/1000.0      
        new_x_labels.append( str(xlocs_km[i]).split(".")[0] )

    #print xlocs
    #print new_x_labels
    #print new_y_labels

    new_xlocs = []
    new_xlocs_km = []
    x_labels = []

    # Now loop through these to get rid of those not in range
    for index,xloc in enumerate(xlocs):
        #print xloc
        if (xloc <= x_max and xloc >= x_min):
            new_xlocs.append(xloc)
            new_xlocs_km.append(xlocs_km[index])
            x_labels.append(new_x_labels[index])
 
            
    #print "======================================="
    #print "I am getting the tick marks now"    
    #print "X extent: " + str(x_min)+ " " +str(x_max)
    #print "x ticks: "
    #print new_xlocs
    #print x_labels
   
    #return xlocs,ylocs,new_x_labels,new_y_labels
    return new_xlocs,x_labels


#==============================================================================
# Formats ticks for an imshow plot in UTM
# Filename is the name of the file with full path
# x_max, x_min, y_max, y_min are the extent of the plotting area (NOT the DEM)
# n_target ticks are the number of ticks for plotting
#------------------------------------------------------------------------------
def GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics):  
   
    
    CellSize,XMin,XMax,YMin,YMax = LSDMap_IO.GetUTMMaxMin(FileName)
    NDV, xsize, ysize, GeoT, Projection, DataType = LSDMap_IO.GetGeoInfo(FileName)    
   
    xmax_UTM = XMax
    xmin_UTM = XMin
      
    ymax_UTM = YMax
    ymin_UTM = YMin
    
    #print "now UTM, xmax: " +str(xmax_UTM)+" x_min: " +str(xmin_UTM)+" y_maxb: " +str(ymax_UTM)+" y_minb: " +str(ymin_UTM)
    
    dy_fig = ymax_UTM-ymin_UTM
    dx_fig = xmax_UTM-xmin_UTM
    
    dx_spacing = dx_fig/n_target_tics
    dy_spacing = dy_fig/n_target_tics
    
    if (dx_spacing>dy_spacing):
        dy_spacing = dx_spacing
    
    str_dy = str(dy_spacing)
    str_dy = str_dy.split('.')[0]
    n_digits = str_dy.__len__()
    nd = int(n_digits)
        
    first_digit = float(str_dy[0])
    
    #print "str_dy: " +str_dy+ " n_digits: " +str(nd)+" first_digit: " + str(first_digit)    
    
    dy_spacing_rounded = first_digit*pow(10,(nd-1))
    #print "n_digits: "+str(n_digits)+" dy_spacing: " +str(dy_spacing) + " and rounded: "+str(dy_spacing_rounded)
 
    str_xmin = str(xmin_UTM)
    str_ymin = str(ymin_UTM)
    #print "before split str_xmin: "+ str_xmin + " str ymin: " + str_ymin
    str_xmin = str_xmin.split('.')[0]
    str_ymin = str_ymin.split('.')[0]
    #print "after split str_xmin: "+ str_xmin + " str ymin: " + str_ymin
    xmin_UTM = float(str_xmin)
    ymin_UTM = float(str_ymin)
    #print "UTM: "+ str(xmin_UTM) + " str ymin: " + str(ymin_UTM)
    
    n_digx = str_xmin.__len__() 
    n_digy = str_ymin.__len__() 
    
    #print "n_dig_x: " + str(n_digx)+ " nd: " + str(nd)  
       
    if (n_digx-nd+1) >= 1:
        front_x = str_xmin[:(n_digx-nd+1)]
    else:
        front_x = str_xmin
        
    if (n_digy-nd+1) >= 1: 
        front_y = str_ymin[:(n_digy-nd+1)]
    else:
        front_y = str_ymin
    
    
      
    #print "xmin: " + str_xmin + " ymin: " + str_ymin + " n_digx: " + str(n_digx)+ " n_digy: " + str(n_digy)
    #print "frontx: " +front_x+" and fronty: "+ front_y
     
    round_xmin = float(front_x)*pow(10,nd-1)
    round_ymin = float(front_y)*pow(10,nd-1)
    
    #print "x_min: " +str(xmin_UTM)+ " round xmin: " +str(round_xmin)+ " y_min: " +str(ymin_UTM)+" round y_min: " + str(round_ymin)
    
    # now we need to figure out where the xllocs and ylocs are
    xUTMlocs = np.zeros(2*n_target_tics)
    yUTMlocs = np.zeros(2*n_target_tics)
    xlocs = np.zeros(2*n_target_tics)
    ylocs = np.zeros(2*n_target_tics)
    
    new_x_labels = []
    new_y_labels = []
    
    round_ymax = round_ymin+dy_spacing_rounded*(2*n_target_tics-1)
    #print "n_target_tics: " + str(n_target_tics) + " round_ymax: " +str(round_ymax)
    
    for i in range(0,2*n_target_tics):
        xUTMlocs[i] = round_xmin+(i)*dy_spacing_rounded
        yUTMlocs[i] = round_ymin+(i)*dy_spacing_rounded
        
                          
        #xlocs[i] = (xUTMlocs[i]-XMin)
        xlocs[i] = xUTMlocs[i]
        
        # need to account for the rows starting at the upper boundary
        ylocs[i] = round_ymax-(yUTMlocs[i]-round_ymin)
        
        #print "i: " + str(i) +" yUTM: " + str(yUTMlocs[i])+ " rounded, reversed: " +str( ylocs[i])
        
        
        new_x_labels.append( str(xUTMlocs[i]).split(".")[0] )
        new_y_labels.append( str(yUTMlocs[i]).split(".")[0] )

    #print xUTMlocs
    #print xlocs
    #print yUTMlocs
    #print ylocs
    #print new_x_labels
    #print new_y_labels

    new_xlocs = []
    new_xUTMlocs = []
    x_labels = []

    # Now loop through these to get rid of those not in range
    for index,xloc in enumerate(xlocs):
        #print xloc
        if (xloc < XMax and xloc > XMin):
            new_xlocs.append(xloc)
            new_xUTMlocs.append(xUTMlocs[index])
            x_labels.append(new_x_labels[index])
    
    new_ylocs = []
    new_yUTMlocs = []
    y_labels = []

    # Now loop through these to get rid of those not in range
    for index,yloc in enumerate(ylocs):
        #print yloc
        if (yloc < YMax and yloc > YMin):
            new_ylocs.append(yloc)
            new_yUTMlocs.append(yUTMlocs[index])
            y_labels.append(new_y_labels[index])    
            
    #print "======================================="
    #print "I am getting the tick marks now"    
    #print "X extent: " + str(XMin)+ " " +str(XMax)
    #print "Y extent: " + str(YMin)+ " " +str(YMax)
    #print "x ticks: "
    #print new_xlocs
    #print "y ticks: "
    #print new_ylocs
    
   
    #return xlocs,ylocs,new_x_labels,new_y_labels
    return new_xlocs,new_ylocs,x_labels,y_labels
#==============================================================================

#==============================================================================
def LogStretchDensityPlot(FileName, thiscmap='gray',colorbarlabel='Elevation in meters',clim_val = (0,0)):
    
    import matplotlib.pyplot as plt
    import matplotlib.lines as mpllines

    label_size = 20
    #title_size = 30
    axis_size = 28

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size 

    # get the data
    raster = LSDMap_IO.ReadRasterArrayBlocks(FileName)
    
    # get the log of the raster
    raster = np.log10(raster)
    
    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(FileName)
    
    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # make a figure, sized for a ppt slide
    fig = plt.figure(1, facecolor='white',figsize=(10,7.5))

    # make room for the colorbar
    fig.subplots_adjust(bottom=0.2)
    fig.subplots_adjust(top=0.9)
    #fig.subplots_adjust(left=0.2)
    #fig.subplots_adjust(right=0.8)
    
    ax1 =  fig.add_subplot(1,1,1)
    im = ax1.imshow(raster[::-1], thiscmap, extent = extent_raster)
    
    print "The is the extent raster data element"
    print extent_raster

    print "now I am in the mapping routine"
    print "x_min: " + str(x_min)
    print "x_max: " + str(x_max)
    print "y_min: " + str(y_min)
    print "y_max: " + str(y_max)

    # now get the tick marks    
    n_target_tics = 5
    xlocs,ylocs,new_x_labels,new_y_labels = GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics)  

    plt.xticks(xlocs, new_x_labels, rotation=60)  #[1:-1] skips ticks where we have no data
    plt.yticks(ylocs, new_y_labels) 
    
    print "The x locs are: " 
    print xlocs
    
    print "The x labels are: "
    print new_x_labels
    
    # some formatting to make some of the ticks point outward    
    for line in ax1.get_xticklines():
        line.set_marker(mpllines.TICKDOWN)
        #line.set_markeredgewidth(3)

    for line in ax1.get_yticklines():
        line.set_marker(mpllines.TICKLEFT)
        #line.set_markeredgewidth(3)  
    
    #plt.xlim(x_min,x_max)    
    #plt.ylim(y_max,y_min)   
   
    plt.xlabel('Easting (m)',fontsize = axis_size)
    plt.ylabel('Northing (m)', fontsize = axis_size)  

    ax1.set_xlabel("Easting (m)")
    ax1.set_ylabel("Northing (m)")
    
    # set the colour limits
    print "Setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1])
    if (clim_val == (0,0)):
        print "I don't think I should be here"
        im.set_clim(0, np.nanmax(raster))
    else:
        print "Now setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1])
        im.set_clim(clim_val[0],clim_val[1])
    
    
    cbar = fig.colorbar(im, orientation='horizontal')
    cbar.set_label(colorbarlabel)  
    
    #plt.tight_layout()

    plt.show()

#==============================================================================

#==============================================================================
def BasicDensityPlot(FileName, thiscmap='gray',colorbarlabel='Elevation in meters',clim_val = (0,0)):
    
    import matplotlib.pyplot as plt
    import matplotlib.lines as mpllines

    label_size = 20
    #title_size = 30
    axis_size = 28

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size 

    # get the data
    raster = LSDMap_IO.ReadRasterArrayBlocks(FileName)
    
    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(FileName)
    
    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # make a figure, sized for a ppt slide
    fig = plt.figure(1, facecolor='white',figsize=(10,7.5))

    # make room for the colorbar
    #fig.subplots_adjust(bottom=0.1)
    #fig.subplots_adjust(top=0.9)
    #fig.subplots_adjust(left=0.2)
    #fig.subplots_adjust(right=0.8)
    
    ax1 =  fig.add_subplot(1,1,1)
    im = ax1.imshow(raster[::-1], thiscmap, extent = extent_raster)
    
    print "The is the extent raster data element"
    print extent_raster

    print "now I am in the mapping routine"
    print "x_min: " + str(x_min)
    print "x_max: " + str(x_max)
    print "y_min: " + str(y_min)
    print "y_max: " + str(y_max)

    # now get the tick marks    
    n_target_tics = 5
    xlocs,ylocs,new_x_labels,new_y_labels = GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics)  

    plt.xticks(xlocs, new_x_labels, rotation=60)  #[1:-1] skips ticks where we have no data
    plt.yticks(ylocs, new_y_labels) 
    
    print "The x locs are: " 
    print xlocs
    
    print "The x labels are: "
    print new_x_labels
    
    # some formatting to make some of the ticks point outward    
    for line in ax1.get_xticklines():
        line.set_marker(mpllines.TICKDOWN)
        #line.set_markeredgewidth(3)

    for line in ax1.get_yticklines():
        line.set_marker(mpllines.TICKLEFT)
        #line.set_markeredgewidth(3)  
    
    plt.xlim(x_min,x_max)    
    plt.ylim(y_max,y_min)   
   
    plt.xlabel('Easting (m)',fontsize = axis_size)
    plt.ylabel('Northing (m)', fontsize = axis_size)  

    ax1.set_xlabel("Easting (m)")
    ax1.set_ylabel("Northing (m)")
    
    # set the colour limits
    print "Setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1])
    if (clim_val == (0,0)):
        print "I don't think I should be here"
        im.set_clim(0, np.nanmax(raster))
    else:
        print "Now setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1])
        im.set_clim(clim_val[0],clim_val[1])
    
    
    cbar = fig.colorbar(im, orientation='vertical')
    cbar.set_label(colorbarlabel)  
    
    #plt.tight_layout()

    plt.show()

#==============================================================================

#==============================================================================
def BasicDensityPlotGridPlot(FileName, thiscmap='gray',colorbarlabel='Elevation in meters',
                             clim_val = (0,0),FigFileName = 'Image.pdf', FigFormat = 'show'):
    

    print "======================================"
    print "Yo, I'm doing a draped plot"
    print FigFileName
    print FigFormat
    print "======================================"

    import matplotlib.pyplot as plt
    import matplotlib.lines as mpllines
    from mpl_toolkits.axes_grid1 import AxesGrid

    label_size = 20
    #title_size = 30
    axis_size = 28

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size 

    # get the data
    raster = LSDMap_IO.ReadRasterArrayBlocks(FileName)
    
    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(FileName)
    
    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # make a figure, sized for a ppt slide
    fig = plt.figure(1, facecolor='white',figsize=(10,7.5))

    gs = plt.GridSpec(100,75,bottom=0.1,left=0.1,right=0.9,top=1.0)
    ax = fig.add_subplot(gs[10:100,10:75])

    #grid = AxesGrid(fig, 111,  # similar to subplot(144)
    #                nrows_ncols=(1, 1),
    #                axes_pad=(0.8, 0.3),
    #                label_mode="1",
    #                share_all=True,
    #                cbar_location="right",
    #                cbar_mode="each",
    #                cbar_size="7%",
    #                cbar_pad="2%",
    #                )


    # now get the tick marks    
    n_target_tics = 5
    xlocs,ylocs,new_x_labels,new_y_labels = GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics)  


    #print "xmax: " + str(x_max)
    #print "xmin: " + str(x_min)
    #print "ymax: " + str(y_max)
    #print "ymin: " + str(y_min)


    im = ax.imshow(raster[::-1], thiscmap, extent = extent_raster, interpolation="nearest")
    #im = grid[0].imshow(raster, thiscmap, interpolation="nearest")

    cbar = plt.colorbar(im)
    cbar.set_label(colorbarlabel) 
    
    # set the colour limits
    #print "Setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1])
    if (clim_val == (0,0)):
        #print "I don't think I should be here"
        im.set_clim(0, np.nanmax(raster))
    else:
        print "Now setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1])
        im.set_clim(clim_val[0],clim_val[1])
    
    # go through the ticks     
    ax.spines['top'].set_linewidth(2.5)
    ax.spines['left'].set_linewidth(2.5)
    ax.spines['right'].set_linewidth(2.5)
    ax.spines['bottom'].set_linewidth(2.5) 
        
    # This affects all axes because we set share_all = True.
    ax.set_xlim(x_min,x_max)    
    ax.set_ylim(y_max,y_min)     

    ax.set_xticks(xlocs)
    ax.set_yticks(ylocs)
    
    ax.set_xticklabels(new_x_labels,rotation=60)
    ax.set_yticklabels(new_y_labels)  
    
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")  

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap        
    ax.tick_params(axis='both', width=2.5, pad = 10)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(10)        

    print "The figure format is: " + FigFormat
    if FigFormat == 'show':    
        plt.show()
    elif FigFormat == 'return':
        return fig
    else:
        plt.savefig(FigFileName,format=FigFormat)
        fig.clf()
    

#==============================================================================

#==============================================================================
def BasicDrapedPlotGridPlot(FileName, DrapeName, thiscmap='gray',drape_cmap='gray',
                            colorbarlabel='Elevation in meters',clim_val = (0,0),
                            drape_alpha = 0.6,FigFileName = 'Image.pdf',FigFormat = 'show'):
    
    
    print "======================================"
    print "Yo, I'm doing a draped plot"
    print FigFileName
    print FigFormat
    print "======================================"
    
    
    import matplotlib.pyplot as plt
    import matplotlib.lines as mpllines
    from mpl_toolkits.axes_grid1 import AxesGrid

    label_size = 20
    #title_size = 30
    axis_size = 28

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size 

    # get the data
    raster = LSDMap_IO.ReadRasterArrayBlocks(FileName)
    raster_drape = LSDMap_IO.ReadRasterArrayBlocks(DrapeName)
    
    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(FileName)
    
    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # make a figure, sized for a ppt slide
    fig = plt.figure(1, facecolor='white',figsize=(10,7.5))

    gs = plt.GridSpec(100,75,bottom=0.1,left=0.1,right=0.9,top=1.0)
    ax = fig.add_subplot(gs[10:100,10:75])

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
    xlocs,ylocs,new_x_labels,new_y_labels = GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics)  


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
    
    cbar = plt.colorbar(im1)
    cbar.set_label(colorbarlabel) 

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
    
    
    ax.spines['top'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5) 
     
    #ax.spines['bottom'].set_capstyle('projecting')

    #for spine in ax.spines.values():
    #    spine.set_capstyle('projecting')

    # This affects all axes because we set share_all = True.
    ax.set_xlim(x_min,x_max)    
    ax.set_ylim(y_max,y_min)     

    ax.set_xticks(xlocs)
    ax.set_yticks(ylocs)
    
    ax.set_xticklabels(new_x_labels,rotation=60)
    ax.set_yticklabels(new_y_labels)  
    
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")  

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap        
    ax.tick_params(axis='both', width=1.5, pad = 10)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(10)    

    print "The figure format is: " + FigFormat
    if FigFormat == 'show':    
        plt.show()
    elif FigFormat == 'return':
        return fig 
    else:
        plt.savefig(FigFileName,format=FigFormat,dpi=250)
        fig.clf()

#==============================================================================


#==============================================================================
def DrapedOverHillshade(FileName, DrapeName, thiscmap='gray',drape_cmap='gray',
                            colorbarlabel='Elevation in meters',clim_val = (0,0),
                            drape_alpha = 0.6, ShowColorbar = False, 
                            ShowDrapeColorbar=False, drape_cbarlabel=None):
    
    import matplotlib.pyplot as plt
    import matplotlib.lines as mpllines
    from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable

    label_size = 20
    #title_size = 30
    axis_size = 28

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size 

    hillshade = Hillshade(FileName) 
    #hillshade = LSDMap_IO.ReadRasterArrayBlocks(DrapeName)
    
    # DAV - option to supply array directly (after masking for example, rather
    # than reading directly from a file. Should not break anyone's code)
    # (You can't overload functions in Python...)
    if isinstance(DrapeName, str):
      raster_drape = LSDMap_IO.ReadRasterArrayBlocks(DrapeName)
    elif isinstance(DrapeName, np.ndarray):
      raster_drape = DrapeName
    else:
      print "DrapeName supplied is of type: ", type(DrapeName)
      raise ValueError('DrapeName must either be a string to a filename, \
      or a numpy ndarray type. Please try again.')
    
    
    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(FileName)
    
    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # make a figure, sized for a ppt slide
    fig = plt.figure(1, facecolor='white',figsize=(10,7.5))

    if ShowColorbar:    
        grid = AxesGrid(fig, 111, 
                        nrows_ncols=(1, 1),
                        axes_pad=(0.45, 0.15),
                        label_mode="1",
                        share_all=True,
                        cbar_location="right",
                        cbar_mode="each",
                        cbar_size="7%",
                        cbar_pad="2%",
                        )
                        
    #ShowDrapeColorbar = True                    
    #if ShowDrapeColorbar and ShowColorbar:
      #divider = make_axes_locatable(fig)
      #cax = divider.append_axes("right", size = "5%", pad=0.05)
      
    else:
        grid = AxesGrid(fig, 111, 
                        nrows_ncols=(1, 1),
                        axes_pad=(0.45, 0.15),
                        label_mode="1",
                        share_all=True,
                        )        


    # now get the tick marks    
    n_target_tics = 5
    xlocs,ylocs,new_x_labels,new_y_labels = GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics)  


    print "xmax: " + str(x_max)
    print "xmin: " + str(x_min)
    print "ymax: " + str(y_max)
    print "ymin: " + str(y_min)


    im = grid[0].imshow(hillshade[::-1], thiscmap, extent = extent_raster, interpolation="nearest")
    #im = grid[0].imshow(raster, thiscmap, interpolation="nearest")
    if ShowColorbar:
        cbar = grid.cbar_axes[0].colorbar(im)
        cbar.set_label_text(colorbarlabel)
    
    # set the colour limits
    print "Setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1])
    if (clim_val == (0,0)):
        print "I don't think I should be here"
        im.set_clim(0, np.nanmax(hillshade))
    else:
        print "Now setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1])
        im.set_clim(clim_val[0],clim_val[1])
    
    # Now for the drape: it is in grayscape
    im2 = grid[0].imshow(raster_drape[::-1], drape_cmap, extent = extent_raster, alpha = drape_alpha, interpolation="none")
      
    if ShowDrapeColorbar:
      cbar2 = grid.cbar_axes[0].colorbar(im2)
      cbar2.set_label_text(drape_cbarlabel)
      
      #plt.colorbar(im)

    # This affects all axes because we set share_all = True.
    grid.axes_llc.set_xlim(x_min,x_max)    
    grid.axes_llc.set_ylim(y_max,y_min)     

    grid.axes_llc.set_xticks(xlocs)
    grid.axes_llc.set_yticks(ylocs)
    
    grid.axes_llc.set_xticklabels(new_x_labels,rotation=60)
    grid.axes_llc.set_yticklabels(new_y_labels)  
    
    grid.axes_llc.set_xlabel("Easting (m)")
    grid.axes_llc.set_ylabel("Northing (m)")      

    plt.show()

#==============================================================================

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## This function plots the chi slope on a shaded relief map
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def DrapedOverFancyHillshade(FileName, HSName, DrapeName, thiscmap='gray',drape_cmap='gray',
                            colorbarlabel='Basin number',clim_val = (0,0),
                            drape_alpha = 0.6,FigFileName = 'Image.pdf',FigFormat = 'show',
                            elevation_threshold = 0):
    label_size = 10

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size 

    # get the data
    raster = LSDMap_IO.ReadRasterArrayBlocks(FileName)
    raster_HS = LSDMap_IO.ReadRasterArrayBlocks(HSName)
    raster_drape = LSDMap_IO.ReadRasterArrayBlocks(DrapeName)

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

    # now get the tick marks    
    n_target_tics = 5
    xlocs,ylocs,new_x_labels,new_y_labels = GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics)  
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
    HS_cmap = 'gray'
    im2 = ax.imshow(raster_HS[::-1], HS_cmap, extent = extent_raster, alpha = 0.4, interpolation="nearest")

    # Set the colour limits of the drape
    im2.set_clim(0,np.nanmax(raster_HS))

    # Now for the drape: it is in grayscale
    #print "drape_cmap is: "+drape_cmap
    im3 = ax.imshow(raster_drape[::-1], drape_cmap, extent = extent_raster, alpha = drape_alpha, interpolation="nearest")    

    # Set the colour limits of the drape
    im3.set_clim(0,np.nanmax(raster_drape))
        
    ax.spines['top'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1) 
   
    ax.set_xticklabels(new_x_labels,rotation=60)
    ax.set_yticklabels(new_y_labels)  
    
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")  

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap        
    ax.tick_params(axis='both', width=1, pad = 2)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(2)    

    # This affects all axes because we set share_all = True.
    ax.set_xlim(x_min,x_max)    
    ax.set_ylim(y_max,y_min)     

    ax.set_xticks(xlocs)
    ax.set_yticks(ylocs)   
 
    cbar = plt.colorbar(im3,cmap=drape_cmap,spacing='uniform', orientation='horizontal',cax=ax2)
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
## This function plots the chi slope on a shaded relief map
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def BasinsOverFancyHillshade(FileName, HSName, BasinName, Basin_csv_name, thiscmap='gray',drape_cmap='gray',
                             clim_val = (0,0), drape_alpha = 0.6,FigFileName = 'Image.pdf',
                             FigFormat = 'show',elevation_threshold = 0):

    label_size = 10

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size 

    # get the data
    raster = LSDMap_IO.ReadRasterArrayBlocks(FileName)
    raster_HS = LSDMap_IO.ReadRasterArrayBlocks(HSName)
    raster_drape = LSDMap_IO.ReadRasterArrayBlocks(BasinName)

    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(FileName)
    
    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]
    
    x_range = x_max-x_min
    y_range = y_max-y_min
    
    

    # make a figure, sized for a ppt slide
    fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))

    gs = plt.GridSpec(100,100,bottom=0.25,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])
    
    # now get the tick marks    
    n_target_tics = 5
    xlocs,ylocs,new_x_labels,new_y_labels = GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics)  
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
    HS_cmap = 'gray'
    im2 = ax.imshow(raster_HS[::-1], HS_cmap, extent = extent_raster, alpha = 0.4, interpolation="nearest")

    # Set the colour limits of the drape
    im2.set_clim(0,np.nanmax(raster_HS))

    # Now for the drape: it is in grayscale
    #print "drape_cmap is: "+drape_cmap
    im3 = ax.imshow(raster_drape[::-1], drape_cmap, extent = extent_raster, alpha = drape_alpha, interpolation="nearest")    

    # Set the colour limits of the drape
    im3.set_clim(0,np.nanmax(raster_drape))

    #========================================================
    # now we need to label the basins
    # Now we get the chi points
    EPSG_string = LSDMap_IO.GetUTMEPSG(FileName)
    print "EPSG string is: " + EPSG_string
    
    thisPointData = LSDMap_PD.LSDMap_PointData(Basin_csv_name) 
    
    # convert to easting and northing
    [easting,northing] = thisPointData.GetUTMEastingNorthingFromQuery(EPSG_string,"outlet_latitude","outlet_longitude")
 
    # The image is inverted so we have to invert the northing coordinate
    Ncoord = np.asarray(northing)
    Ncoord = np.subtract(extent_raster[3],Ncoord)
    Ncoord = np.add(Ncoord,extent_raster[2])
    
    these_data = thisPointData.QueryData("outlet_junction")
    #print M_chi
    these_data = [int(x) for x in these_data]

    # plot the centroids
    ax.scatter(easting,Ncoord,s=1, marker='o',color= 'r',alpha=0.7)
    
    # add a load of xmin, xmax points
    buffered_east = []
    buffered_north = []
    
    for loc in easting:
        buffered_east.append(loc)
        buffered_east.append(loc)
        buffered_east.append(loc)
        buffered_east.append(extent_raster[0])
        buffered_east.append(extent_raster[1])
    
    for loc in Ncoord:
        buffered_north.append(loc)
        buffered_north.append(extent_raster[2]) 
        buffered_north.append(extent_raster[3]) 
        buffered_north.append(loc) 
        buffered_north.append(loc) 
        
    
    
    # add text
    texts = []
    for index, datum in enumerate(these_data):
        this_easting = easting[index]
        this_northing = Ncoord[index]
        texts.append(ax.text(this_easting,this_northing, str(index),fontsize = 8, color= "r",alpha=0.7))
        #ax.annotate(str(datum), xy=(this_easting, this_northing),  xycoords='data',
        #    xytext=(this_easting+x_range/20, this_northing+y_range/20), textcoords='data', fontsize = 8,
        #    arrowprops=dict(facecolor='black', shrink=0.1,headlength = 2,
        #                    headwidth = 3,width=1),
        #    horizontalalignment='right', verticalalignment='top',
        #    )
    #iterations = adjust_text(texts,x=easting,y=Ncoord,autoalign='x',ax=ax,force_points=0.9)
    adjust_text(texts,x=buffered_east,y=buffered_north,autoalign='xy',ax=ax)
    #iterations = adjust_text(texts, )
    #plt.title(str(iterations)+ " iterations")
    

        
    ax.spines['top'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1) 
   
    ax.set_xticklabels(new_x_labels,rotation=60)
    ax.set_yticklabels(new_y_labels)  
    
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")  

    # This gets all the ticks, and pads them away from the axis so that the corners don't overlap        
    ax.tick_params(axis='both', width=1, pad = 2)
    for tick in ax.xaxis.get_major_ticks():
        tick.set_pad(2)    

    # This affects all axes because we set share_all = True.
    ax.set_xlim(x_min,x_max)    
    ax.set_ylim(y_max,y_min)     

    ax.set_xticks(xlocs)
    ax.set_yticks(ylocs)   
    
    print "The figure format is: " + FigFormat
    if FigFormat == 'show':    
        plt.show()
    elif FigFormat == 'return':
        return fig 
    else:
        plt.savefig(FigFileName,format=FigFormat,dpi=500)
        fig.clf()


#==============================================================================
# Make a simple hillshade plot
def Hillshade(raster_file, azimuth = 315, angle_altitude = 45, NoDataValue = -9999): 
    
    array = LSDMap_IO.ReadRasterArrayBlocks(raster_file,raster_band=1)    
    
    # DAV attempting mask nodata vals
    nodata_mask = array == NoDataValue
    array[nodata_mask] = np.nan
    
    x, y = np.gradient(array)
    slope = np.pi/2. - np.arctan(np.sqrt(x*x + y*y))
    aspect = np.arctan2(-x, y)
    azimuthrad = azimuth*np.pi / 180.
    altituderad = angle_altitude*np.pi / 180.
     
 
    shaded = np.sin(altituderad) * np.sin(slope)\
     + np.cos(altituderad) * np.cos(slope)\
     * np.cos(azimuthrad - aspect)
     
     
     
    #this_array = 255*(shaded + 1)/2 
    return 255*(shaded + 1)/2
#==============================================================================


#==============================================================================
def SwathPlot(path, filename, axis):

    # get the path to the raster file
    NewPath = LSDOst.AppendSepToDirectoryPath(path)    
    FileName = NewPath+filename
        
    # get the data vectors
    means,medians,std_deviations,twentyfifth_percentile,seventyfifth_percentile = LSDMap_BM.SimpleSwath(path, filename, axis)
    
    print "Means shape is: "
    print means.shape    
    
    x_vec,y_vec = LSDMap_IO.GetLocationVectors(FileName)
    
    
    print "X shape is: "
    print x_vec.shape
    
    print "Y shape is: "
    print y_vec.shape
    
    import matplotlib.pyplot as plt
    import matplotlib.lines as mpllines
    from mpl_toolkits.axes_grid1 import AxesGrid

    label_size = 20
    #title_size = 30
    axis_size = 28

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure, sized for a ppt slide
    fig = plt.figure(1, facecolor='white',figsize=(10,7.5)) 

    gs = plt.GridSpec(100,75,bottom=0.1,left=0.1,right=0.9,top=1.0)
    ax = fig.add_subplot(gs[10:100,10:75])
    
    if axis == 0:
        dir_vec = x_vec
    else:
        dir_vec = y_vec
        
    min_sd = np.subtract(means,std_deviations)
    plus_sd = np.add(means,std_deviations) 
        
    ax.plot(dir_vec,means, linewidth = 2, color = "red")
    #ax.fill_between(dir_vec, twentyfifth_percentile, seventyfifth_percentile, facecolor='green', alpha = 0.7, interpolate=True)
    ax.fill_between(dir_vec, min_sd, plus_sd, facecolor='blue', alpha = 0.5, interpolate=True)  
    
    ax.set_xlim(dir_vec[0],dir_vec[-1])
    
    plt.show()
#==============================================================================    
    


#==============================================================================


#==============================================================================
def round_to_n(x, n):
    if n < 1:
        raise ValueError("number of significant digits must be >= 1")
    # Use %e format to get the n most significant digits, as a string.
    format = "%." + str(n-1) + "e"
    as_string = format % x
    return float(as_string)
#==============================================================================



def init_plotting_DV():
    plt.rcParams['figure.figsize'] = (8, 8)
    plt.rcParams['font.size'] = 17
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['axes.labelsize'] = 1.2*plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = 1.2*plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['savefig.dpi'] = 2*plt.rcParams['savefig.dpi']
    plt.rcParams['xtick.major.size'] = 3
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.major.size'] = 3
    plt.rcParams['ytick.minor.size'] = 3
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['ytick.minor.width'] = 1
    plt.rcParams['legend.frameon'] = True
    plt.rcParams['legend.loc'] = 'center left'
    plt.rcParams['axes.linewidth'] = 1
    plt.rcParams['xtick.minor.visible'] = False
    plt.rcParams['ytick.minor.visible'] = False
    plt.rcParams['lines.linewidth'] = 2