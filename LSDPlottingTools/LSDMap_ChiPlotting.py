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
from matplotlib import rcParams
import LSDMap_GDALIO as LSDMap_IO
import LSDMap_BasicManipulation as LSDMap_BM
import LSDOSystemTools as LSDOst
import LSDMap_BasicPlotting as LSDMap_BP
import LSDMap_PointData as LSDMap_PD

def BasicChiPlotGridPlot(FileName, DrapeName, chi_csv_fname, thiscmap='gray',drape_cmap='gray',
                            colorbarlabel='Elevation in meters',clim_val = (0,0),
                            drape_alpha = 0.6,FigFileName = 'Image.pdf',FigFormat = 'show'):
    
    import matplotlib.pyplot as plt
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
    
    # now get the extent
    extent_raster = LSDMap_IO.GetRasterExtent(FileName)
    
    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # make a figure, sized for a ppt slide
    fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))

    gs = plt.GridSpec(100,100,bottom=0.3,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:100])
    
    # This is the axis for the colorbar
    ax2 = fig.add_subplot(gs[5:10,10:60])

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
        im1.set_clim(0, np.max(raster))
    else:
        print "Now setting colour limits to "+str(clim_val[0])+" and "+str(clim_val[1])
        im1.set_clim(clim_val[0],clim_val[1])
   
    plt.hold(True)
   
    # Now for the drape: it is in grayscale
    #print "drape_cmap is: "+drape_cmap
    im3 = ax.imshow(raster_drape[::-1], drape_cmap, extent = extent_raster, alpha = drape_alpha, interpolation="nearest")

    # Set the colour limits of the drape
    im3.set_clim(0,np.max(raster_drape))
    
    
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
    
    # convert to easting and northing
    [easting,northing] = thisPointData.GetUTMEastingNorthing(EPSG_string)
    
    
    # The image is inverted so we have to invert the northing coordinate
    Ncoord = np.asarray(northing)
    Ncoord = np.subtract(extent_raster[3],Ncoord)
    Ncoord = np.add(Ncoord,extent_raster[2])
    
    M_chi = thisPointData.QueryData('m_chi')
    M_chi = [float(x) for x in M_chi]
    #print M_chi
    
    # make a color map of fixed colors
    this_cmap = colors.ListedColormap(['#2c7bb6','#abd9e9','#ffffbf','#fdae61','#d7191c'])
    bounds=[0,50,100,175,250,1205]
    norm = colors.BoundaryNorm(bounds, this_cmap.N)
    
    sc = ax.scatter(easting,Ncoord,s=1, c=M_chi,cmap=this_cmap,norm=norm,edgecolors='none')

    # This affects all axes because we set share_all = True.
    ax.set_xlim(x_min,x_max)    
    ax.set_ylim(y_max,y_min)     

    ax.set_xticks(xlocs)
    ax.set_yticks(ylocs)   
    
    cbar = plt.colorbar(sc,cmap=this_cmap,norm=norm,spacing='uniform', ticks=bounds, boundaries=bounds,orientation='horizontal',cax=ax2)
    cbar.set_label(colorbarlabel, fontsize=10)
    
    print "The figure format is: " + FigFormat
    if FigFormat == 'show':    
        plt.show()
    elif FigFormat == 'return':
        return fig 
    else:
        plt.savefig(FigFileName,format=FigFormat,dpi=250)
        fig.clf()