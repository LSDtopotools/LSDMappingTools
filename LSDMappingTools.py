## raster_tools.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## These functions are tools to deal with rasters
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## SMM
## 26/07/2014
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

import osgeo.gdal as gdal
import numpy as np
import numpy.ma as ma
from osgeo import osr
from os.path import exists
from osgeo.gdalconst import GA_ReadOnly
from numpy import uint8
from matplotlib import rcParams

#==============================================================================
def GetUTMMaxMin(FileName):

    if exists(FileName) is False:
            raise Exception('[Errno 2] No such file or directory: \'' + FileName + '\'')    
    
    NDV, xsize, ysize, GeoT, Projection, DataType = GetGeoInfo(FileName)
    CellSize = GeoT[1]
    XMin = GeoT[0]
    XMax = XMin+CellSize*xsize

    YMax = GeoT[3]
    YMin = YMax-CellSize*ysize
    
    return CellSize,XMin,XMax,YMin,YMax
#==============================================================================    

#==============================================================================
# This gets the extent of the raster
def GetRasterExtent(FileName):
    
    CellSize,XMin,XMax,YMin,YMax = GetUTMMaxMin(FileName)
    extent = [XMin,XMax,YMin,YMax]
    return extent    

#==============================================================================
# Function to read the original file's projection:
def GetGeoInfo(FileName):

    if exists(FileName) is False:
            raise Exception('[Errno 2] No such file or directory: \'' + FileName + '\'')    
    
    
    SourceDS = gdal.Open(FileName, gdal.GA_ReadOnly)
    if SourceDS == None:
        raise Exception("Unable to read the data file")
    
    NDV = SourceDS.GetRasterBand(1).GetNoDataValue()
    xsize = SourceDS.RasterXSize
    ysize = SourceDS.RasterYSize
    GeoT = SourceDS.GetGeoTransform()
    Projection = osr.SpatialReference()
    Projection.ImportFromWkt(SourceDS.GetProjectionRef())
    DataType = SourceDS.GetRasterBand(1).DataType
    DataType = gdal.GetDataTypeName(DataType)
    
    return NDV, xsize, ysize, GeoT, Projection, DataType
#==============================================================================

#==============================================================================
def ReadRasterArrayBlocks(raster_file,raster_band=1):
    
    if exists(raster_file) is False:
            raise Exception('[Errno 2] No such file or directory: \'' + raster_file + '\'')    
    
    dataset = gdal.Open(raster_file, GA_ReadOnly )
    if dataset == None:
        raise Exception("Unable to read the data file")
    
    band = dataset.GetRasterBand(raster_band)

    block_sizes = band.GetBlockSize()
    x_block_size = block_sizes[0]
    y_block_size = block_sizes[1]

    #If the block y size is 1, as in a GeoTIFF image, the gradient can't be calculated, 
    #so more than one block is used. In this case, using8 lines gives a similar 
    #result as taking the whole array.
    if y_block_size < 8:
        y_block_size = 8

    xsize = band.XSize
    ysize = band.YSize
    
    print "xsize: " +str(xsize)+" and y size: " + str(ysize)

    max_value = band.GetMaximum()
    min_value = band.GetMinimum()
    
    # now initiate the array
    data_array = np.zeros((ysize,xsize))
    
    #print "data shape is: " 
    #print data_array.shape

    if max_value == None or min_value == None:
        stats = band.GetStatistics(0, 1)
        max_value = stats[1]
        min_value = stats[0]
        
    for i in range(0, ysize, y_block_size):
        if i + y_block_size < ysize:
            rows = y_block_size
        else:
            rows = ysize - i
        
        for j in range(0, xsize, x_block_size):
            if j + x_block_size < xsize:
                cols = x_block_size
            else:
                cols = xsize - j
            
            # get the values for this block
            values = band.ReadAsArray(j, i, cols, rows)
            
            # move these values to the data array
            data_array[i:i+rows,j:j+cols] = values
            
    return data_array
#==============================================================================

#==============================================================================
# Formats ticks for an imshow plot in UTM
# Filename is the name of the file with full path
# x_max, x_min, y_max, y_min are the extent of the plotting area (NOT the DEM)
# n_target ticks are the number of ticks for plotting
#------------------------------------------------------------------------------
def GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics):  
   
    
    CellSize,XMin,XMax,YMin,YMax = GetUTMMaxMin(FileName)
    NDV, xsize, ysize, GeoT, Projection, DataType = GetGeoInfo(FileName)    
    
    xmax_UTM = XMin+XMax*CellSize
    xmin_UTM = XMin+XMin*CellSize
      
    # need to be careful with the ymax_UTM since the rows go from the top
    # but the header index is to bottom corner    
    
    print "yll: "+str(YMin)+" and nrows: " +str(ysize) + " dx: "+str(CellSize)   
    
    ymax_from_bottom = ysize-YMin
    ymin_from_bottom = ysize-YMax
    ymax_UTM = YMin+ymax_from_bottom*CellSize
    ymin_UTM = YMin+ymin_from_bottom*CellSize
    
    print "now UTM, xmax: " +str(xmax_UTM)+" x_min: " +str(xmin_UTM)+" y_maxb: " +str(ymax_UTM)+" y_minb: " +str(ymin_UTM)
    
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
    
    print "str_dy: " +str_dy+ " n_digits: " +str(nd)+" first_digit: " + str(first_digit)    
    
    dy_spacing_rounded = first_digit*pow(10,(nd-1))
    print "n_digits: "+str(n_digits)+" dy_spacing: " +str(dy_spacing) + " and rounded: "+str(dy_spacing_rounded)
 
    str_xmin = str(xmin_UTM)
    str_ymin = str(ymin_UTM)
    str_xmin = str_xmin.split('.')[0]
    str_ymin = str_ymin.split('.')[0]
    xmin_UTM = float(str_xmin)
    ymin_UTM = float(str_ymin)
    
    n_digx = str_xmin.__len__() 
    n_digy = str_ymin.__len__() 
    
    front_x = str_xmin[:(n_digx-nd+1)]
    front_y = str_ymin[:(n_digy-nd+1)]
      
    print "xmin: " + str_xmin + " ymin: " + str_ymin + " n_digx: " + str(n_digx)+ " n_digy: " + str(n_digy)
    print "frontx: " +front_x+" and fronty: "+ front_y
     
    round_xmin = float(front_x)*pow(10,nd-1)
    round_ymin = float(front_y)*pow(10,nd-1)
    
    print "x_min: " +str(xmin_UTM)+ " round xmin: " +str(round_xmin)+ " y_min: " +str(ymin_UTM)+" round y_min: " + str(round_ymin)
    
    # now we need to figure out where the xllocs and ylocs are
    xUTMlocs = np.zeros(2*n_target_tics)
    yUTMlocs = np.zeros(2*n_target_tics)
    xlocs = np.zeros(2*n_target_tics)
    ylocs = np.zeros(2*n_target_tics)
    
    new_x_labels = []
    new_y_labels = []
    
    for i in range(0,2*n_target_tics):
        xUTMlocs[i] = round_xmin+(i)*dy_spacing_rounded
        yUTMlocs[i] = round_ymin+(i)*dy_spacing_rounded
                  
        xlocs[i] = (xUTMlocs[i]-XMin)/CellSize
        
        # need to account for the rows starting at the upper boundary
        ylocs[i] = ysize-((yUTMlocs[i]-YMin)/CellSize)
        
        new_x_labels.append( str(xUTMlocs[i]).split(".")[0] )
        new_y_labels.append( str(yUTMlocs[i]).split(".")[0] )
   
    return xlocs,ylocs,new_x_labels,new_y_labels
#==============================================================================



#==============================================================================
def BasicDensityPlot(FileName):
    
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
    raster = ReadRasterArrayBlocks(FileName)
    
    # now get the extent
    extent_raster = GetRasterExtent(FileName)
    
    x_min = extent_raster[0]
    x_max = extent_raster[1]
    y_min = extent_raster[2]
    y_max = extent_raster[3]

    # make a figure, sized for a ppt slide
    fig = plt.figure(1, facecolor='white',figsize=(10,7.5))
    ax1 =  fig.add_subplot(1,1,1)
    im = ax1.imshow(raster, cmap='gray', extent = extent_raster)

    # now get the tick marks    
    n_target_tics = 5
    xlocs,ylocs,new_x_labels,new_y_labels = GetTicksForUTM(FileName,x_max,x_min,y_max,y_min,n_target_tics)  
    plt.xticks(xlocs, new_x_labels, rotation=60)  #[1:-1] skips ticks where we have no data
    plt.yticks(ylocs, new_y_labels) 
    
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
    im.set_clim(0, np.max(raster))
    cbar = fig.colorbar(im, orientation='horizontal')
    cbar.set_label("Elevation in meters")  

    plt.show()

#==============================================================================








def round_to_n(x, n):
    if n < 1:
        raise ValueError("number of significant digits must be >= 1")
    # Use %e format to get the n most significant digits, as a string.
    format = "%." + str(n-1) + "e"
    as_string = format % x
    return float(as_string)


def read_headers(input_file):

    with open(input_file+'.hdr','r') as f:   
        return [float(h) if not h.isalpha() else h for h in [l.split()[1] for l in f.readlines()]]  #isdigit() does not catch floats      

def read_bin(filename):
    import sys
    import numpy as np

    with open(filename + '.flt', "rb") as f:
        raster_data = np.fromstring(f.read(), 'f')

    if sys.byteorder == 'big':
        raster_data = raster_data.byteswap()  #ensures data is little endian

    return raster_data
    
def read_flt(input_file):

    if input_file.endswith('.flt') or input_file.endswith('.hdr'):
        input_file = input_file[:-4]    
    else:
        print 'Incorrect filename'
        return 0,0 #exits module gracefully
    
    headers = read_headers(input_file)
    
    #read the data as a 1D array and reshape it to the dimensions in the header
    raster_array = read_bin(input_file).reshape(headers[1], headers[0]) 
    raster_array = raster_array.reshape(headers[1], headers[0]) #rows, columns

    return raster_array, headers

def read_ascii_raster(ascii_raster_file):
    import numpy as np
    
    with open(ascii_raster_file) as f:
        header_data = [float(f.next().split()[1]) for x in xrange(6)] #read the first 6 lines
         
    raster_data = np.genfromtxt(ascii_raster_file, delimiter=' ', skip_header=6)
    raster_data = raster_data.reshape(header_data[1], header_data[0]) #rows, columns
    
    return raster_data, header_data

# this gets the extent of the asc for use with plotting
# It returns a list with 4 elements, x_min, x_max, y_min,y_max
def get_raster_extent_asc(header):


  x_min = header[2]
  y_min = header[3]
  spacing = header[4]
  n_cols = header[0]
  n_rows = header[1]

  x_max = x_min+n_cols*spacing
  y_max = y_min+n_rows*spacing

  extent = [x_min,x_max,y_min,y_max]
  return extent
  

# This function makes a simple density plot of a raster
def simple_density_plot_asc(rfname):

  import numpy as np, matplotlib.pyplot as plt
  from matplotlib import rcParams
  import matplotlib.colors as colors
  import matplotlib.cm as cmx

  label_size = 20
  #title_size = 30
  axis_size = 28

  # Set up fonts for plots
  rcParams['font.family'] = 'sans-serif'
  rcParams['font.sans-serif'] = ['arial']
  rcParams['font.size'] = label_size 

  # get the data
  raster,header = read_ascii_raster(rfname)

  # now get the extent
  extent_raster = get_raster_extent_asc(header)

  #print extent_raster

  # make a figure, sized for a ppt slide
  fig = plt.figure(1, facecolor='white',figsize=(10,7.5))
  ax1 =  fig.add_subplot(1,1,1)
  im = ax1.imshow(raster, cmap='gray', extent = extent_raster)
  ax1.set_xlabel("Easting (m)")
  ax1.set_ylabel("Northing (m)")
  im.set_clim(0, np.max(raster))
  cbar = fig.colorbar(im, orientation='horizontal')
  cbar.set_label("Elevation in meters")  

  plt.show()

def format_ticks_for_UTM_imshow(hillshade_header,x_max,x_min,y_max,y_min,n_target_tics):
    import numpy as np    
   
    xmax_UTM = hillshade_header[2]+x_max*hillshade_header[4]
    xmin_UTM = hillshade_header[2]+x_min*hillshade_header[4]
      
    # need to be careful with the ymax_UTM since the rows go from the top
    # but the header index is to bottom corner    
    
    print "yll: "+str(hillshade_header[3])+" and nrows: " +str(hillshade_header[1]) + " dx: "+str(hillshade_header[4])   
    
    ymax_from_bottom = hillshade_header[1]-y_min
    ymin_from_bottom = hillshade_header[1]-y_max
    ymax_UTM = hillshade_header[3]+ymax_from_bottom*hillshade_header[4]
    ymin_UTM = hillshade_header[3]+ymin_from_bottom*hillshade_header[4]
    
    print "now UTM, xmax: " +str(xmax_UTM)+" x_min: " +str(xmin_UTM)+" y_maxb: " +str(ymax_UTM)+" y_minb: " +str(ymin_UTM)
    
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
    
    print "str_dy: " +str_dy+ " n_digits: " +str(nd)+" first_digit: " + str(first_digit)    
    
    dy_spacing_rounded = first_digit*pow(10,(nd-1))
    print "n_digits: "+str(n_digits)+" dy_spacing: " +str(dy_spacing) + " and rounded: "+str(dy_spacing_rounded)
 
    str_xmin = str(xmin_UTM)
    str_ymin = str(ymin_UTM)
    str_xmin = str_xmin.split('.')[0]
    str_ymin = str_ymin.split('.')[0]
    xmin_UTM = float(str_xmin)
    ymin_UTM = float(str_ymin)
    
    n_digx = str_xmin.__len__() 
    n_digy = str_ymin.__len__() 
    
    front_x = str_xmin[:(n_digx-nd+1)]
    front_y = str_ymin[:(n_digy-nd+1)]
      
    print "xmin: " + str_xmin + " ymin: " + str_ymin + " n_digx: " + str(n_digx)+ " n_digy: " + str(n_digy)
    print "frontx: " +front_x+" and fronty: "+ front_y
     
    round_xmin = float(front_x)*pow(10,nd-1)
    round_ymin = float(front_y)*pow(10,nd-1)
    
    print "x_min: " +str(xmin_UTM)+ " round xmin: " +str(round_xmin)+ " y_min: " +str(ymin_UTM)+" round y_min: " + str(round_ymin)
    
    # now we need to figure out where the xllocs and ylocs are
    xUTMlocs = np.zeros(2*n_target_tics)
    yUTMlocs = np.zeros(2*n_target_tics)
    xlocs = np.zeros(2*n_target_tics)
    ylocs = np.zeros(2*n_target_tics)
    
    new_x_labels = []
    new_y_labels = []
    
    for i in range(0,2*n_target_tics):
        xUTMlocs[i] = round_xmin+(i)*dy_spacing_rounded
        yUTMlocs[i] = round_ymin+(i)*dy_spacing_rounded
                  
        xlocs[i] = (xUTMlocs[i]-hillshade_header[2])/hillshade_header[4]
        
        # need to account for the rows starting at the upper boundary
        ylocs[i] = hillshade_header[1]-((yUTMlocs[i]-hillshade_header[3])/hillshade_header[4])
        
        new_x_labels.append( str(xUTMlocs[i]).split(".")[0] )
        new_y_labels.append( str(yUTMlocs[i]).split(".")[0] )
   
    return xlocs,ylocs,new_x_labels,new_y_labels




def vectorize(hillshade_file, m_value_file):
    
    import matplotlib.pyplot as pp
    import numpy as np
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    from matplotlib import rcParams
    
    #get data
    hillshade, hillshade_header = read_flt(hillshade_file)
    m_values, m_values_header = read_flt(m_value_file)
    
    #handle plotting hillshades which are larger than the m_value raster
    #cannot cope with m_value raster larger than the hillshade
    corrected_x = 0    
    corrected_y = 0
    if (hillshade_header[0] != m_values_header[0]) or (hillshade_header[1] != m_values_header[1]):
         corrected_x = (m_values_header[2] - hillshade_header[2]) / hillshade_header[4]
         corrected_y = (((m_values_header[3] / m_values_header[4] )+ m_values_header[1]) - ((hillshade_header[3] / hillshade_header[4]) + hillshade_header[1])) * -1
    
    #ignore nodata values    
    hillshade = np.ma.masked_where(hillshade == -9999, hillshade)    
    m_values = np.ma.masked_where(m_values == -9999, m_values)
    
    #fonts
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = 12  
    
    fig, ax = pp.subplots()
    
    ax.imshow(hillshade, vmin=0, vmax=255, cmap=cmx.gray)
          
    xlocs, xlabels = pp.xticks()
    ylocs, ylabels = pp.yticks()
   
    new_x_labels = np.linspace(hillshade_header[2],hillshade_header[2]+(hillshade_header[1]*hillshade_header[4]), len(xlocs))
    new_y_labels = np.linspace(hillshade_header[3],hillshade_header[3]+(hillshade_header[0]*hillshade_header[4]), len(ylocs))        
    
    new_x_labels = [str(x).split('.')[0] for x in new_x_labels] #get rid of decimal places in axis ticks
    new_y_labels = [str(y).split('.')[0] for y in new_y_labels][::-1] #invert y axis
    pp.xticks(xlocs[1:-1], new_x_labels[1:-1], rotation=30)  #[1:-1] skips ticks where we have no data
    pp.yticks(ylocs[1:-1], new_y_labels[1:-1])
    
    pp.xlabel('Easting (m)')
    pp.ylabel('Northing (m)')    
    
    # SET UP COLOURMAPS
    jet = pp.get_cmap('jet')
    
    m_MIN = np.min(m_values)
    m_MAX = np.max(m_values)
    cNorm_m_values  = colors.Normalize(vmin=m_MIN, vmax=m_MAX)
    scalarMap_m_values = cmx.ScalarMappable(norm=cNorm_m_values, cmap=jet)    
    
    for i in xrange(len(m_values)):
        for j in xrange(len(m_values[0])):
            if m_values[i][j] > 0:
                colorVal = scalarMap_m_values.to_rgba(m_values[i][j])
                pp.scatter(j + corrected_x, i + corrected_y, marker=".", color=colorVal,edgecolors='none')               
                  
    # Configure final plot
    sm = pp.cm.ScalarMappable(cmap=jet,norm=pp.normalize(vmin=m_MIN, vmax=m_MAX))
    sm._A = []
    cbar = pp.colorbar(sm)
    cbar.set_label('M Values')
    
    pp.show()
    
