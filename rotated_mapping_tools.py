# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#import modules
import numpy as np
from shapely.geometry import shape, MultiPolygon, Polygon
import numpy.ma as ma

# import the basemap library
from mpl_toolkits.basemap import Basemap
from osgeo import gdal
import fiona
from pyproj import Proj, transform

# import plotting tools and set the back end for running on server
import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams
import matplotlib.pyplot as plt
#plt.style.use('ggplot')

## modules
from shapely.ops import cascaded_union
from descartes import PolygonPatch

def HaversineDistance(lon1,lat1,lon2,lat2):
    """
    Function to calculate the great circle distance between two points
    using the Haversine formula
    """
    R = 6371. #Mean radius of the Earth

    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2.)**2. + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.)**2.
    c = 2.*np.arcsin(np.sqrt(a))
    distance = R * c

    return distance

def CreateFigure(FigSizeFormat="default", AspectRatio=16./9.):
    """
    This function creates a default matplotlib figure object

    Args:
        FigSizeFormat: the figure size format according to journal for which the figure is intended
            values are geomorphology,ESURF, ESPL, EPSL, JGR, big
            default is ESURF

        AspectRatio: The shape of the figure determined by the aspect ratio, default is 16./9.

    Returns:
        matplotlib figure object

    Author: MDH
    """
    # set figure sizes (in inches) based on format
    if FigSizeFormat == "geomorphology":
        FigWidth_Inches = 6.25
    elif FigSizeFormat == "big":
        FigWidth_Inches = 16
    elif FigSizeFormat == "ESURF":
        FigWidth_Inches = 4.92
    elif FigSizeFormat == "ESPL":
        FigWidth_Inches = 7.08
    elif FigSizeFormat == "EPSL":
        FigWidth_Inches = 7.48
    elif FigSizeFormat == "JGR":
        FigWidth_Inches = 6.6

    else:
        FigWidth_Inches = 4.92126

    MinimumAspectRatio=0.5
    MaximumAspectRatio=2.

    if AspectRatio < MinimumAspectRatio:
        AspectRatio = MinimumAspectRatio
    elif AspectRatio > MaximumAspectRatio:
        AspectRatio = MaximumAspectRatio

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = 10
    rcParams['text.usetex'] = False

    Fig = plt.figure(figsize=(FigWidth_Inches,FigWidth_Inches/AspectRatio))

    return Fig



def CreateMapFigure(BasemapExtentDataset, AspectRatio=0, FigSizeFormat="default", Rotation=0):
    """
    This function creates a default matplotlib basemap figure object

    Args:
        BasemapExtentDataset: a raster or shapefile dataset on which to base the figure extent
            raster must be a geotif in lat long (for now) or a shapefile
        AspectRatio: the aspect ratio of the resulting plot. If not specified it will be set by the
            rectangular extent of the dataset (convex hull)
        FigSizeFormat: the figure size format according to journal for which the figure is intended
            values are geomorphology,ESURF, ESPL, EPSL, JGR, big
            default is ESURF
        Rotation: The amount to rotate the lat-long coordinate system by in order to produce maps
            that maximise the use of rectangular space. Default is zero which will result in no rotation?

    Returns:
        matplotlib figure object

    Author: MDH
    """

    print ("Creating map figure...")
    print(BasemapExtentDataset.split(".")[-1])

    #Check data format is tif or shapefile and open, transform coordinates if necessary
    if (BasemapExtentDataset.split(".")[-1] == "tif"):
        #load the dataset
        Dataset = gdal.Open(BasemapExtentDataset, gdal.GA_ReadOnly)
        Geo = Dataset.GetGeoTransform()
        #get dataset extent
        xres = Geo[1]
        yres = Geo[5]
        xmin = Geo[0] + xres * 0.5
        xmax = Geo[0] + (xres * Dataset.RasterXSize) - xres * 0.5
        ymin = Geo[3] + (yres * Dataset.RasterYSize) + yres * 0.5
        ymax = Geo[3] - yres * 0.5

    #Check data format is tif or shapefile and open, transform coordinates if necessary
    elif (BasemapExtentDataset.split(".")[-1] == "bil"):
        #load the dataset
        Dataset = gdal.Open(BasemapExtentDataset, gdal.GA_ReadOnly)
        Geo = Dataset.GetGeoTransform()
        #get dataset extent
        xres = Geo[1]
        yres = Geo[5]
        xmin = Geo[0] + xres * 0.5
        xmax = Geo[0] + (xres * Dataset.RasterXSize) - xres * 0.5
        ymin = Geo[3] + (yres * Dataset.RasterYSize) + yres * 0.5
        ymax = Geo[3] - yres * 0.5

    elif(BasemapExtentDataset.split(".")[-1] == "shp"):

        # Read polygons from shapefile in projected coordinates
        Polygons = ConvertShapefile2LatLong(BasemapExtentDataset)

        # Combine intoa multipolygon and get the minimum rotated rectangle
        PolygonList = [Poly for Key,Poly in Polygons.iteritems()]
        MP = cascaded_union(PolygonList)
        MinRotatedRectangle = MP.minimum_rotated_rectangle
        [xmin,ymin,xmax,ymax] = MinRotatedRectangle.bounds

        # use the minimum rotated rectangle to set the plotting extent
        lon,lat = MinRotatedRectangle.exterior.coords.xy
        Length = HaversineDistance(lon[0],lat[0],lon[1],lat[1])
        Width = HaversineDistance(lon[-1],lat[-1],lon[-2],lat[-2])
        Offset = HaversineDistance(lon[0],lat[0],lon[0],lat[1])
        #Rotation = np.abs(np.degrees(np.arcsin(Offset/Length)))
    else:
        print "Can't handle file format"
        return

    if (Rotation > 45):
        Rotation = -(90-Rotation)
        llx = lon[1]
        lly = lat[1]
        urx = lon[3]
        ury = lat[3]
        if (AspectRatio == 0):
            AspectRatio = Width/Length

    else:
        Pad = 1
        llx = lon[0]
        lly = lat[0]
        urx = lon[2]
        ury = lat[2]
        if (AspectRatio == 0):
            AspectRatio = Length/Width

    #create the figure
    Fig = CreateFigure(AspectRatio=AspectRatio,FigSizeFormat="JGR")
    Ax = Fig.add_axes([0.15,0.15,0.75,0.75])

    #create the map object
    # lat_0 needs to be set to 90 in order for rotation values to be the true rotation of the plot.
    # use a lambert conformal conic projection
    latpad=np.abs(lly-ury)*0.1
    lonpad=np.abs(llx-urx)*0.1
    Map = Basemap(llcrnrlon=llx-lonpad,urcrnrlon=urx+lonpad, llcrnrlat=lly-latpad,  urcrnrlat=ury+latpad,
            projection='lcc', resolution='l', lat_0=90, lon_0=llx-Rotation)

#    # convert to map coordinates
#    RectDict = {}
#    RectDict[0] = MinRotatedRectangle
#    Rect = ConvertPolygonsLatLong2MapCoords(RectDict, Map)
#    #plot the rectangle extent
#    Patch = PolygonPatch(Rect[0],ec='k',fc='r')
#    Ax.add_patch(Patch)

    # setup meridian and parallel printing
    if (np.ceil(xmax)-np.floor(xmin) <= 2): dx = 0.1
    elif (np.ceil(xmax)-np.floor(xmin) <= 5): dx = 0.2
    elif (np.ceil(xmax)-np.floor(xmin) <= 10): dx = 0.5
    else: dx = 1.

    if (np.ceil(ymax)-np.floor(ymin) <= 2): dy = 0.1
    elif (np.ceil(ymax)-np.floor(ymin) <= 5): dy = 0.2
    elif (np.ceil(ymax)-np.floor(ymin) <= 10): dy = 0.5
    else: dy = 1.

    Meridians = np.arange(np.floor(xmin),np.ceil(xmax),dx)
    Parallels = np.arange(np.floor(ymin),np.ceil(ymax),dy)

    Offset = np.abs(Map.xmax-Map.xmin)*0.05
    Meridians = Map.drawmeridians(Meridians,color='r',linewidth=1.5,labels=[1,0,0,1],fontsize=10,xoffset=Offset,yoffset=Offset,rotation=-Rotation)
    Parallels = Map.drawparallels(Parallels,color='r',linewidth=1.5,labels=[0,1,1,0],fontsize=10,xoffset=Offset,yoffset=Offset,rotation=Rotation)

#    for m in Meridians:
#        x,y = Meridians[m][0][0].get_xydata()[100]
#        plt.text(x,y,"$"+str(m)+"\/^{\circ}}$",rotation=-Rotation)
#    x1,y1 = Map(xmin-londiff,ymin-latdiff)
#    for m in Meridians:
#        x,y = Meridians[m][0][0].get_xydata()[0]
#
#        for i in range(0,len(Meridians[m][1])):
#            print Meridians[m][1][i].get_text()
#            print x, y
#            Meridians[m][1][i].set_position((-x,y))
#            Meridians[m][1][i].set_rotation(-Rotation)

    for p in Parallels:
        try:
            Parallels[p][1][0].set_rotation(Rotation)
        except:
            pass

    return Fig, Ax, Map

def ResampleRaster(InputRasterFile,OutputRasterFile,XResolution,YResolution=None,Format="ENVI"):

    """
    Description goes here...

    MDH

    """

    # import modules
    import rasterio, affine
    from rasterio.warp import reproject, Resampling

    # read the source raster
    with rasterio.open(InputRasterFile) as src:
        Array = src.read()
        OldResolution = src.res

        #setup output resolution
        if YResolution == None:
            YResolution = XResolution
        NewResolution = (XResolution,YResolution)


        # setup the transform to change the resolution
        XResRatio = OldResolution[0]/NewResolution[0]
        YResRatio = OldResolution[1]/NewResolution[1]
        NewArray = np.empty(shape=(Array.shape[0], int(round(Array.shape[1] * XResRatio)), int(round(Array.shape[2] * YResRatio))))
        Aff = src.affine
        NewAff = affine.Affine(Aff.a/XResRatio, Aff.b, Aff.c, Aff.d, Aff.e/YResRatio, Aff.f)

        # reproject the raster
        reproject(Array, NewArray, src_transform=Aff, dst_transform=NewAff, src_crs = src.crs, dst_crs = src.crs, resample=Resampling.bilinear)

        # write results to file
        with rasterio.open(OutputRasterFile, 'w', driver=src.driver, \
                            height=NewArray.shape[1],width=NewArray.shape[2], \
                            nodata=src.nodata,dtype=str(NewArray.dtype), \
                            count=src.count,crs=src.crs,transform=NewAff) as dst:
            dst.write(NewArray)

def ConvertRaster2LatLong(InputRasterFile,OutputRasterFile):

    """
    Convert a raster to lat long WGS1984 EPSG:4326 coordinates for global plotting

    MDH

    """

    # import modules
    import rasterio
    from rasterio.warp import reproject, calculate_default_transform as cdt, Resampling

    # read the source raster
    with rasterio.open(InputRasterFile) as src:
        #get input coordinate system
        Input_CRS = src.crs
        # define the output coordinate system
        Output_CRS = {'init': "epsg:4326"}
        # set up the transform
        Affine, Width, Height = cdt(Input_CRS,Output_CRS,src.width,src.height,*src.bounds)
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': Output_CRS,
            'transform': Affine,
            'affine': Affine,
            'width': Width,
            'height': Height
        })

        with rasterio.open(OutputRasterFile, 'w', **kwargs) as dst:
            for i in range(1, src.count+1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.affine,
                    src_crs=src.crs,
                    dst_transform=Affine,
                    dst_crs=Output_CRS,
                    resampling=Resampling.bilinear)

def PlotRaster(RasterFile, Map, alpha=1.):

    """
    Description goes here...

    MDH
    """

    print "Plotting raster..."

    #Read data
    gdata = gdal.Open(RasterFile, gdal.GA_ReadOnly)
    geo = gdata.GetGeoTransform()
    Data = gdata.ReadAsArray()

    # make topodat a masked array, masking values lower than sea level.
    Data = ma.masked_where(Data < 0, Data)

    #setup meshgrid for raster plotting
    xres = geo[1]
    yres = geo[5]
    xmin = geo[0] + xres * 0.5
    xmax = geo[0] + (xres * gdata.RasterXSize) - xres * 0.5
    ymin = geo[3] + (yres * gdata.RasterYSize) + yres * 0.5
    ymax = geo[3] - yres * 0.5

    x,y = np.mgrid[xmin:xmax+xres:xres, ymax+yres:ymin:yres]
    x,y = Map(x,y)
    Map.pcolormesh(x, y, Data.T, cmap=plt.cm.Greys, alpha=alpha)

def ReadShapeFile(ShapeFile):

    """
    Open shapefile and create Polygon dictionary
    returns dictionary of shapely Polygons

    MDH

    """

    #open shapefile and read shapes
    Shapes = fiona.open(ShapeFile)

    # get the input coordinate system
    Input_CRS = Proj(Shapes.crs)

    # Create a dictionary of shapely polygons
    PolygonDict = {}

    # loop through shapes and add to dictionary
    for Feature in Shapes:
        if Feature['geometry']['type'] == 'MultiPolygon':
            Shape = MultiPolygon(shape(Feature['geometry']))
            Value = float(Feature['properties']['ID'])
        elif Feature['geometry']['type'] == 'Polygon':
            Shape = Polygon(shape(Feature['geometry']))
            Value = float(Feature['properties']['ID'])

        #check for multipolygons
        if Value in PolygonDict:
            Polygons = [Shape, PolygonDict[Value]]
            Shape = cascaded_union(Polygons)
        #update dictionary
        PolygonDict[Value] = Shape

    return PolygonDict, Input_CRS

def ConvertPolygons2LatLong(PolygonDict, CRS):

    """
    Convert coordinates to latlong using pyproj
    https://gis.stackexchange.com/questions/127427/transforming-shapely-polygon-and-multipolygon-objects

    returns shapely Polygon/Multipolygon Dictionary

    MDH

    """

    # define the output coordinate system
    Output_CRS = Proj({'init': "epsg:4326"})

    NewPolygonDict = {}

    #loop through polygon dict and convert to latlong
    for Key, Poly in PolygonDict.iteritems():

        if Poly.geom_type == 'Polygon':
            # get x and y as arrays
            x,y = Poly.exterior.coords.xy
            lon,lat = transform(CRS, Output_CRS, x, y)
            Shape = Polygon(zip(lon,lat))

            #check for multipolygons
            if Key in NewPolygonDict:
                Polygons = [Shape, NewPolygonDict[Key]]
                Shape = cascaded_union(Polygons)

            #add converted polygon to new dictionary
            NewPolygonDict[Key] = Shape

        elif Poly.geom_type == 'MultiPolygon':
            for SinglePoly in Poly:
                # get x and y as arrays
                x,y = SinglePoly.exterior.coords.xy
                lon,lat = transform(CRS, Output_CRS, x, y)
                Shape = Polygon(zip(lon,lat))

                #check for multipolygons
                if Key in NewPolygonDict:
                    Polygons = [Shape, NewPolygonDict[Key]]
                    Shape = cascaded_union(Polygons)

                #add converted polygon to new dictionary
                NewPolygonDict[Key] = Shape

    return NewPolygonDict

def ConvertPolygonsLatLong2MapCoords(PolygonDict, Map):

    """
    Convert latlong coordinates to map coordinates

    returns shapely Polygon/Multipolygon Dictionary in map coordinates

    MDH

    """
    NewPolygonDict = {}

    #loop through polygon dict and convert to latlong
    for Key, Poly in PolygonDict.iteritems():

        if Poly.geom_type == 'Polygon':
            # get x and y as arrays
            lon,lat = Poly.exterior.coords.xy
            x,y = Map(lon,lat)
            Shape = Polygon(zip(x,y))

            #check for multipolygons
            if Key in NewPolygonDict:
                Polygons = [Shape, NewPolygonDict[Key]]
                Shape = cascaded_union(Polygons)

            #add converted polygon to new dictionary
            NewPolygonDict[Key] = Shape

        elif Poly.geom_type == 'MultiPolygon':
            for SinglePoly in Poly:
                # get x and y as arrays
                lon,lat = SinglePoly.exterior.coords.xy
                x,y = Map(lon,lat)
                Shape = Polygon(zip(x,y))

                #check for multipolygons
                if Key in NewPolygonDict:
                    Polygons = [Shape, NewPolygonDict[Key]]
                    Shape = cascaded_union(Polygons)

                #add converted polygon to new dictionary
                NewPolygonDict[Key] = Shape

    return NewPolygonDict

def ConvertShapefile2LatLong(ShapeFile):
    PolygonDict, CRS = ReadShapeFile(ShapeFile)
    PolygonDict = ConvertPolygons2LatLong(PolygonDict, CRS)
    return PolygonDict

def PlotPolygons(Polygons, Map=None, Ax=None, OutlineColour='k', FillColour='w', ColourMap="None", alpha=0.5):

    """
    Function to plot polygons from a shapely Polygon Dictionary

    Modified from PlottingRaster.py code by FJC

    Outline colour can be name, tuple or range of value to shade

    MDH

    """

    #create a figure if one doesnt already exist?
    if Ax == None:
        print("PlotPolygons: Warning, no axes provided, creating new figure and axes")
        Fig = plt.figure()
        Ax = plt.gca()
        plt.axis('equal')
        plt.xlabel('Longitude ($^o$)')
        plt.ylabel('Latitude ($^o$)')

    # convert to map coordinates
    if Map != None:
       Polygons = ConvertPolygonsLatLong2MapCoords(Polygons, Map)

    # loop through shapes in polygons and plot patches
    for Key, Poly in Polygons.iteritems():
        if Poly.geom_type == 'Polygon':
            Patch = PolygonPatch(Poly,fc=FillColour,ec=OutlineColour,alpha=alpha)
            Ax.add_patch(Patch)
        elif Poly.geom_type == 'MultiPolygon':
            for singlepoly in Poly:
                Patch = PolygonPatch(singlepoly,fc=FillColour,ec=OutlineColour,alpha=alpha)
                Ax.add_patch(Patch)

    if Ax == None:
        Ax.autoscale_view()


def PlotShapefile(ShapeFile, Map=None, Ax=None, OutlineColour='k', FillColour='w', ColourMap="None", alpha=0.5):
    print "Plotting shapefile..."
    Polygons = ConvertShapefile2LatLong(ShapeFile)
    PlotPolygons(Polygons, Map=Map, Ax=Ax, OutlineColour=OutlineColour, FillColour=FillColour, ColourMap=ColourMap, alpha=alpha)

def PlotHillshadeBasins():

    #Filename for hillshade raster
    BasinsFile = Directory+FilenamePrefix+"_AllBasins.shp"

    # setup the figure
    Fig, Ax, Map = CreateMapFigure(BasinsFile, AspectRatio=0, FigSizeFormat="default")

if __name__ == "__main__":
    #do something
    Directory = "/home/martin/bolinas/"
    DataDirectory = "/home/martin/bolinas/data/"
    PlotDirectory = "/home/martin/bolinas/plots/"
    RasterExtension = "bil"
    FilenamePrefix = "bolinas"

    # Load Basins and convert to shapefile, returning a shapely dict
    BasinsFile = FilenamePrefix+"_AllBasins.bil"
    ShapeFile = Directory+FilenamePrefix+"_AllBasins.shp"
    #Polygons = LSDMap_IO.PolygoniseRasterMerge(Directory, BasinsFile, ShapeFile)

    #ResampleRaster(Directory+"bolinas_hs.bil",Directory+"bolinas_hs_resample.bil",10)
    #ConvertRaster2LatLong(Directory+"bolinas_hs_resample.bil",Directory+"bolinas_hs_resample_latlong.bil")
    Fig, Ax, Map = CreateMapFigure(Directory+"bolinas_AllBasins.shp",FigSizeFormat="JGR",Rotation=45)
    PlotRaster(Directory+"bolinas_hs_resample_latlong.bil",Map)
    PlotShapefile(ShapeFile,Map,Ax,'k','w')
    plt.savefig(Directory+"example_map.png")
    print "Done."
