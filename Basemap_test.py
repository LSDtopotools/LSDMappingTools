import matplotlib
matplotlib.use('Agg')


from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import matplotlib.cm
 
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize


def FullGlobePerspective():

    fig, ax = plt.subplots(figsize=(12, 10))

    # setup Lambert Conformal basemap.
    m = Basemap(width=1800000,height=1500000,projection='nsper',area_thresh = 100000,
                resolution='h',lat_1=45.,lat_2=55,lat_0=25.7,lon_0=91.5, satellite_height = 10000000)
    # draw coastlines.
    m.drawcoastlines()
    # draw a boundary around the map, fill the background.
    # this background will end up being the ocean color, since
    # the continents will be drawn on top.
    m.drawmapboundary(fill_color='aqua')
    # fill continents, set lake color same as ocean color.
    m.fillcontinents(color='peru',lake_color='aqua')
    # draw parallels and meridians.
    # label parallels on right and top
    # meridians on bottom and left
    parallels = np.arange(0.,81,5.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels,labels=[False,True,True,False])
    meridians = np.arange(10.,351.,5.)
    m.drawmeridians(meridians,labels=[True,False,False,True])
    m.drawcountries()
    m.readshapefile('Mega_divide_footprint', 'MDF')

    # All this stuff from:
    # http://www.datadependence.com/2016/06/creating-map-visualisations-in-python/
    df_poly = pd.DataFrame({
            'shapes': [Polygon(np.array(shape), True) for shape in m.MDF]})

    print(df_poly)

    #df_poly = df_poly.merge(new_areas, on='area', how='left')
    cmap = plt.get_cmap('Oranges')   
    pc = PatchCollection(df_poly.shapes, zorder=2)
    pc.set_facecolor(cmap(1))
    ax.add_collection(pc)  

    FigFileName = "Test.png"
    FigFormat = "png"

    plt.savefig(FigFileName,format=FigFormat,dpi=200)
    
def ZoomPerspective():
    
    fig, ax = plt.subplots(figsize=(6, 5))

    # setup Lambert Conformal basemap.
    print("I am trying to get rid of these bloody lakes")
    m = Basemap(width=1800000,height=1500000,projection='lcc',area_thresh=10000000.,
                resolution='l',lat_1=45.,lat_2=55,lat_0=25.7,lon_0=91.5)
    # draw coastlines.
    m.drawcoastlines()
    # draw a boundary around the map, fill the background.
    # this background will end up being the ocean color, since
    # the continents will be drawn on top.
    m.drawmapboundary(fill_color='aqua')
    # fill continents, set lake color same as ocean color.
    m.fillcontinents(color='peru',lake_color='aqua')
    # draw parallels and meridians.
    # label parallels on right and top
    # meridians on bottom and left
    parallels = np.arange(0.,81,5.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels,labels=[False,True,True,False])
    meridians = np.arange(10.,351.,5.)
    m.drawmeridians(meridians,labels=[True,False,False,True])
    m.drawcountries()
    m.readshapefile('Mega_divide_footprint', 'MDF')

    # All this stuff from:
    # http://www.datadependence.com/2016/06/creating-map-visualisations-in-python/
    df_poly = pd.DataFrame({
            'shapes': [Polygon(np.array(shape), True) for shape in m.MDF]})

    print(df_poly)

    #df_poly = df_poly.merge(new_areas, on='area', how='left')
    cmap = plt.get_cmap('Oranges')   
    pc = PatchCollection(df_poly.shapes, zorder=2, alpha = 0.5)
    pc.set_facecolor("r")
    ax.add_collection(pc)  

    FigFileName = "Test.png"
    FigFormat = "png"

    plt.savefig(FigFileName,format=FigFormat,dpi=200)   

    
def ZoomPerspectiveEtopo():
    
    fig, ax = plt.subplots(figsize=(3, 2.5))

    # setup Lambert Conformal basemap.
    print("I am trying to get rid of these bloody lakes")
    m = Basemap(width=1800000,height=1500000,projection='lcc',area_thresh=10000000.,
                resolution='h',lat_1=45.,lat_2=55,lat_0=25.7,lon_0=91.5)
    # draw coastlines.
    #m.drawrivers()
    m.drawcoastlines()
    # draw a boundary around the map, fill the background.
    # this background will end up being the ocean color, since
    # the continents will be drawn on top.
    m.drawmapboundary(fill_color='aqua')
    # fill continents, set lake color same as ocean color.
    m.fillcontinents(color='peru',lake_color='white')
    # draw parallels and meridians.
    # label parallels on right and top
    # meridians on bottom and left
    parallels = np.arange(0.,81,5.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels,labels=[False,True,True,False])
    meridians = np.arange(10.,351.,5.)
    m.drawmeridians(meridians,labels=[True,False,False,True])
    m.drawcountries()
    m.readshapefile('Mega_divide_footprint', 'MDF')

    # All this stuff from:
    # http://www.datadependence.com/2016/06/creating-map-visualisations-in-python/
    df_poly = pd.DataFrame({
            'shapes': [Polygon(np.array(shape), True) for shape in m.MDF]})

    print(df_poly)

    #df_poly = df_poly.merge(new_areas, on='area', how='left')
    cmap = plt.get_cmap('Oranges')   
    pc = PatchCollection(df_poly.shapes, zorder=2, alpha = 0.5)
    pc.set_facecolor("gray")
    ax.add_collection(pc)  

    FigFileName = "Test.svg"
    FigFormat = "svg"

    plt.savefig(FigFileName,format=FigFormat,dpi=200)       
    
    
if __name__ == "__main__":
     ZoomPerspectiveEtopo()
