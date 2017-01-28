# -*- coding: utf-8 -*-
"""
Location map.

Plots a location map using the Cartopy package

Install cartopy first for this to work.

http://scitools.org.uk/cartopy/docs/v0.13/index.html

Add annotations for locations using their lon/lats.

Author: DAV

"""

import cartopy.feature as cfeature
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.patches as patches


from shapely.geometry.polygon import LinearRing
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

def location_map(extent, gazetter, offset=0.):
    """Plots a series of points marking towns/sample sites/locations etc 
       given a dictionary of places and lat/lons.
       
    Arguments:
        extent (list): A list of the cooridinates of the bounding extent of the
               location map in format: [West, East, South, North] e.g.:
                   [lonW, lonE, latS, latN]
                   
        gazetter (dict): A dictionary of arbitrary length of the format:
            { 'Placename1' : (LATITUDE, LONGITUDE),
              'Placename2' : (LATITUDE, LONGITUDE).
              and so on...}
            
        offset (float): Offset of the text label to the marker points, 
               in degrees.
              
    Todo: Greying out the landmass appears to block out the marker points
          (weird...) so this has been commented out and you just get an 
          outline map for now.
               
    Author: DAV
        
    """
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Put a background image on for nice sea rendering.

    ax.coastlines('50m')
    # Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
    # DAV - 50m scale not enough for UK
    # Using '10m' works but too many boundaries as it defaults to unitary authorities
    # in England and Council areas in Wales/Scotland so appears way too cluttered.
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')
    
    # Add grey for the land
    land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])
    
    ax.add_feature(states_provinces, edgecolor='gray')
    #ax.xaxis.set_visible(True)
    #ax.yaxis.set_visible(True)
    
    ax.set_yticks([50,54,58], crs=ccrs.PlateCarree())
    ax.set_xticks([-6, -4, -2, 0], crs=ccrs.PlateCarree())
    
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    
    #ax.add_feature(land_50m, edgecolor='gray')
    
    for location in gazetter:
        plt.text(gazetter[location][1] - offset, gazetter[location][0] + offset, location,
                 horizontalalignment='right',
                 transform=ccrs.PlateCarree())
        
        plt.scatter(gazetter[location][1], gazetter[location][0],
                    color='black', marker='o',
                    transform=ccrs.PlateCarree())
    
    """
    plt.text(boscastle_lon - 0.1, boscastle_lat + 0.1, 'Boscastle',
     horizontalalignment='right',
     transform=ccrs.PlateCarree())
    
    plt.text(ryedale_lon + 0.3, ryedale_lat + 0.3, 'Helmsley',
     horizontalalignment='left',
     transform=ccrs.PlateCarree())
    
    # Bug to fix. Markers do not appear when add_feature used above??
    plt.scatter(boscastle_lon, boscastle_lat,
         color='black', marker='o',
         transform=ccrs.PlateCarree(),
         )
    
    plt.scatter(ryedale_lon, ryedale_lat,
         color='black', marker='o',
         transform=ccrs.PlateCarree(),
         )
    """
    # Set the longitude and latitude extents [West, East, South, North]
    ax.set_extent(extent)
    
    # Draw a bounding box square
    """
    lons = [-5.8, -5.8, -5.5, -5.5]
    lats = [50.27, 50.48, 50.48, 50.27]
    ring = LinearRing(list(zip(lons, lats)))
    ax.add_geometries([ring], ccrs.PlateCarree(), facecolor='none', edgecolor='blue')
    """

    """This method didn't work
    # Add location rectangle
    rect = patches.Rectangle((-4,1),2,2,linewidth=1,edgecolor='r',facecolor='none')
    ax.add_patch(patches.Rectangle(xy=[1, 1], width=2, height=2,
                                    facecolor='blue',
                                    alpha=0.2,
                                    transform=ccrs.Robinson()))
    """

    plt.show()
