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

def location_map():

    ax = plt.axes(projection=ccrs.PlateCarree())
    
    boscastle_lon, boscastle_lat = -4.692, 50.686
    ryedale_lon, ryedale_lat = -1.056, 54.246
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

    ax.add_feature(land_50m, edgecolor='gray')
    
    plt.text(boscastle_lon , boscastle_lat, 'Boscastle',
     horizontalalignment='right',
     transform=ccrs.PlateCarree())
    
    plt.text(ryedale_lon, ryedale_lat, 'Helmsley',
     horizontalalignment='right',
     transform=ccrs.PlateCarree())
    
    # Bug to fix. Markers do not appear when add_feature used above??
    plt.scatter(boscastle_lon, boscastle_lat,
         color='blue', marker='o',
         transform=ccrs.PlateCarree(),
         )
    
    plt.scatter(ryedale_lon, ryedale_lat,
         color='blue', marker='o',
         transform=ccrs.PlateCarree(),
         )

    # Set the longitude and latitude extents [West, East, South, North]
    ax.set_extent([-7, 2, 50, 59])
    
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


if __name__ == '__main__':
    location_map()