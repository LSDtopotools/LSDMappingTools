# -*- coding: utf-8 -*-
"""
Location map test script

Plots a location map using the Cartopy package

Install cartopy first for this to work.

http://scitools.org.uk/cartopy/docs/v0.13/index.html

Add annotations for locations using their lon/lats.

Author: DAV

"""

import LSDPlottingTools.locationmap as locmap

extent = [-7.5, 2, 50, 59]
# Placename1 : (LATITUDE, LONGITUDE), Placename2...
gazetter = { 'Boscastle' : (50.686, -4.692),
            'Helmsley'  : (54.246, -1.056)
           }

locmap.location_map(extent, gazetter, offset = 0.1)