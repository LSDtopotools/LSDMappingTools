#=============================================================================
# These functions create figures for visualising the data from the drainage
# capture metrics
#
# It creates separate plots for each basin in the DEM.
#
# Authors:
#     Simon M. Mudd
#     Fiona J. Clubb
#=============================================================================
#=============================================================================
# IMPORT MODULES
#=============================================================================
import numpy as np
import LSDPlottingTools as LSDP
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.ticker as ticker
import pandas as pd
import math

#=============================================================================
# ANALYSIS FUNCTIONS
#=============================================================================

def ClockwiseAngleAndDistance(point):
    """
    Function to calculate clockwise angle and distance
    between two points.
    From https://stackoverflow.com/questions/41855695/sorting-list-of-two-dimensional-coordinates-by-clockwise-angle-using-python

    Author: FJC
    """
    # Vector between point and the origin: v = p - o
    vector = [point[0]-origin[0], point[1]-origin[1]]
    # Length of vector: ||v||
    lenvector = math.hypot(vector[0], vector[1])
    # If length is zero there is no angle
    if lenvector == 0:
        return -math.pi, 0
    # Normalize vector: v/||v||
    normalized = [vector[0]/lenvector, vector[1]/lenvector]
    dotprod  = normalized[0]*refvec[0] + normalized[1]*refvec[1]     # x1*x2 + y1*y2
    diffprod = refvec[1]*normalized[0] - refvec[0]*normalized[1]     # x1*y2 - y1*x2
    angle = math.atan2(diffprod, dotprod)
    # Negative angles represent counter-clockwise angles so we need to subtract them
    # from 2*pi (360 degrees)
    if angle < 0:
        return 2*math.pi+angle, lenvector
    # I return first the angle because that's the primary sorting criterium
    # but if two vectors have the same angle then the shorter distance should come first.
    return angle, lenvector

def SortBasinPerimeter(PerimeterDF, OutletNode):
    """
    This function takes the dataframe of perimeter nodes
    and sorts into a clockwise order based on the angle
    and distance between each node and the outlet node.

    Args:
        PerimeterDF: pandas dataframe of the perimeter nodes
        OutletNode (int): specify what the outlet node is.

    Author: FJC
    """
    df
    #get x and y coords of the perimeter
    perimeter_x = PerimeterDF.x
    perimeter_y = PerimeterDF.y

    df.assign(zip(perimeter_x, perimeter_y))
    print df

    OutletRow = df.loc[df['node'] == OutletNode]
    origin = [OutletRow.x, OutletRow.y]

    # get the coordinates of the outlet junction
    refvec = [0, 1]
    df.assign()
    sorted_points = sorted(points, key=ClockwiseAngleAndDistance)

    sorted_points = zip(*sorted_points)
    sorted_x = sorted_points[0]
    sorted_y = sorted_points[1]

    # now sort the dataframe based on this
