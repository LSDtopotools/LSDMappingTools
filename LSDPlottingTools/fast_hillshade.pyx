#fast_hillshade.pyx
# This is cython version of the nicer looking hillshade function in the
# LSDTopoTools core C++ libraries.

# Cython rule of thumb no 1. If there are equivalent C-libraries for
# numpy stuff, use them. (E.g. math functions)
# Let's use the native C-libraries for math functions.
from libc.math cimport sin, cos, sqrt, M_PI, atan, atan2

import cython
cimport cython

# Author: DAV, 2017
import numpy as np
cimport numpy as np

# Fix a data type for our arrays.
DTYPE = np.float64
# Define a compile type to DTYPE_t
ctypedef np.float64_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def Hillshade(np.ndarray[DTYPE_t, ndim=2] terrain_array,
              float DataResolution, float azimuth = 315,
              float angle_altitude = 45,
              float NoDataValue = -9999, float z_factor = 1):
  """Creates a hillshade raster

  Args:
      raster_file (str): The name of the raster file with path and extension
        (Can also be numpy array of loaded data.)
      azimuth (float): Azimuth of sunlight
      angle_altitude (float): Angle altitude of sun
      NoDataValue (float): The nodata value of the raster

  Returns:
      HSArray (numpy.array): The hillshade array

  Author:
      DAV, SWDG, SMM

  """
  #assert terrain_array.dtype == DTYPE_t

  cdef unsigned int ncols = terrain_array.shape[0]
  cdef unsigned int nrows = terrain_array.shape[1]

  # DAV attempting mask nodata vals
  #nodata_mask = terrain_array == NoDataValue
  #terrain_array[nodata_mask] = np.nan

  # Ndarray best choice? Will revisit later...
  cdef np.ndarray[DTYPE_t, ndim=2] HSarray = np.empty((ncols,nrows))
  HSarray.fill(np.nan)

  #cdef float M_PI = np.pi

  cdef float zenith_rad = (90 - angle_altitude) * M_PI / 180.0
  cdef float azimuth_math = 360-azimuth + 90
  if (azimuth_math >= 360.0):
    azimuth_math = azimuth_math - 360
  cdef float azimuth_rad = azimuth_math * M_PI  / 180.0

  cdef float slope_rad = 0
  cdef float aspect_rad = 0
  cdef float dzdx = 0
  cdef float dzdy = 0
  cdef int i, j = 0
  # For loop necessary? Revisit...

  for i in range(0, ncols - 1):
    for j in range(0, nrows - 1):
        if terrain_array[i, j] != NoDataValue:
          # No need to check nodata value every iter, they are nans handled by numpy?
          dzdx = (((terrain_array[i, j+1] + 2*terrain_array[i+1, j] + terrain_array[i+1, j+1]) -
                  (terrain_array[i-1, j-1] + 2*terrain_array[i-1, j] + terrain_array[i-1, j+1]))
                  / (8 * DataResolution))
          dzdy = (((terrain_array[i-1, j+1] + 2*terrain_array[i, j+1] + terrain_array[i+1, j+1]) -
                  (terrain_array[i-1, j-1] + 2*terrain_array[i, j-1] + terrain_array[i+1, j-1]))
                  / (8 * DataResolution))

          slope_rad = atan(z_factor * sqrt((dzdx * dzdx) + (dzdy * dzdy)))

          if (dzdx != 0):
            aspect_rad = atan2(dzdy, (dzdx*-1))
            if (aspect_rad < 0):
              aspect_rad = 2 * M_PI + aspect_rad

          else:
            if (dzdy > 0):
              aspect_rad = M_PI / 2
            elif (dzdy < 0):
              aspect_rad = 2 * M_PI - M_PI / 2
            else:
              aspect_rad = aspect_rad
          # Same comment as above...for loop assignment? Use array instead of ndarray
          # and then whole  array operation
          HSarray[i, j] = 255.0 * ((cos(zenith_rad) * cos(slope_rad)) +
                          (sin(zenith_rad) * sin(slope_rad) *
                          cos(azimuth_rad - aspect_rad)))
          # Necessary?
          if (HSarray[i, j] < 0):
            HSarray[i, j] = 0

  return HSarray
