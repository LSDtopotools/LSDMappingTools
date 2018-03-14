#from __future__ import print_function
import sys
import os
from LSDPlottingTools import LSDMap_MOverNPlotting as MN
from LSDPlottingTools import LSDMap_SAPlotting as SA
from LSDMapFigure import PlottingHelpers as Helper

#let's test the basin appenders
DataDirectory = '/home/mhurst/movern_analysis/Sierras/'

# Helper.AppendBasinCSVs(DataDirectory)
Helper.AppendChiResidualsCSVs(DataDirectory)

