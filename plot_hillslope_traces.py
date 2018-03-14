#import modules
from __future__ import print_function
from geopandas import GeoDataFrame
import pandas as pd
import numpy as np
from sys import platform, stdout

# import plotting tools and set the back end for running on server
import matplotlib
matplotlib.use('Agg')

from matplotlib import rcParams, ticker, gridspec, cm
import matplotlib.pyplot as plt

# import mapping tools
import rotated_mapping_tools as rmt
import LSDPlottingTools as LSDP
from LSDMapFigure import PlottingHelpers as Helper
from LSDMapFigure.PlottingRaster import MapFigure
from LSDMapFigure.PlottingRaster import BaseRaster

# set figure sizes (in inches) based on format
FigSizeFormat = "EPSL"
if FigSizeFormat == "geomorphology":
    FigWidth_Inches = 6.25
elif FigSizeFormat == "big":
    FigWidth_Inches = 16
elif FigSizeFormat == "small":
    FigWidth_Inches = 3.3
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
 
# Set up fonts for plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 8
rcParams['text.usetex'] = True
   
# define filenames and relative workspace
fname_prefix = "bolinas"
HillshadeName = "bolinas_hs.bil"
Directory = "/home/mhurst/bolinas_paper/"
DataDirectory = Directory+"data/"
ChannelDataDirectory = DataDirectory+"channel_data/"
HillslopeDataDirectory = DataDirectory+"hillslope_data/"
HilltopPointsData = HillslopeDataDirectory+"bolinas_HilltopData.csv"
ChannelPointsDat = ChannelDataDirectory+"bolinas_MChiSegmented.csv"
ChannelHeadPointsData = ChannelDataDirectory+"bolinas_CH_wiener_nodeindices_for_Arc.csv"

# create the map figure
MF = MapFigure(HillshadeName, DataDirectory, coord_type="UTM_km", colourbar_location='None')

# add hilltops
HilltopPointsDF = pd.read_csv(HilltopPointsData)
HilltopPoints = LSDP.LSDMap_PointData(HilltopPointsDF, data_type = "pandas", PANDEX = True)
MF.add_point_data(HilltopPoints,alpha=0.5,zorder=100,unicolor="blue",manual_size=5)

# add channel heads
#ChannelHeadsDF = pd.read_csv(ChannelHeadPointsData)
#ChannelHeadPoints = LSDP.LSDMap_PointData(ChannelHeadsDF, data_type = "pandas", PANDEX = True)
#MF.add_point_data(ChannelHeadPoints,alpha=0.5,zorder=100,unicolor="blue",manual_size=5)

# add channels
#ChannelDF = Helper.ReadChiDataMapCSV(ChannelDataDirectory,fname_prefix)
#ChannelPoints = LSDP.LSDMap_PointData(ChannelDF, data_type = "pandas", PANDEX = True)
#MF.add_point_data(ChannelPoints,show_colourbar="False", scale_points=True, column_for_scaling='drainage_area',alpha=0.5,zorder=90)

# add hillslope traces    
#Plot HillslopeTraces():

# Save the figure
ImageName = Directory+"plots/bolinas_traces.png"
MF.save_fig(fig_width_inches = FigWidth_Inches, FigFileName = ImageName, FigFormat="png", Fig_dpi = 300)
