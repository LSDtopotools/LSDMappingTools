#=============================================================================
# These functions create figures for visualising the hillslope-channel data
# produced by the hillslope morphology driver
#
# It creates separate plots for each basin in the DEM.
#
# Authors:
#     Martin D. Hurst
#     Fiona J. Clubb
#=============================================================================
#=============================================================================
# IMPORT MODULES
#=============================================================================
import geopandas as gpd
from geopandas import GeoDataFrame
import pandas as pd
import numpy as np
from shapely.geometry import LineString, shape, Point, MultiPolygon, Polygon
import numpy.ma as ma
from scipy import stats
import os.path, sys
import math

# import the basemap library
from mpl_toolkits.basemap import Basemap
from osgeo import gdal
import fiona
from pyproj import Proj, transform

# import plotting tools and set the back end for running on server
import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams, ticker
import matplotlib.pyplot as plt
from matplotlib import cm
from shapely.ops import cascaded_union
import matplotlib.colors as colors
from descartes import PolygonPatch

# plotting tools for using LSDMapFigure
import LSDPlottingTools as LSDP
from LSDMapFigure import PlottingHelpers as Helper
from LSDMapFigure.PlottingRaster import MapFigure
from LSDMapFigure.PlottingRaster import BaseRaster

#=============================================================================
# PRELIMINARY FUNCTIONS
#=============================================================================

def Get_FigWidth_Inches(FigSizeFormat="default"):

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

    return FigWidth_Inches

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

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = 8
    rcParams['text.usetex'] = True

    FigWidth_Inches = Get_FigWidth_Inches(FigSizeFormat)
    Fig = plt.figure(figsize=(FigWidth_Inches,FigWidth_Inches/AspectRatio))

    return Fig

def ReadHillslopeData(DataDirectory, FilenamePrefix):
    """
    This function reads in the file with the suffix '_HilltopData.csv'
    to a pandas dataframe

    Args:
        DataDirectory: the data directory
        FilenamePrefix: the file name prefix

    Returns:
        pandas dataframe with data from the csv file

    Author: MDH
    """
    # get the csv filename
    Suffix = '_HilltopData.csv'
    Filename = FilenamePrefix+Suffix

    # read in the dataframe using pandas
    HillslopeData = pd.read_csv(DataDirectory+Filename)

    # drop any rows with no data (hillslope traces to outside the study domain)
    # or with value of -9999 for Basin ID
    HillslopeData = HillslopeData.dropna()
    HillslopeData = HillslopeData[HillslopeData.BasinID != -9999]

    #return the hillslope data
    return HillslopeData

def ReadChannelData(DataDirectory, FilenamePrefix):
    """
    This function reads in the file with the suffix '_MChiSegmented.csv'
    or _MChiSegmented.geojson to a pandas dataframe

    Args:
        DataDirectory: the data directory
        FilenamePrefix: the file name prefix

    Returns:
        pandas dataframe with data from the csv file

    Author: MDH
    """
    # get the filename and open either csv or geojson
    Suffix = '_MChiSegmented'
    Filename = FilenamePrefix+Suffix

    if os.path.isfile(DataDirectory+Filename+".csv"):
        # read in the dataframe using pandas
        ChannelData = pd.read_csv(DataDirectory+Filename+".csv")

    elif os.path.isfile(DataDirectory+Filename+".geojson"):
        # read in the dataframe using pandas
        ChannelData = gpd.read_file(DataDirectory+Filename+".geojson")
    else:
        print("No file named "+DataDirectory+Filename+".* found")
        sys.exit()

    # If there is no chi values due to threshold then chi will be -9999
    # throw out these segments
    Segments2Remove = ChannelData[ChannelData.chi == -9999].segment_number.unique()
    ChannelData = ChannelData[~ChannelData.segment_number.isin(Segments2Remove)]

    #return the hillslope data
    return ChannelData

def ReadHillslopeTraces(DataDirectory, FilenamePrefix,ThinningFactor=1,CustomExtent=[-9999]):
    """
    This function reads in the file with the suffix '_hillslope_traces.csv'
    and creates a geopandas GeoDataFrame

    Args:
        DataDirectory: the data directory
        FilenamePrefix: the file name prefix
        ThinningFactor: An integer to skip every X traces for speed and clarity
        CustomExtent: A list containing [xmin, xmax, ymin, ymax] so as only to consider the traces in the plotting area, default is plot all traces

    Returns:
        geopandas GeoDataFrame with data from the csv file spatially organised

    Author: MDH
    """

    # get the csv filename
    Suffix = '_hillslope_traces'
    Extension = '.csv'
    ReadFilename = DataDirectory+FilenamePrefix+Suffix+Extension

    # read in the dataframe using pandas and convert to geopandas geodataframe
    df = pd.read_csv(ReadFilename)

    # thin the data
    df = df.iloc[::ThinningFactor,:]

    # clip to custom extent_raster
    if len(CustomExtent) == 4:
      df.drop(df[df.Easting < CustomExtent[0]].index, inplace=True)
      df.drop(df[df.Easting > CustomExtent[1]].index, inplace=True)
      df.drop(df[df.Northing < CustomExtent[2]].index, inplace=True)
      df.drop(df[df.Northing > CustomExtent[3]].index, inplace=True)

    # check for and delete any traces tat are only 1 point long since these wont plot
    df['is_unique'] = ~df['HilltopID'].duplicated(keep=False)
    temp = ~df['HilltopID'].duplicated(keep=False)
    df.drop(df[df.is_unique == True].index, inplace=True)

    geometry = [Point(xy) for xy in zip(df.Easting, df.Northing)]
    df = df.drop(['Easting','Northing','Longitude', 'Latitude'], axis=1)
    crs = {'init': 'epsg:4326'}
    geo_df = GeoDataFrame(df, crs=crs, geometry=geometry)

    return geo_df

def ReadTerraceData(DataDirectory,FilenamePrefix):
    """
    This function reads in the file with the suffix '_terrace_info.csv'
    and creates a geopandas GeoDataFrame

    Args:
        DataDirectory: the data directory
        FilenamePrefix: the file name prefix

    Returns:
        geopandas GeoDataFrame with data from the csv file spatially organised

    Author: FJC
    """
    # get the csv filename
    Suffix = '_terrace_info'
    Extension = '.csv'
    ReadFilename = DataDirectory+FilenamePrefix+Suffix+Extension

    # read in the dataframe using pandas and convert to geopandas geodataframe
    df = pd.read_csv(ReadFilename)
    geometry = [Point(xy) for xy in zip(df.Longitude, df.Latitude)]
    df = df.drop(['X','Y','Longitude', 'Latitude'], axis=1)
    crs = {'init': 'epsg:4326'}
    geo_df = GeoDataFrame(df, crs=crs, geometry=geometry)

    return geo_df

def ReadBaselineData(DataDirectory,fname_prefix):
    """
    This function reads in the csv file with the extension "_baseline_channel_info.csv"
    and returns it as a pandas dataframe. Used for comparing
    terrace locations to the main channel.

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the name of the DEM

    Returns:
        pandas dataframe with the channel info

    Author: FJC
    """
    csv_suffix = '_baseline_channel_info.csv'
    fname = DataDirectory+fname_prefix+csv_suffix

    df = pd.read_csv(fname)

    return df

def MapBasinKeysToJunctions(DataDirectory,FilenamePrefix):
    """
    Function to write a dict of basin keys vs junctions

    Author: FJC
    """
    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)
    #print BasinChannelData

    # load the hillslopes data
    HillslopeData = ReadHillslopeData(DataDirectory, FilenamePrefix)

    basin_keys = ChannelData.basin_key.unique()
    basin_junctions = HillslopeData.BasinID.unique()

    basin_dict = {}

    for i, key in enumerate(basin_keys):
        basin_dict[key] = basin_junctions[i]

    print basin_dict
    return basin_dict

def WriteHillslopeTracesShp(DataDirectory,FilenamePrefix,ThinningFactor=1, CustomExtent=[-9999]):
    """
    This function writes a shapefile of hillslope traces

    Args:
        DataDirectory: the data directory
        FilenamePrefix: the file name prefix


    Author: MDH
    """

    #read the raw data to geodataframe
    geo_df = ReadHillslopeTraces(DataDirectory,FilenamePrefix,ThinningFactor, CustomExtent)
    Suffix = '_hillslope_traces'
    WriteFilename = DataDirectory+FilenamePrefix+Suffix+'.shp'

    #aggregate these points with group by
    geo_df = geo_df.groupby(['HilltopID'])['geometry'].apply(lambda x: LineString(x.tolist()))
    geo_df = GeoDataFrame(geo_df, geometry='geometry')
    geo_df.to_file(WriteFilename, driver='ESRI Shapefile')

def SaveHillslopeDataByBasin(DataDirectory,FilenamePrefix):
    """
    This function organises hillslope data by basin number
    and writes the results to a new, numbered file

    Args:
        DataDirectory: the data directory
        FilenamePrefix: the file name prefix

    Returns:
        writes new files to data directory

    Author: MDH
    """

    # load the hillslope data
    HillslopeData = ReadHillslopeData(DataDirectory, FilenamePrefix)

    # get a list of basins
    Basins = HillslopeData.BasinID.unique()

    # get the csv filename
    Suffix = '_HilltopData.csv'

    # loop through basins
    for i in range(0,len(Basins)):
        #isolate basin data
        BasinHillslopeData = HillslopeData[HillslopeData.BasinID == Basins[i]]
        #setup an output file
        OutputFilename = DataDirectory + FilenamePrefix + "_" + str(i) + Suffix
        #write to file
        BasinHillslopeData.to_csv(OutputFilename, index=False)

def SaveChannelDataByBasin(DataDirectory,FilenamePrefix):
    """
    This function organises channel data by basin number
    and writes the results to a new, numbered file

    Args:
        DataDirectory: the data directory
        FilenamePrefix: the file name prefix

    Returns:
        writes new files to data directory

    Author: MDH
    """

    # load the hillslope data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)

    # get a list of basins
    Basins = ChannelData.basin_key.unique()

    # get the csv filename
    Suffix = '_MChiSegmented.csv'

    # loop through basins
    for i in range(0,len(Basins)):
        #isolate basin data
        BasinChannelData = ChannelData[ChannelData.basin_key == Basins[i]]
        #setup an output file
        OutputFilename = DataDirectory + FilenamePrefix + "_" + str(i) + Suffix
        #write to file
        BasinChannelData.to_csv(OutputFilename, index=False)

def PlotChiElevationSegments(DataDirectory, FilenamePrefix, PlotDirectory, BasinID):

    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)

    # isolate basin data
    BasinChannelData = ChannelData[ChannelData.basin_key == BasinID]
    MinimumChi = BasinChannelData.chi.min()

    # how many segments are we dealing with?
    Segments = BasinChannelData.segment_number.unique()

    # setup the figure
    CreateFigure()

    # Get the data columns for plotting
    for i in range(0, len(Segments)):
        #get data arrays
        Chi = ChannelData.chi[ChannelData.segment_number == Segments[i]]
        Elevation = ChannelData.elevation[ChannelData.segment_number == Segments[i]]
        SegmentedElevation = ChannelData.segmented_elevation[ChannelData.segment_number == Segments[i]]
        #normalise chi by outlet chi
        Chi = Chi-MinimumChi
        #plot, colouring segments
        Colour = np.random.rand()
        plt.plot(Chi,Elevation,'k--',dashes=(2,2), lw=0.5,zorder=10)
        plt.plot(Chi, SegmentedElevation, '-', lw=2, c=plt.cm.Paired(Colour),zorder=9)

    # Finalise the figure
    plt.xlabel(r'$\chi$ (m)')
    plt.ylabel('Elevation (m)')
    plt.title('Basin ID ' + str(BasinID))
    plt.tight_layout()
    plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(BasinID) + "_ChiElevSeg.png", dpi=300)

def PlotLongProfileSegments(DataDirectory, FilenamePrefix, PlotDirectory, BasinID):

    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)

    # isolate basin data
    BasinChannelData = ChannelData[ChannelData.basin_key == BasinID]
    MinimumDistance = BasinChannelData.flow_distance.min()

    # how many segments are we dealing with?
    Segments = BasinChannelData.segment_number.unique()

    # setup the figure
    CreateFigure()

    # Get the data columns for plotting
    for i in range(0, len(Segments)):
        #get data arrays
        Dist = ChannelData.flow_distance[ChannelData.segment_number == Segments[i]]
        Elevation = ChannelData.elevation[ChannelData.segment_number == Segments[i]]
        SegmentedElevation = ChannelData.segmented_elevation[ChannelData.segment_number == Segments[i]]
        #normalise distance by outlet distance
        Dist = Dist-MinimumDistance
        #plot, colouring segments
        Colour = np.random.rand()
        plt.plot(Dist/1000,Elevation,'k--',dashes=(2,2), lw=0.5,zorder=10)
        plt.plot(Dist/1000, SegmentedElevation, '-', lw=2, c=plt.cm.Paired(Colour),zorder=9)

    # Finalise the figure
    plt.xlabel('Distance (km)')
    plt.ylabel('Elevation (m)')
    plt.title('Basin ID ' + str(BasinID))
    plt.tight_layout()
    plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(BasinID) + "_LongProfSeg.png", dpi=300)

def PlotChiElevationMChi(DataDirectory, FilenamePrefix, PlotDirectory, BasinID):

    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)

    # isolate basin data
    BasinChannelData = ChannelData[ChannelData.basin_key == BasinID]
    MinimumChi = BasinChannelData.chi.min()
    MaximumMChi = BasinChannelData.m_chi.max()

    # how many segments are we dealing with?
    Segments = BasinChannelData.segment_number.unique()

    # separate into main stem and trib data
    MainStemChannelData = BasinChannelData[BasinChannelData.source_key == 0]
    MainStemSegments = MainStemChannelData.segment_number.unique()

    # setup the figure
    Fig = CreateFigure(AspectRatio=4./3.)

    #choose colormap
    ColourMap = cm.coolwarm

    # Get the data columns for plotting
    for i in range(0, len(Segments)):
        if Segments[i] in MainStemSegments:
            #get data arrays
            Chi = ChannelData.chi[ChannelData.segment_number == Segments[i]]
            Elevation = ChannelData.elevation[ChannelData.segment_number == Segments[i]]
            SegmentedElevation = ChannelData.segmented_elevation[ChannelData.segment_number == Segments[i]]
            MChi = ChannelData.m_chi[ChannelData.segment_number == Segments[i]].unique()[0]

            #normalise chi by outlet chi
            Chi = Chi-MinimumChi
            #plot, colouring segments
            Colour = MChi/MaximumMChi
            plt.plot(Chi,Elevation,'k--',dashes=(2,2), lw=0.5,zorder=10)
            plt.plot(Chi, SegmentedElevation, '-', lw=2, c=ColourMap(Colour),zorder=9)

    # Finalise the figure
    plt.xlabel(r'$\chi$ (m)')
    plt.ylabel('Elevation (m)')
    plt.title('Basin ID ' + str(BasinID))
    plt.tight_layout()
    #add colourbar
    CAx = Fig.add_axes([0.15,0.8,0.4,0.05])
    m = cm.ScalarMappable(cmap=ColourMap,norm=colors.Normalize(vmin=0, vmax=MaximumMChi))
    m.set_array(ChannelData.m_chi)
    plt.colorbar(m, cax=CAx,orientation='horizontal')
    CAx.tick_params(labelsize=6)
    plt.xlabel('$k_{sn}$')
    #save output
    plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(BasinID) + "_ChiElevMChi.png", dpi=300)

def PlotLongProfileMChi(DataDirectory, FilenamePrefix, PlotDirectory, BasinID):

    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)

    # isolate basin data
    BasinChannelData = ChannelData[ChannelData.basin_key == BasinID]
    MinimumDistance = BasinChannelData.flow_distance.min()
    MaximumMChi = BasinChannelData.m_chi.max()

    # how many segments are we dealing with?
    Segments = BasinChannelData.segment_number.unique()

    # setup the figure
    Fig = CreateFigure()

    #choose colormap
    ColourMap = cm.viridis

    # Get the data columns for plotting
    for i in range(0, len(Segments)):
        #get data arrays
        Dist = ChannelData.flow_distance[ChannelData.segment_number == Segments[i]]
        Elevation = ChannelData.elevation[ChannelData.segment_number == Segments[i]]
        SegmentedElevation = ChannelData.segmented_elevation[ChannelData.segment_number == Segments[i]]
        MChi = ChannelData.m_chi[ChannelData.segment_number == Segments[i]].unique()[0]

        #normalise distance by outlet distance
        Dist = Dist-MinimumDistance
        #plot, colouring segments
        Colour = MChi/MaximumMChi
        plt.plot(Dist/1000,Elevation,'k--',dashes=(2,2), lw=0.5,zorder=10)
        plt.plot(Dist/1000, SegmentedElevation, '-', lw=2, c=ColourMap(Colour),zorder=9)

    # Finalise the figure
    plt.xlabel('Distance (km)')
    plt.ylabel('Elevation (m)')
    plt.title('Basin ID ' + str(BasinID))
    plt.tight_layout()
    #add colourbar
    CAx = Fig.add_axes([0.15,0.8,0.4,0.05])
    m = cm.ScalarMappable(cmap=ColourMap)
    m.set_array(ChannelData.m_chi)
    plt.colorbar(m, cax=CAx,orientation='horizontal')
    plt.xlabel('$M_{\chi}$ m$^{0.64}$')
    #save output
    plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(BasinID) + "_LongProfMChi.png", dpi=300)

def PlotCHTAgainstChannelData(DataDirectory, FilenamePrefix, PlotDirectory, BasinID):
    """
    This function makes a plot of hilltop curavture against data
    from the channel segments.

    Args:
        BasinID: id of basin for analysis

    Returns: plot of CHT

    Author: FJC
    """
    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)

    # load the hillslopes data
    HillslopeData = ReadHillslopeData(DataDirectory, FilenamePrefix)
    # isolate basin data
    BasinChannelData = ChannelData[ChannelData.basin_key == BasinID]
    BasinJunctions = HillslopeData.BasinID.unique()
    BasinHillslopeData = HillslopeData[HillslopeData.BasinID == BasinJunctions[BasinID]]
    MinimumDistance = BasinChannelData.flow_distance.min()
    MaximumDistance = BasinChannelData.flow_distance.max()
    MaximumMChi = BasinChannelData.m_chi.max()

    # how many segments are we dealing with?
    Segments = BasinChannelData.segment_number.unique()

    # separate into main stem and trib data
    MainStemChannelData = BasinChannelData[BasinChannelData.source_key == 0]
    MainStemSegments = MainStemChannelData.segment_number.unique()

    # read in the terrace data
    TerraceData = ReadTerraceData(DataDirectory,FilenamePrefix)


    # set up the figure
    Fig = CreateFigure()
    Ax = plt.subplot(111)

    #choose colormap
    ColourMap = cm.viridis

    # loop through the channel data and get the CHT for this distance upstream.
    MainStemMeanCHT = []
    MainStemDist = []
    MainStemLH = []
    TribsMeanCHT = []
    TribsDist = []
    TribsLH = []
    TribsCHTData = []
    MainStemCHTData = []
    MainStemSTDev = []

    for i in range (0, len(Segments)):
        SegmentHillslopeData = BasinHillslopeData[BasinHillslopeData.StreamID == Segments[i]]
        SegmentChannelData = BasinChannelData[BasinChannelData.segment_number == Segments[i]]
        if SegmentHillslopeData.size != 0:
            if Segments[i] in MainStemSegments:
                MainStemMeanCHT.append(SegmentHillslopeData.Cht.mean())
                MainStemDist.append(SegmentChannelData.flow_distance.median()/1000)
                MainStemLH.append(SegmentHillslopeData.Lh)
                MainStemCHTData.append(np.array(SegmentHillslopeData.Cht))
                MainStemSTDev.append(SegmentHillslopeData.Cht.std())
            else:
                TribsMeanCHT.append(SegmentHillslopeData.Cht.mean())
                TribsDist.append(SegmentChannelData.flow_distance.median()/1000)
                TribsCHTData.append(np.array(SegmentHillslopeData.Cht))

    print TribsCHTData

    # # bin the data by distance upstream
    # bin_width = 5
    # n_bins = int((MaximumDistance-MinimumDistance)/bin_width)
    # bins = np.linspace(MinimumDistance,MaximumDistance,n_bins)
    # print bins
    #
    # n, dist_bins = np.histogram(TribsDist, bins=n_bins)
    # print dist_bins
    # s_TribsCHT, dist_bins = np.histogram(TribsDist, bins=n_bins, weights=TribsMeanCHT)
    # mean_CHTs = s_TribsCHT/n
    #
    # print mean_CHTs

    #Ax.violinplot(TribsCHTData, TribsDist, points=20, widths=1,showmeans=True, showextrema=False, showmedians=True)
    # now make the plot of the channel profile and the cht data
    #Ax.scatter(TribsDist, TribsMeanCHT, s=1, c='0.5')
    # Ax.scatter(bins/1000, mean_CHTs, s=1)
    Ax.errorbar(MainStemDist, MainStemMeanCHT, yerr=MainStemSTDev, mec='k',mfc='white',zorder=2, ecolor='k',fmt='o', ms=5)
    # Ax.violinplot(MainStemCHTData, MainStemDist, points=20, widths=1, showmeans=True, showmedians=True)
    #Ax.boxplot(MainStemCHTData, positions=MainStemDist, manage_xticks=False, widths=1, sym='')
    plt.xlabel('Distance from outlet ($km$)')
    plt.ylabel('Mean hilltop curvature ($m^{-1}$)')
    #Ax.set_xlim(0,50)
    #plt.text(-0.2,-0.3,'Basin ID ' + str(BasinID),transform = Ax.transAxes,color=[0.35,0.35,0.35])
    plt.tight_layout()

    #save output
    plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(BasinID) + "_CHT_flowdist.png", dpi=300)

def PlotEStarRStar(DataDirectory, FilenamePrefix, PlotDirectory, BasinID):
    """
    Makes a plot of E* against R* where the points are coloured by
    their distance from the outlet of the basin

    Author: FJC
    """
    import math
    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)

    # load the hillslopes data
    HillslopeData = ReadHillslopeData(DataDirectory, FilenamePrefix)

    # isolate basin data
    BasinChannelData = ChannelData[ChannelData.basin_key == BasinID]
    BasinJunctions = HillslopeData.BasinID.unique()
    BasinHillslopeData = HillslopeData[HillslopeData.BasinID == BasinJunctions[BasinID]]
    MinimumDistance = BasinChannelData.flow_distance.min()
    MaximumDistance = BasinChannelData.flow_distance.max()
    MaximumMChi = BasinChannelData.m_chi.max()

    # how many segments are we dealing with?
    Segments = BasinChannelData.segment_number.unique()

    # separate into main stem and trib data
    MainStemChannelData = BasinChannelData[BasinChannelData.source_key == 0]
    MainStemSegments = MainStemChannelData.segment_number.unique()

    # set up the figure
    Fig = CreateFigure()
    Ax = plt.subplot(111)

    #choose colormap
    ColourMap = cm.viridis

    # loop through the channel data and get the E* and R* for this distance upstream.
    MainStemEStar = []
    MainStemRStar = []
    MainStemDist = []
    TribsEStar = []
    TribsRStar = []
    TribsDist = []

    for i in range (0, len(Segments)):
        SegmentHillslopeData = BasinHillslopeData[BasinHillslopeData.StreamID == Segments[i]]
        SegmentChannelData = BasinChannelData[BasinChannelData.segment_number == Segments[i]]
        if SegmentHillslopeData.size != 0:
            if Segments[i] in MainStemSegments:
                MainStemEStar.append(SegmentHillslopeData.E_Star.mean())
                MainStemRStar.append(SegmentHillslopeData.R_Star.mean())
                MainStemDist.append(SegmentChannelData.flow_distance.median()/1000)
            else:
                TribsEStar.append(SegmentHillslopeData.E_Star.mean())
                TribsRStar.append(SegmentHillslopeData.R_Star.mean())
                TribsDist.append(SegmentChannelData.flow_distance.median()/1000)

    # get the model data for this E_Star
    ModelRStar = []
    for x in MainStemEStar:
        ModelRStar.append((1./x) * (np.sqrt(1.+(x*x)) - np.log(0.5*(1. + np.sqrt(1.+(x*x)))) - 1.))

    ModelEStar = [x for x,_ in sorted(zip(MainStemEStar,ModelRStar))]
    ModelRStar = [y for _,y in sorted(zip(MainStemEStar,ModelRStar))]


    Ax.plot(ModelEStar,ModelRStar, c='k')
    Ax.scatter(MainStemEStar,MainStemRStar,c=MainStemDist,s=10, edgecolors='k', lw=0.1,cmap=ColourMap)
    # Ax.scatter(TribsEStar,TribsRStar,c=TribsDist,s=10, edgecolors='k', lw=0.1,cmap=ColourMap)
    # Ax.set_xscale('log')
    # Ax.set_yscale('log')
    plt.xlabel('$E*$')
    plt.ylabel('$R*$')



    plt.subplots_adjust(left=0.18,right=0.85, bottom=0.2, top=0.9)
    CAx = Fig.add_axes([0.86,0.2,0.02,0.7])
    m = cm.ScalarMappable(cmap=ColourMap)
    m.set_array(MainStemDist)
    plt.colorbar(m, cax=CAx,orientation='vertical', label='Distance from outlet (km)')

    #Ax.set_ylim(-20,1)
    #plt.text(-0.2,-0.3,'Basin ID ' + str(BasinID),transform = Ax.transAxes,color=[0.35,0.35,0.35])
    #plt.tight_layout()

    #save output
    plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(BasinID) + "_EStar_RStar.png", dpi=300)
    plt.clf()

def PlotRStarDistance(DataDirectory, FilenamePrefix, PlotDirectory, BasinID):
    """
    Makes a plot of R* against distance from outlet

    Author: FJC
    """
    import math
    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)

    # load the hillslopes data
    HillslopeData = ReadHillslopeData(DataDirectory, FilenamePrefix)

    # isolate basin data
    BasinChannelData = ChannelData[ChannelData.basin_key == BasinID]
    BasinJunctions = HillslopeData.BasinID.unique()
    BasinHillslopeData = HillslopeData[HillslopeData.BasinID == BasinJunctions[BasinID]]
    MinimumDistance = BasinChannelData.flow_distance.min()
    MaximumDistance = BasinChannelData.flow_distance.max()
    MaximumMChi = BasinChannelData.m_chi.max()

    # how many segments are we dealing with?
    Segments = BasinChannelData.segment_number.unique()

    # separate into main stem and trib data
    MainStemChannelData = BasinChannelData[BasinChannelData.source_key == 0]
    MainStemSegments = MainStemChannelData.segment_number.unique()

    # set up the figure
    Fig = CreateFigure()
    Ax = plt.subplot(111)

    #choose colormap
    ColourMap = cm.viridis

    # loop through the channel data and get the E* and R* for this distance upstream.
    MainStemEStar = []
    MainStemRStar = []
    MainStemDist = []
    TribsEStar = []
    TribsRStar = []
    TribsDist = []

    for i in range (0, len(Segments)):
        SegmentHillslopeData = BasinHillslopeData[BasinHillslopeData.StreamID == Segments[i]]
        SegmentChannelData = BasinChannelData[BasinChannelData.segment_number == Segments[i]]
        if SegmentHillslopeData.size != 0:
            if Segments[i] in MainStemSegments:
                MainStemEStar.append(SegmentHillslopeData.E_Star.mean())
                MainStemRStar.append(SegmentHillslopeData.R_Star.mean())
                MainStemDist.append(SegmentChannelData.flow_distance.median()/1000)
            else:
                TribsEStar.append(SegmentHillslopeData.E_Star.mean())
                TribsRStar.append(SegmentHillslopeData.R_Star.mean())
                TribsDist.append(SegmentChannelData.flow_distance.median()/1000)

    # get the model data for this E_Star
    ModelRStar = []
    for x in MainStemEStar:
        ModelRStar.append((1./x) * (np.sqrt(1.+(x*x)) - np.log(0.5*(1. + np.sqrt(1.+(x*x)))) - 1.))

    ModelEStar = [x for x,_ in sorted(zip(MainStemEStar,ModelRStar))]
    ModelRStar = [y for _,y in sorted(zip(MainStemEStar,ModelRStar))]


    #Ax.plot(ModelEStar,ModelRStar, c='k')
    Ax.scatter(MainStemDist,MainStemRStar,c=MainStemEStar,s=10, edgecolors='k', lw=0.1,cmap=ColourMap)
    # Ax.scatter(TribsEStar,TribsRStar,c=TribsDist,s=10, edgecolors='k', lw=0.1,cmap=ColourMap)
    # Ax.set_xscale('log')
    # Ax.set_yscale('log')
    plt.xlabel('Distance from outlet (km)')
    plt.ylabel('$R*$')



    plt.subplots_adjust(left=0.18,right=0.85, bottom=0.2, top=0.9)
    CAx = Fig.add_axes([0.86,0.2,0.02,0.7])
    m = cm.ScalarMappable(cmap=ColourMap)
    m.set_array(MainStemDist)
    plt.colorbar(m, cax=CAx,orientation='vertical', label='$E*$')

    #Ax.set_ylim(-20,1)
    #plt.text(-0.2,-0.3,'Basin ID ' + str(BasinID),transform = Ax.transAxes,color=[0.35,0.35,0.35])
    #plt.tight_layout()

    #save output
    plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(BasinID) + "_RStar_Dist.png", dpi=300)
    plt.clf()

def PlotLHDistance(DataDirectory, FilenamePrefix, PlotDirectory, BasinID):
    """
    Makes a plot of Lh against distance from outlet

    Author: FJC
    """

    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)

    # load the hillslopes data
    HillslopeData = ReadHillslopeData(DataDirectory, FilenamePrefix)

    # isolate basin data
    BasinChannelData = ChannelData[ChannelData.basin_key == BasinID]
    BasinJunctions = HillslopeData.BasinID.unique()
    BasinHillslopeData = HillslopeData[HillslopeData.BasinID == BasinJunctions[BasinID]]
    MinimumDistance = BasinChannelData.flow_distance.min()
    MaximumDistance = BasinChannelData.flow_distance.max()
    MaximumMChi = BasinChannelData.m_chi.max()

    # how many segments are we dealing with?
    Segments = BasinChannelData.segment_number.unique()

    # separate into main stem and trib data
    MainStemChannelData = BasinChannelData[BasinChannelData.source_key == 0]
    MainStemSegments = MainStemChannelData.segment_number.unique()

    # set up the figure
    Fig = CreateFigure()
    Ax = plt.subplot(111)

    #choose colormap
    ColourMap = cm.viridis

    # loop through the channel data and get the E* and R* for this distance upstream.
    MainStemLh = []
    MainStemDist = []
    TribsLh = []
    TribsRStar = []
    TribsDist = []

    for i in range (0, len(Segments)):
        SegmentHillslopeData = BasinHillslopeData[BasinHillslopeData.StreamID == Segments[i]]
        SegmentChannelData = BasinChannelData[BasinChannelData.segment_number == Segments[i]]
        if SegmentHillslopeData.size != 0:
            if Segments[i] in MainStemSegments:
                MainStemLh.append(SegmentHillslopeData.Lh.mean())
                MainStemDist.append(SegmentChannelData.flow_distance.median()/1000)
            else:
                TribsLh.append(SegmentHillslopeData.Lh.mean())
                TribsDist.append(SegmentChannelData.flow_distance.median()/1000)


    #Ax.plot(ModelEStar,ModelRStar, c='k')
    Ax.scatter(MainStemDist,MainStemLh, c = 'k', s=10, edgecolors='k', lw=0.1,cmap=ColourMap)
    # Ax.scatter(TribsEStar,TribsRStar,c=TribsDist,s=10, edgecolors='k', lw=0.1,cmap=ColourMap)
    # Ax.set_xscale('log')
    # Ax.set_yscale('log')
    plt.xlabel('Distance from outlet (km)')
    plt.ylabel('$L_h$')



    # plt.subplots_adjust(left=0.18,right=0.85, bottom=0.2, top=0.9)
    # CAx = Fig.add_axes([0.86,0.2,0.02,0.7])
    # m = cm.ScalarMappable(cmap=ColourMap)
    # m.set_array(MainStemDist)
    # plt.colorbar(m, cax=CAx,orientation='vertical', label='$E*$')

    #Ax.set_ylim(-20,1)
    #plt.text(-0.2,-0.3,'Basin ID ' + str(BasinID),transform = Ax.transAxes,color=[0.35,0.35,0.35])
    plt.tight_layout()

    #save output
    plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(BasinID) + "_Lh_dist.png", dpi=300)
    plt.clf()

def PlotHillslopeDataVsDistance(DataDirectory, FilenamePrefix, PlotDirectory, BasinID):
    """
    This function makes some composite plots of the hillslope data vs
    distance upstream from the outlet

    Author: FJC
    """
    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)
    #print BasinChannelData

    # load the hillslopes data
    HillslopeData = ReadHillslopeData(DataDirectory, FilenamePrefix)

    # isolate basin data
    BasinChannelData = ChannelData[ChannelData.basin_key == BasinID]
    BasinJunctions = HillslopeData.BasinID.unique()
    BasinHillslopeData = HillslopeData[HillslopeData.BasinID == BasinJunctions[BasinID]]
    MinimumDistance = BasinChannelData.flow_distance.min()
    MaximumDistance = BasinChannelData.flow_distance.max()
    MaximumMChi = BasinChannelData.m_chi.max()

    # how many segments are we dealing with?
    Segments = BasinChannelData.segment_number.unique()

    # separate into main stem and trib data
    MainStemChannelData = BasinChannelData[BasinChannelData.source_key == 0]
    MainStemSegments = MainStemChannelData.segment_number.unique()

    #choose colormap
    ColourMap = cm.viridis

    # loop through the channel data and get the E* and R* for this distance upstream.
    DistanceFromOutlet = []
    Lh = []
    Lh_std = []
    Cht = []
    Cht_std = []
    R_star = []
    R_star_std = []
    E_star = []
    E_star_std = []
    M_chi = []
    M_chi_std = []

    for i in range (0, len(Segments)):
        SegmentHillslopeData = BasinHillslopeData[BasinHillslopeData.StreamID == Segments[i]]
        SegmentChannelData = BasinChannelData[BasinChannelData.segment_number == Segments[i]]
        if SegmentHillslopeData.size != 0:
            # only analysing segments directly connected to the main stem
            if Segments[i] in MainStemSegments:
                DistanceFromOutlet.append(SegmentChannelData.flow_distance.median()/1000)
                Lh.append(SegmentHillslopeData.Lh.mean())
                Lh_std.append(SegmentHillslopeData.Lh.std())
                Cht.append(SegmentHillslopeData.Cht.mean())
                Cht_std.append(SegmentHillslopeData.Cht.std())
                R_star.append(SegmentHillslopeData.R_Star.mean())
                R_star_std.append(SegmentHillslopeData.R_Star.std())
                E_star.append(SegmentHillslopeData.E_Star.mean())
                E_star_std.append(SegmentHillslopeData.E_Star.std())
                M_chi.append(SegmentChannelData.m_chi.mean())
                M_chi_std.append(SegmentChannelData.m_chi.std())

    # set up the figure
    fig, ax = plt.subplots(nrows = 4, ncols=1, sharex=True, figsize=(6,7))
    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=0)

    # plot the hillslope length
    print DistanceFromOutlet
    ax[0].errorbar(DistanceFromOutlet,Lh,yerr=Lh_std,fmt='o', ecolor='0.5',markersize=5,mec='k')
    ax[0].set_ylabel('$L_h$')

    #plot the cht
    ax[1].errorbar(DistanceFromOutlet,Cht,yerr=Cht_std,fmt='o', ecolor='0.5',markersize=5,mfc='red',mec='k')
    ax[1].set_ylabel('$C_{HT}$')

    #plot the R*
    ax[2].errorbar(DistanceFromOutlet,R_star,yerr=R_star_std,fmt='o', ecolor='0.5',markersize=5,mfc='orange',mec='k')
    ax[2].set_ylabel('$R*$')

    #plot the R*
    ax[3].errorbar(DistanceFromOutlet,M_chi,yerr=M_chi_std,fmt='o', ecolor='0.5',markersize=5,mfc='purple',mec='k')
    ax[3].set_ylabel('$k_{sn}$')

    # set the axes labels
    ax[3].set_xlabel('Distance from outlet (km)')
    plt.tight_layout()

    #save output
    plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(BasinID) + "_hillslopes_distance.png", dpi=300)
    #plt.clf()

def PlotHillslopeDataWithBasins(DataDirectory,FilenamePrefix,PlotDirectory):
    """
    Function to make plots of hillslope data vs basin id.
    Martin probably has nice versions of this but I need something
    quick for my poster

    Author: FJC
    """

    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)
    #print BasinChannelData

    # load the hillslopes data
    HillslopeData = ReadHillslopeData(DataDirectory, FilenamePrefix)

    basin_dict = MapBasinKeysToJunctions(DataDirectory,FilenamePrefix)
    basin_keys = basin_dict.keys()

    median_cht = []
    cht_lower_err = []
    cht_upper_err = []

    median_Lh = []
    Lh_lower_err = []
    Lh_upper_err = []

    median_Rstar = []
    Rstar_lower_err = []
    Rstar_upper_err = []

    median_Estar = []
    Estar_lower_err = []
    Estar_upper_err = []

    median_mchi = []
    mchi_lower_err = []
    mchi_upper_err = []

    for key, jn in basin_dict.iteritems():
        BasinHillslopeData = HillslopeData[HillslopeData.BasinID == jn]
        BasinChannelData = ChannelData[ChannelData.basin_key == key]

        # now get all the hillslope data for this basin
        this_median = abs(BasinHillslopeData.Cht.median())
        median_cht.append(this_median)
        cht_lowerP = np.percentile(BasinHillslopeData.Cht, 16)
        cht_upperP = np.percentile(BasinHillslopeData.Cht, 84)
        cht_lower_err.append(this_median-abs(cht_upperP)) # these are the opposite way round because
        cht_upper_err.append(abs(cht_lowerP)-this_median) # I am inverting the axis to show positive Chts

        this_median = BasinHillslopeData.Lh.median()
        median_Lh.append(this_median)
        Lh_lowerP = np.percentile(BasinHillslopeData.Lh, 16)
        Lh_upperP = np.percentile(BasinHillslopeData.Lh, 84)
        Lh_lower_err.append(this_median-Lh_lowerP)
        Lh_upper_err.append(Lh_upperP-this_median)

        this_median = BasinHillslopeData.R_Star.median()
        median_Rstar.append(this_median)
        Rstar_lowerP = np.percentile(BasinHillslopeData.R_Star, 16)
        Rstar_upperP = np.percentile(BasinHillslopeData.R_Star, 84)
        Rstar_lower_err.append(this_median-Rstar_lowerP)
        Rstar_upper_err.append(Rstar_upperP-this_median)

        this_median = BasinHillslopeData.E_Star.median()
        median_Estar.append(this_median)
        Estar_lowerP = np.percentile(BasinHillslopeData.E_Star, 16)
        Estar_upperP = np.percentile(BasinHillslopeData.E_Star, 84)
        Estar_lower_err.append(this_median-Estar_lowerP)
        Estar_upper_err.append(Estar_upperP-this_median)


        # get the channel data
        this_median = BasinChannelData.m_chi.median()
        median_mchi.append(this_median)
        mchi_lowerP = np.percentile(BasinChannelData.m_chi, 16)
        mchi_upperP = np.percentile(BasinChannelData.m_chi, 84)
        mchi_lower_err.append(this_median-mchi_lowerP)
        mchi_upper_err.append(mchi_upperP-this_median)

    # set up the figure
    fig, ax = plt.subplots(nrows = 6, ncols=1, sharex=True, figsize=(6,10), facecolor='white')
    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=0)

    # plot the hillslope length
    ax[0].errorbar(basin_keys,median_Lh,yerr=[Lh_lower_err, Lh_upper_err],fmt='o', ecolor='0.5',markersize=5,mec='k')
    ax[0].set_ylabel('$L_h$')

    #plot the cht
    ax[1].errorbar(basin_keys,median_cht,yerr=[cht_lower_err, cht_upper_err],fmt='o', ecolor='0.5',markersize=5,mfc='red',mec='k')
    ax[1].set_ylabel('$C_{HT}$')

    #plot the R*
    ax[2].errorbar(basin_keys,median_Rstar,yerr=[Rstar_lower_err, Rstar_upper_err],fmt='o', ecolor='0.5',markersize=5,mfc='orange',mec='k')
    ax[2].set_ylabel('$R*$')

    #plot the Mchi
    ax[3].errorbar(basin_keys,median_mchi,yerr=[mchi_lower_err, mchi_upper_err],fmt='o', ecolor='0.5',markersize=5,mfc='purple',mec='k')
    ax[3].set_ylabel('$k_{sn}$')

    # read the uplift data in
    # read in the csv
    uplift_df = pd.read_csv(DataDirectory+'MTJ_basin_uplift.csv')
    dd_df = pd.read_csv(DataDirectory+FilenamePrefix+'_basin_dd.csv')

    # get the drainage density
    drainage_density = dd_df['drainage_density']
    ax[4].scatter(basin_keys, drainage_density, c='k', edgecolors='k', s=20)
    ax[4].set_ylim(np.min(drainage_density)-0.001, np.max(drainage_density)+0.001)
    ax[4].set_ylabel('Drainage density (m/m$^2$)')

    # get the data
    uplift_rate = uplift_df['Uplift_rate']
    ax[5].plot(basin_keys, uplift_rate, c='k', ls='--')
    ax[5].set_ylabel('Uplift rate (mm/yr)')

    # set the axes labels
    ax[5].set_xlabel('Basin ID')
    plt.xticks(np.arange(min(basin_keys), max(basin_keys)+1, 1))
    plt.tight_layout()

    #save output
    plt.savefig(PlotDirectory+FilenamePrefix +"_basin_hillslope_data.png", dpi=300)
    plt.clf()

    output_list = [('basin_keys', basin_keys),
                   ('uplift_rate', uplift_rate),
                   ('Lh_median', median_Lh),
                   ('Lh_lower_err', Lh_lower_err),
                   ('Lh_upper_err', Lh_upper_err),
                   ('cht_median', median_cht),
                   ('cht_lower_err', cht_lower_err),
                   ('cht_upper_err', cht_upper_err),
                   ('Rstar_median', median_Rstar),
                   ('Rstar_lower_err', Rstar_lower_err),
                   ('Rstar_upper_err', Rstar_upper_err),
                   ('Estar_median', median_Estar),
                   ('Estar_lower_err', Estar_lower_err),
                   ('Estar_upper_err', Estar_upper_err),
                   ('mchi_median', median_mchi),
                   ('mchi_lower_err', mchi_lower_err),
                   ('mchi_upper_err', mchi_upper_err),
                   ('drainage_density', drainage_density)]

    # write output to csv
    OutDF = pd.DataFrame.from_items(output_list)
    csv_outname = PlotDirectory+FilenamePrefix+'_basin_hillslope_data.csv'
    OutDF.to_csv(csv_outname,index=False)

def PlotHillslopeDataWithBasinsFromCSV(DataDirectory, FilenamePrefix):
    """
    Function to make same plot as above but read data from csv
    with the extension '_hillslope_means.csv'

    Author: FJC
    Will clean this up after AGU
    """
    input_csv = DataDirectory+FilenamePrefix+'_hillslope_means.csv'
    df = pd.read_csv(input_csv)

    # set up the figure
    fig, ax = plt.subplots(nrows = 5, ncols=1, sharex=True, figsize=(6,10))
    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=0)

    # plot the hillslope length
    ax[0].errorbar(df['basin_keys'],df['Lh_mean'],yerr=df['Lh_std'],fmt='o', ecolor='0.5',markersize=5,mec='k')
    ax[0].set_ylabel('$L_h$')

    #plot the cht
    ax[1].errorbar(df['basin_keys'],df['Cht_mean'],yerr=df['Cht_std'],fmt='o', ecolor='0.5',markersize=5,mfc='red',mec='k')
    ax[1].set_ylabel('$C_{HT}$')

    #plot the R*
    ax[2].errorbar(df['basin_keys'],df['R_star_mean'],yerr=df['R_star_std'],fmt='o', ecolor='0.5',markersize=5,mfc='orange',mec='k')
    ax[2].set_ylabel('$R*$')

    #plot the Mchi
    ax[3].errorbar(df['basin_keys'],df['M_chi_mean'],yerr=df['M_chi_std'],fmt='o', ecolor='0.5',markersize=5,mfc='purple',mec='k')
    ax[3].set_ylabel('$k_{sn}$')

    ax[4].plot(df['basin_keys'], df['uplift_rate'], c='k', ls='--')
    ax[4].set_ylabel('Uplift rate (mm/yr)')

    # set the axes labels
    ax[4].set_xlabel('Basin ID')
    plt.xticks(np.arange(min(df['basin_keys']), max(df['basin_keys'])+1, 1))
    plt.tight_layout()

    #save output
    plt.savefig(DataDirectory+FilenamePrefix +"_mean_hillslope_data.png", dpi=300)
    plt.clf()

def PlotBasinDataAgainstUplift(DataDirectory, FilenamePrefix, PlotDirectory):
    """
    Function to plot mean basin data as a function of uplift rate
    using the hillslope means csv file

    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): the file name prefix
        PlotDirectory (str): directory to save the plots

    Author: FJC
    """

    input_csv = PlotDirectory+FilenamePrefix+'_basin_hillslope_data.csv'
    df = pd.read_csv(input_csv)
    print df['uplift_rate']

    # set up the figure
    fig, ax = plt.subplots(nrows=3, ncols=2, sharex=True, figsize=(10,10))
    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=0)
    ax = ax.ravel()
    print ax[0]

    # plot the hillslope length
    ax[0].errorbar(df['uplift_rate'],df['Lh_median'],yerr=[df['Lh_lower_err'], df['Lh_upper_err']],fmt='o', ecolor='0.5',markersize=5,mec='k', mfc='blue')
    ax[0].set_ylabel('Hillslope length ($L_h$)')

    #plot the cht
    ax[1].errorbar(df['uplift_rate'],df['cht_median'],yerr=[df['cht_lower_err'], df['cht_upper_err']],fmt='o', ecolor='0.5',markersize=5,mfc='red',mec='k')
    ax[1].set_ylabel('Hilltop curvature ($C_{HT}$)')

    #plot the R*
    ax[2].errorbar(df['uplift_rate'],df['Rstar_median'],yerr=[df['Rstar_lower_err'], df['Rstar_upper_err']],fmt='o', ecolor='0.5',markersize=5,mfc='orange',mec='k')
    ax[2].set_ylabel('Relief ($R*$)')

    #plot the Mchi
    ax[3].errorbar(df['uplift_rate'],df['mchi_median'],yerr=[df['mchi_lower_err'], df['mchi_upper_err']],fmt='o', ecolor='0.5',markersize=5,mfc='purple',mec='k')
    ax[3].set_ylabel('Channel steepness ($k_{sn}$)')

    ax[4].scatter(df['uplift_rate'], df['drainage_density'], c='k')
    ax[4].set_ylabel('Drainage densiy (m/m$^2$)')
    ax[4].set_ylim(np.min(df['drainage_density'])-0.001, np.max(df['drainage_density'])+0.001)

    # set the axes labels
    ax[4].set_xlabel('Uplift rate (mm/yr)')
    #plt.xticks(np.arange(min(df['basin_keys']), max(df['basin_keys'])+1, 1))
    plt.tight_layout()

    #save output
    plt.savefig(PlotDirectory+FilenamePrefix +"_mean_data_fxn_U.png", dpi=300)
    plt.clf()

def PlotKsnAgainstRStar(DataDirectory, FilenamePrefix, PlotDirectory):
    """
    Function to plot median Ksn against R* for a series of basins

    Author: FJC
    """

    input_csv = PlotDirectory+FilenamePrefix+'_basin_hillslope_data.csv'
    df = pd.read_csv(input_csv)

    # linregress
    slope, intercept, r_value, p_value, std_err = stats.linregress(df['mchi_median'],df['Rstar_median'])
    print (slope, intercept, r_value, p_value)
    x = np.linspace(0, 200, 100)
    new_y = slope*x + intercept

    # set up the figure
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=(5,5))

    ax.scatter(df['mchi_median'], df['Rstar_median'], c=df['basin_keys'], s=25, edgecolors='k', zorder=100, cmap=cm.viridis)
    ax.errorbar(df['mchi_median'], df['Rstar_median'], xerr=[df['mchi_lower_err'], df['mchi_upper_err']], yerr=[df['Rstar_lower_err'], df['Rstar_upper_err']], fmt='o', ecolor='0.5',markersize=1,mfc='white',mec='k')
    ax.text(0.55, 0.1, '$y = $'+str(np.round(slope,4))+'$x + $'+str(np.round(intercept,2))+'\n$R^2 = $'+str(np.round(r_value,2))+'\n$p = $'+str(p_value), fontsize=9, color='black', transform=ax.transAxes)
    ax.plot(x, new_y, c='0.5', ls='--')
    ax.set_xlim(0,200)

    ax.set_xlabel('$k_{sn}$')
    ax.set_ylabel('$R*$')

    plt.subplots_adjust(left=0.15,right=0.85, bottom=0.1, top=0.95)
    CAx = fig.add_axes([0.87,0.1,0.02,0.85])
    m = cm.ScalarMappable(cmap=cm.viridis)
    m.set_array(df['basin_keys'])
    plt.colorbar(m, cax=CAx,orientation='vertical', label='Basin key')

    #plt.tight_layout()

    #save output
    plt.savefig(PlotDirectory+FilenamePrefix +"_ksn_vs_rstar.png", dpi=300)
    plt.clf()

def PlotEStarRStarBasins(DataDirectory, FilenamePrefix, PlotDirectory):
    """
    Function to make an E*R* plot for a series of drainage basins
    Author: FJC
    """
    input_csv = PlotDirectory+FilenamePrefix+'_basin_hillslope_data.csv'
    df = pd.read_csv(input_csv)

    # steady state model
    x = np.linspace(0,50,1000)
    predicted_Rstar = ((1. / x) * (np.sqrt(1. + (x * x)) -
            np.log(0.5 * (1. + np.sqrt(1. + (x * x)))) - 1.))

    # set up the figure
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=(5,5))

    ax.scatter(df['Estar_median'], df['Rstar_median'], c=df['basin_keys'], s=25, edgecolors='k', zorder=100, cmap=cm.viridis, label=None)
    ax.errorbar(df['Estar_median'], df['Rstar_median'], xerr=[df['Estar_lower_err'], df['Estar_upper_err']], yerr=[df['Rstar_lower_err'], df['Rstar_upper_err']], fmt='o', ecolor='0.5',markersize=1,mfc='white',mec='k',label=None)
    #ax.text(0.55, 0.1, '$y = $'+str(np.round(slope,4))+'$x + $'+str(np.round(intercept,2))+'\n$R^2 = $'+str(np.round(r_value,2))+'\n$p = $'+str(p_value), fontsize=9, color='black', transform=ax.transAxes)
    ax.plot(x, predicted_Rstar, c='0.5', ls='--', label='Steady state model')
    ax.set_xlim(0,50)

    ax.set_xlabel('$E*$')
    ax.set_ylabel('$R*$')

    ax.legend(loc='lower right')

    # ax.set_xscale('log')
    # ax.set_yscale('log')

    plt.subplots_adjust(left=0.15,right=0.85, bottom=0.1, top=0.95)
    CAx = fig.add_axes([0.87,0.1,0.02,0.85])
    m = cm.ScalarMappable(cmap=cm.viridis)
    m.set_array(df['basin_keys'])
    plt.colorbar(m, cax=CAx,orientation='vertical', label='Basin key')

    #plt.tight_layout()

    #save output
    plt.savefig(PlotDirectory+FilenamePrefix +"_estar_vs_rstar.png", dpi=300)
    plt.clf()

def PlotJunctionAnglesAgainstBasinID(DataDirectory, FilenamePrefix, PlotDirectory):
    """
    Function to make plot of the junction angles vs. stream order for different
    basins.
    May move this to a different plotting script eventually.

    FJC 08/03/18
    """
    df = pd.read_csv(DataDirectory+FilenamePrefix+'_junc_angles.csv')
    basin_dict = MapBasinKeysToJunctions(DataDirectory, FilenamePrefix)
    print basin_dict

    # set up the figure
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=(12,5))

    # separate the series by stream order
    stream_orders = df['StreamOrder'].unique()

    # CHANGE THIS LATER - Fiona
    basin_jns = df['BasinJunction'].unique()
    n_jns = len(basin_jns)
    basin_keys = np.arange(0, n_jns, 1)
    basin_dict = dict(zip(basin_keys,basin_jns))
    temp_df = pd.DataFrame(list(basin_dict.iteritems()), columns=['BasinKey','BasinJunction'])
    df = df.merge(temp_df, left_on="BasinJunction", right_on = "BasinJunction")

    for SO in stream_orders:
        # mask to this SO
        this_df = df[df['StreamOrder'] == SO]

        # get da basin keys - this way SHOULD work when new m/n is finished
        # FIONA - REMEMBER TO CHANGE THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # temp_df = pd.DataFrame(list(basin_dict.iteritems()), columns=['BasinKey','BasinJunction'])
        # this_df = this_df.merge(temp_df, left_on="BasinJunction", right_on = "BasinJunction")

        # these_juncs = this_df['BasinJunction']
        # these_keys =

        ax.errorbar(this_df['BasinKey'], this_df['Median'], yerr=[this_df['Median']-this_df['25th_percentile'], this_df['75th_percentile']-this_df['Median']], fmt='o', label='Stream order = '+str(SO))

    ax.set_xlabel('Basin key')
    ax.set_ylabel('Junction angle (deg)')

    ax.legend(loc='upper right', bbox_to_anchor=(1.25,0.8))
    plt.subplots_adjust(left=0.1,right=0.8)
    plt.xticks(np.arange(min(df['BasinKey']), max(df['BasinKey'])+1, 1))

    #save output
    plt.savefig(PlotDirectory+FilenamePrefix +"_junction_angles.png", dpi=300)
    plt.clf()


def PlotHillslopeTraces(DataDirectory, FilenamePrefix, PlotDirectory, CustomExtent=[-9999],FigSizeFormat="epsl"):
    """
    Function to plot a hillshade image with hilltops, hillslope traces and the channel network superimposed.


    MDH
    """


    # Save the figure
    ImageName = PlotDirectory+"bolinas_traces.png"
    print(ImageName)
    FigWidth_Inches = Get_FigWidth_Inches(FigSizeFormat)

    HillshadeName = FilenamePrefix+"_hs.bil"

    # create the map figure
    MF = MapFigure(HillshadeName, DataDirectory, coord_type="UTM_km", colourbar_location='None')

    #customise the extent of the plot if required
    if len(CustomExtent) == 4:
      xmin = CustomExtent[0]
      xmax = CustomExtent[1]
      ymin = CustomExtent[2]
      ymax = CustomExtent[3]
      MF.SetCustomExtent(xmin,xmax,ymin,ymax)

    # add hilltops
    HilltopPointsDF = ReadHillslopeData(DataDirectory, FilenamePrefix)
    HilltopPoints = LSDP.LSDMap_PointData(HilltopPointsDF, data_type = "pandas", PANDEX = True)
    MF.add_point_data(HilltopPoints,alpha=1,zorder=100,unicolor=[0.8,0,0],manual_size=1)

    # add channel heads
    ChannelHeadsDF = pd.read_csv(DataDirectory+FilenamePrefix+"_CH_wiener_nodeindices_for_Arc.csv")
    ChannelHeadPoints = LSDP.LSDMap_PointData(ChannelHeadsDF, data_type = "pandas", PANDEX = True)
    MF.add_point_data(ChannelHeadPoints,alpha=0.5,zorder=100,unicolor="blue",manual_size=5)

    # add channels
    ChannelDF = Helper.ReadChiDataMapCSV(DataDirectory,FilenamePrefix)
    ChannelPoints = LSDP.LSDMap_PointData(ChannelDF, data_type = "pandas", PANDEX = True)
    MF.add_point_data(ChannelPoints,show_colourbar="False", scale_points=True, max_point_size = 2.5, min_point_size = 0.5, column_for_scaling='drainage_area',alpha=0.5,zorder=90)

    # add hillslope traces
    ThinningFactor=1
    HillslopeTracesShp = DataDirectory+FilenamePrefix+"_hillslope_traces.shp"
    if os.path.exists(HillslopeTracesShp) == False:
      WriteHillslopeTracesShp(DataDirectory,FilenamePrefix,ThinningFactor,CustomExtent)

    MF.add_line_data(DataDirectory+FilenamePrefix+"_hillslope_traces.shp",zorder=80,alpha=0.5,linewidth=0.2)

    #finalise and save figure
    MF.SetRCParams(label_size=8)
    MF.save_fig(fig_width_inches = FigWidth_Inches, FigFileName = ImageName, FigFormat="png", Fig_dpi = 300)
