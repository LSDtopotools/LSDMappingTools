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
from geopandas import GeoDataFrame
import pandas as pd
import numpy as np
from shapely.geometry import LineString, shape, Point, MultiPolygon, Polygon
import numpy.ma as ma

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
from descartes import PolygonPatch

# plotting tools for using LSDMapFigure

#=============================================================================
# PRELIMINARY FUNCTIONS
#=============================================================================

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

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = 10
    #rcParams['text.usetex'] = True

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
    to a pandas dataframe

    Args:
        DataDirectory: the data directory
        FilenamePrefix: the file name prefix

    Returns:
        pandas dataframe with data from the csv file

    Author: MDH
    """
    # get the csv filename
    Suffix = '_MChiSegmented.csv'
    Filename = FilenamePrefix+Suffix

    # read in the dataframe using pandas
    ChannelData = pd.read_csv(DataDirectory+Filename)

    # If there is no chi values due to threshold then chi will be -9999
    # throw out these segments
    Segments2Remove = ChannelData[ChannelData.chi == -9999].segment_number.unique()
    ChannelData = ChannelData[~ChannelData.segment_number.isin(Segments2Remove)]

    #return the hillslope data
    return ChannelData

def ReadHillslopeTraces(DataDirectory, FilenamePrefix):
    """
    This function reads in the file with the suffix '_hillslope_traces.csv'
    and creates a geopandas GeoDataFrame

    Args:
        DataDirectory: the data directory
        FilenamePrefix: the file name prefix

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
    geometry = [Point(xy) for xy in zip(df.Longitude, df.Latitude)]
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

def WriteHillslopeTracesShp(DataDirectory,FilenamePrefix):
    """
    This function writes a shapefile of hillslope traces

    Args:
        DataDirectory: the data directory
        FilenamePrefix: the file name prefix

    Author: MDH
    """

    #read the raw data to geodataframe
    geo_df = ReadHillslopeTraces(DataDirectory,FilenamePrefix)
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

    # setup the figure
    Fig = CreateFigure()

    #choose colormap
    ColourMap = cm.viridis

    # Get the data columns for plotting
    for i in range(0, len(Segments)):
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
    m = cm.ScalarMappable(cmap=ColourMap)
    m.set_array(ChannelData.m_chi)
    plt.colorbar(m, cax=CAx,orientation='horizontal')
    plt.xlabel('$M_{\chi}$ m$^{0.64}$')
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
``
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

def PlotCHTAgainstChannelData(DataDirectory, FilenamePrefix, PlotDirectory, BasinID):
    """
    Makes a plot of E* against R* where the points are coloured by
    their distance from the outlet of the basin

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
                MainStemEStar.append(SegmentHillslopeData.E_star.mean())
                MainStemRStar.append(SegmentHillslopeData.R_star.mean())
                MainStemDist.append(SegmentChannelData.flow_distance.median()/1000)
            else:
                TribsEStar.append(SegmentHillslopeData.E_star.mean())
                TribsRStar.append(SegmentHillslopeData.R_star.mean())
                TribsDist.append(SegmentChannelData.flow_distance.median()/1000)
