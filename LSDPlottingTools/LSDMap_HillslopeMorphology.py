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

# 3d projection
from mpl_toolkits.mplot3d import Axes3D

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
from matplotlib import cm, gridspec
from shapely.ops import cascaded_union
import matplotlib.colors as colors
from descartes import PolygonPatch

# plotting tools for using LSDMapFigure
import LSDPlottingTools as LSDP
from LSDMapFigure import PlottingHelpers as Helper
from LSDMapFigure.PlottingRaster import MapFigure

#=============================================================================
# PRELIMINARY FUNCTIONS
#=============================================================================

def Get_FigWidth_Inches(FigSizeFormat="default"):
    """
    This function gets the figure width in inches for different formats

    Args:
        FigSizeFormat: the figure size format according to journal for which the figure is intended
            values are geomorphology,ESURF, ESPL, EPSL, JGR, big
            default is ESURF

    Returns:
        figure width in inches. Inches is required since matplotlib was written by unreformed yanks.

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
    rcParams['font.sans-serif'] = ['Liberation Sans']
    rcParams['font.size'] = 8
    rcParams['text.usetex'] = False

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

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the name of the DEM

    Returns:
        A dictionary with the basin key as the key and the junction as the value

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

    print(basin_dict)
    return basin_dict

def WriteHillslopeTracesShp(DataDirectory,FilenamePrefix,ThinningFactor=1, CustomExtent=[-9999]):
    """
    This function writes a shapefile of hillslope traces

    Args:
        DataDirectory: the data directory
        FilenamePrefix: the file name prefix
        ThinningFactor (int): This determines how many of the traces are discarded. Higher numbers mean more traces are discarded
        CustomExtent (list): if this is [-9999] the extent is grabbed from the raster. Otherwise you give it a 4 element list giving the extent of the area of interest.

    Returns: None, but writes a shapefile

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

def RemoveSmallSegments(BasinHillslopeData, n_traces=50):
    """
    Remove hilltop segments with less than a specified number of traces in a basin

    Args:
        BasinHillslopeData (pandas dataframe): The dataframe containing the hillslope data (ridgelines). You get this using the ReadHillslopeData function
        n_traces (int) the minimum number of traces

    Author: FJC
    """
    # remove segments shorter than the threshold length
    BasinHillslopeData = BasinHillslopeData.groupby('StreamID').filter(lambda x: x.shape[0] > n_traces)
    return BasinHillslopeData


#---------------------------------------------------------------------------------#
# ANALYSIS FUNCTIONS
#---------------------------------------------------------------------------------#
def DetermineSc(DataDirectory,FilenamePrefix,PlotDirectory):
    """
    Determines the critical slope following Grieve et al. 2016 How Long is a Hillslope?

    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): the file name prefix
        PlotDirectory (str): The directory into which the plots are saved

    MDH

    """

    # load the hillslopes data and isolate segments in basin
    HillslopeData = ReadHillslopeData(DataDirectory, FilenamePrefix)

    # get all the segments
    Segments = HillslopeData.StreamID.unique()

    Lh_segments = np.zeros(len(Segments))
    R_segments = np.zeros(len(Segments))

    # For each segment get collect hillslope data
    # record the number of traces in each segment into a new dataframe
    Data = pd.DataFrame(columns=['SegmentNo','Lh','LhLower','LhUpper','R','RLower','RUpper','NTraces'])

    for i in range(0,len(Segments)):

        #Get segment hillslope dataMChi,FlowLength,SegmentLength,
        SegmentHillslopeData = HillslopeData[HillslopeData.StreamID == float(Segments[i])]

        #hillslopes
        Lh = SegmentHillslopeData.Lh.quantile(0.5)
        LhUpper = SegmentHillslopeData.Lh.quantile(0.75)
        LhLower = SegmentHillslopeData.Lh.quantile(0.25)

        R = SegmentHillslopeData.R.quantile(0.5)
        RLower = SegmentHillslopeData.R.quantile(0.25)
        RUpper = SegmentHillslopeData.R.quantile(0.75)

        NTraces = SegmentHillslopeData.size

        #add to data frame
        Data.loc[i] = [Segments[i],Lh,LhLower,LhUpper,R,RLower,RUpper,NTraces]

    # remove rows with no data (i.e. no hillslope traces)
    Data = Data.dropna()
    Data = Data[Data.NTraces > 50]

    # plot theoretical relationship
    LH = np.arange(0.,200.,1.)
    #Sc = 0.8
    #pr = 2400 #kg/m3
    #ps = 1400 #kg/m3
    #K = 0.01 #m2/y
    #E = 0.2 #m/y
    #k = (pr*E)/(2.*ps*K)

    #R = (Sc * (-1. + np.sqrt(1 + k**2. * LH**2.) + np.log(3.) - np.log(2. + np.sqrt(1. + k**2. * LH**2.))))/k
    #R = LH*Sc

    # declare colour map
    #ColourMap = cm.plasma_r

    #Create a figure instance for plotting Length vs Relief
    CreateFigure(AspectRatio=3.)

    # setup subplots
    gs = gridspec.GridSpec(1, 2, width_ratios=[2.5, 1])

    a0 = plt.subplot(gs[0])

    # plot the raw data
    a0.plot(HillslopeData.Lh,HillslopeData.R,'.',ms=2,color=[0.8,0.8,0.8])

    # plot the segmented data
    # Error bars with colours but faded (alpha)
    for i, row in Data.iterrows():
        #LhErr = np.array([[row.LhLower],[row.LhUpper]])
        #RErr = np.array([[row.RLower],[row.RUpper]])
        #plt.plot([row.Lh,row.Lh],RErr,'-', lw=1, color=[0.25,0.25,0.25], alpha=0.5,zorder=9)
        #plt.plot(LhErr,[row.R,row.R],'-', lw=1, color=[0.25,0.25,0.25], alpha=0.5,zorder=9)
        a0.plot(row.Lh,row.R,'.',ms=2,color=[0.25,0.25,0.25],zorder=32)

    # set up range of S_c values to test and empty array for results
    Sc_test = np.arange(0.5,1.,0.01)
    Percent_Less_Than_Segs = np.zeros(len(Sc_test))
    Percent_Less_Than_All = np.zeros(len(Sc_test))

    NoSegments = len(Data)
    NoHillslopes = len(HillslopeData)

    for i in range(0,len(Sc_test)):
        #print(Sc_test[i])

        # get max hillslpoe relief
        R_test = Data.Lh*Sc_test[i]

        #compare to actual hillslope relief and count how many fall below max line
        #first for segments
        BelowMask = Data.R < R_test
        NumberBelow = np.sum(BelowMask)
        Percent_Less_Than_Segs[i] = 100.*float(NumberBelow)/float(NoSegments)

        # get max hillslpoe relief
        R_test = HillslopeData.Lh*Sc_test[i]

        #and for all data
        BelowMask = HillslopeData.R < R_test
        NumberBelow = np.sum(BelowMask)
        Percent_Less_Than_All[i] = 100.*float(NumberBelow)/float(NoHillslopes)

    ind = np.argmin(np.abs(99.-Percent_Less_Than_Segs))
    print("Sc = "+str(Sc_test[ind]))
    Sc = np.around(Sc_test[ind],decimals=2)
    R = LH*Sc

    a0.plot(LH,R,'--', color='r',zorder=36)
    a0.text(10,50,"$S_C = $"+str(Sc),rotation=20.)

    plt.xlabel("Hillslope Length (m)")
    plt.ylabel("Hillslope Relief (m)")

    a0.set_xlim(0, 200)
    a0.set_ylim(0, 200)

    a0.text(0.05*2., 0.95, "(a)", transform=a0.transAxes, va='top', ha='right')

    # Add axes for plotting Sc vs % less than
    a1 = plt.subplot(gs[1])

    a1.plot(Sc_test,Percent_Less_Than_Segs,'k-')
    a1.plot(Sc_test,Percent_Less_Than_All,'-',color=[0.5,0.5,0.5])
    a1.plot([0.5,Sc_test[ind],Sc_test[ind]],[99.,99.,0.],'r--',lw=0.5)
    a1.text(Sc_test[ind],85,"$S_C = $"+str(Sc),rotation=-90)
    plt.xlabel("$S_C$")
    plt.ylabel("Percent Lower Relief")
    plt.ylim(70,100)
    a1.set_xlim(0.5,1.0)
    a1.yaxis.tick_right()
    a1.yaxis.set_label_position("right")

    a1.text(0.95, 0.95, "(b)", transform=a1.transAxes, va='top', ha='right')

    plt.tight_layout()
    plt.savefig(PlotDirectory+"Determine_Sc.png",dpi=300)
    plt.savefig(PlotDirectory+"Determine_Sc.pdf")

def CalculateEStarRStar(DataDirectory,FilenamePrefix,Basin,Sc=0.71):

    """
    Calculate EStar and RStar here so that you can change the critical slope
    Calculate for a specific basin.

    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): the file name prefix
        Sc (float): The critical slope to use

    returns: pandas data frame with Estar Rstar data and quantiles for hillslopes
        organised by channel segments for the specified basin

    MDH, Septmeber 2017

    """

    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)

    #load the hillslopes data
    HillslopeData = ReadHillslopeData(DataDirectory, FilenamePrefix)

    # isolate basin data
    BasinChannelData = ChannelData[ChannelData.basin_key == Basin]

    # segments in the hillslope data
    #Segments = BasinHillslopeData.StreamID.unique()
    Segments = BasinChannelData.segment_number.unique()

    # For each segment get the MChi value and collect dimensionless hillslope data
    # record the number of traces in each segment inbto a new dataframe
    Data = pd.DataFrame(columns=['SegmentNo','MChi','FlowLength','SegmentLength','EStar','EStarLower','EStarUpper','RStar','RStarLower','RStarUpper','NTraces'])

    # Loop through the segments
    for i in range(0,len(Segments)):

        #Get segment hillslope data
        SegmentHillslopeData = HillslopeData[HillslopeData.StreamID == float(Segments[i])]

        #Get segment channel data and calculate flow length
        SegmentChannelData = BasinChannelData[BasinChannelData.segment_number == Segments[i]]

        #channels
        MChi = SegmentChannelData.m_chi.unique()[0]
        TempFL = SegmentChannelData.flow_distance
        FlowLength = np.median(TempFL)
        SegmentLength = np.max(TempFL)-np.min(TempFL)

        #hillslopes
        TempEs = (-2.*SegmentHillslopeData.Cht*SegmentHillslopeData.Lh)/Sc
        TempRs = SegmentHillslopeData.S/Sc


        #get the stats to plot
        EStar = TempEs.quantile(0.5)
        EStarUpper = TempEs.quantile(0.75)
        EStarLower = TempEs.quantile(0.25)
        RStar = TempRs.quantile(0.5)
        RStarUpper = TempRs.quantile(0.75)
        RStarLower = TempRs.quantile(0.25)

        NTraces = SegmentHillslopeData.size

        #add to data frame
        Data.loc[i] = [Segments[i],MChi,FlowLength,SegmentLength,EStar,EStarLower,EStarUpper,RStar,RStarLower,RStarUpper,NTraces]

    # remove rows with no data (i.e. no hillslope traces)
    Data = Data.dropna(0,'any')

    # only keep segments with more than 50 hillslope traces
    #Data = Data[Data.NTraces > 50]

    return Data


def CalculateRStar(EStar):
    """
    Calculates E*

    returns: an array (?) of RStar values.

    MDH

    """

    RStar = (1./EStar)*(np.sqrt(1.+(EStar**2.)) - np.log(0.5*(1.+np.sqrt(1+EStar**2.))) - 1.)
    return RStar

def PlotEStarRStarTheoretical():
    """
    This makes the theoretical E* vs R* plot. It prints to the current open figure.
    SMM Note: This would be better if it used a supploed figure. Can the default be get_clf()?

    MDH

    """
    # Calculate analytical relationship
    EStar = np.logspace(-1,3,1000)
    RStar = CalculateRStar(EStar)

    # Plot with open figure
    plt.plot(EStar,RStar,'k--')


#-------------------------------------------------------------------------------#
# PLOTTING FUNCTIONS
#-------------------------------------------------------------------------------#
# SMM: Checked and working 13/06/2018
def PlotChiElevationSegments(DataDirectory, FilenamePrefix, PlotDirectory, BasinID):
    """
    This plots the chi--elevation prfile with the segments used in the hilltop analyses plotted in random colours.
    The segments are not the same as the ones determined by the segmentation algorithm. Instead they are bits of the chi
    profile split up to see the correspondence between channel and hillslope data.

    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): the file name prefix
        PlotDirectory (str): The directory into which the plots are saved
        BasinID (int): The basin to be plotted

    Author: MDH

    """

    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)

    # isolate basin data
    BasinChannelData = ChannelData[ChannelData.basin_key == BasinID]
    MinimumChi = BasinChannelData.chi.min()

    # how many segments are we dealing with?
    Segments = BasinChannelData.segment_number.unique()

    # setup the figure
    Fig = CreateFigure()

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
    plt.close(Fig)

# SMM: Checked and working 13/06/2018
def PlotLongProfileSegments(DataDirectory, FilenamePrefix, PlotDirectory, BasinID):
    """
    This plots the chi--elevation prfile with the segments used in the hilltop analyses plotted in random colours.
    The segments are not the same as the ones determined by the segmentation algorithm. Instead they are bits of the chi
    profile split up to see the correspondence between channel and hillslope data.

    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): the file name prefix
        PlotDirectory (str): The directory into which the plots are saved
        BasinID (int): The basin to be plotted

    Author: MDH


    """


    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)

    # isolate basin data
    BasinChannelData = ChannelData[ChannelData.basin_key == BasinID]
    MinimumDistance = BasinChannelData.flow_distance.min()

    # how many segments are we dealing with?
    Segments = BasinChannelData.segment_number.unique()

    # setup the figure
    Fig = CreateFigure()

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
    plt.close(Fig)

# SMM: Checked and working 13/06/2018
def PlotChiElevationMChi(DataDirectory, FilenamePrefix, PlotDirectory, BasinID):
    """
    This function reads the channel data file and plots the chi-elevation profile along with the segments extracted from the segmentation algorithm.
    It also colours the plot with the M_chi value (or k_sn if A_0 = 1).

    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): the file name prefix
        PlotDirectory (str): The directory into which the plots are saved
        BasinID (int): The basin to be plotted

    Author: MDH
    """

    print("Plotting the chi-elevation plot for basin: " +str(BasinID))
    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)

    # isolate basin data
    BasinChannelData = ChannelData[ChannelData.basin_key == BasinID]
    MinimumChi = BasinChannelData.chi.min()
    MaximumMChi = BasinChannelData.m_chi.max()

    # how many segments are we dealing with?
    Segments = BasinChannelData.segment_number.unique()

    # separate into main stem and trib data
    # SMM this isn't used
    #MainStemChannelData = BasinChannelData[BasinChannelData.source_key == 0]
    #MainStemSegments = MainStemChannelData.segment_number.unique()


    # Get the minimum and maximum chi
    MinimumChi = BasinChannelData.chi.min()
    MaximumMChi = BasinChannelData.m_chi.max()

    # how many segments are we dealing with?
    Segments = BasinChannelData.segment_number.unique()


    # setup the figure
    Fig = CreateFigure(AspectRatio=4./3.)

    #choose colormap
    ColourMap = cm.viridis

    # Get the data columns for plotting
    for i in range(0, len(Segments)):
        #if Segments[i] in MainStemSegments:
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
    plt.close(Fig)

# SMM: Checked and working 13/06/2018
def PlotLongProfileMChi(DataDirectory, FilenamePrefix, PlotDirectory, BasinID):
    """
    This function reads the channel data file and plots the long profile along with the segments extracted from the segmentation algorithm.
    It also colours the plot with the M_chi value (or k_sn if A_0 = 1).

    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): the file name prefix
        PlotDirectory (str): The directory into which the plots are saved
        BasinID (int): The basin to be plotted

    Author: MDH


    """

    print("Plotting the distance-elevation plot for basin: " +str(BasinID))
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
    plt.close(Fig)

# SMM: Checked and working 13/06/2018
# However it has a number of things that need to be revised: we should calculate E* R* directly (I think)
# Also error measurements need to be in quartiles rather than means and standard deviations
def PlotHillslopeDataVsDistance(DataDirectory, FilenamePrefix, PlotDirectory, BasinID, plot_vs_chi = False, minimum_traces = 50):
    """
    This function makes some composite plots of the hillslope data vs
    distance upstream from the outlet.
    It requires the hillslope and channel data csv files.
    Note: It reads E* R* data from file so S_c is not tuned. It must be entered in the c++ code.
    Note2: Uses standard deviation for errors so should be used with extreme caution

    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): the file name prefix
        PlotDirectory (str): The directory into which the plots are saved
        BasinID (int): The basin to be plotted

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

    # Get the minimum and maximum chi
    MinimumChi = BasinChannelData.chi.min()
    MaximumMChi = BasinChannelData.m_chi.max()

    # how many segments are we dealing with?
    Segments = BasinChannelData.segment_number.unique()

    # separate into main stem and trib data
    # (SMM): This is hard coded to take source key == 0 but this is only the source key for the first basin

    # try to figure out the source key
    mainstem_source_key = BasinChannelData.source_key.iloc[0]
    print("The mainstem source key is: " +str(mainstem_source_key))

    MainStemChannelData = BasinChannelData[BasinChannelData.source_key == mainstem_source_key]
    MainStemSegments = MainStemChannelData.segment_number.unique()

    print("The main stem segment numbers are: ")
    print(MainStemSegments)

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


    # COMMENT (SMM): Why does this not use chi coordinate? I will try to add it
    chi_coord = []

    for i in range (0, len(MainStemSegments)):

        # Isolate the correct hillslope and channel segments
        SegmentHillslopeData = BasinHillslopeData[BasinHillslopeData.StreamID == MainStemSegments[i]]
        SegmentChannelData = BasinChannelData[BasinChannelData.segment_number == MainStemSegments[i]]
        N_traces = len(SegmentHillslopeData["i"].tolist())
        #print("Sid: "+str(MainStemSegments[i])+" Bjunc: "+str(BasinJunctions[BasinID])+" Number of hilltops are: "+ str(len(SegmentHillslopeData["i"].tolist())))


        if N_traces > minimum_traces: # remove traces less than the minimum number of traces
            chi_coord.append(SegmentChannelData.chi.median())
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
        #else:
        #    print("Sid: "+str(MainStemSegments[i])+" Bjunc: "+str(BasinJunctions[BasinID])+" # Traces = "+str(N_traces)+ ", not enough at chi = "+ str(SegmentChannelData.chi.median()-MinimumChi))


    # normalise by the outlet chi
    chi_coord = chi_coord-MinimumChi

    # set up the figure
    fig, ax = plt.subplots(nrows = 4, ncols=1, sharex=True, figsize=(6,7))
    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=0)


    if(plot_vs_chi):
        print(chi_coord)
        ax[0].errorbar(chi_coord,Lh,yerr=Lh_std,fmt='o', ecolor='0.5',markersize=5,mec='k')

        #plot the cht
        ax[1].errorbar(chi_coord,Cht,yerr=Cht_std,fmt='o', ecolor='0.5',markersize=5,mfc='red',mec='k')

        #plot the R*
        ax[2].errorbar(chi_coord,R_star,yerr=R_star_std,fmt='o', ecolor='0.5',markersize=5,mfc='orange',mec='k')

        #plot the M_chi
        ax[3].errorbar(chi_coord,M_chi,yerr=M_chi_std,fmt='o', ecolor='0.5',markersize=5,mfc='purple',mec='k')

        # set the axes labels
        ax[3].set_xlabel('Chi (m)')
    else:
        print(DistanceFromOutlet)

        # plot the hillslope length
        ax[0].errorbar(DistanceFromOutlet,Lh,yerr=Lh_std,fmt='o', ecolor='0.5',markersize=5,mec='k')

        #plot the cht
        ax[1].errorbar(DistanceFromOutlet,Cht,yerr=Cht_std,fmt='o', ecolor='0.5',markersize=5,mfc='red',mec='k')

        #plot the R*
        ax[2].errorbar(DistanceFromOutlet,R_star,yerr=R_star_std,fmt='o', ecolor='0.5',markersize=5,mfc='orange',mec='k')

        #plot the M_chi
        ax[3].errorbar(DistanceFromOutlet,M_chi,yerr=M_chi_std,fmt='o', ecolor='0.5',markersize=5,mfc='purple',mec='k')

        # set the axes labels
        ax[3].set_xlabel('Distance from outlet (km)')

    # remaining axis labels
    ax[0].set_ylabel('$L_h$')
    ax[1].set_ylabel('$C_{HT}$')
    ax[2].set_ylabel('$R*$')
    ax[3].set_ylabel('$k_{sn}$')


    plt.tight_layout()

    #save output
    plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(BasinID) + "_hillslopes_distance.png", dpi=300)
    # close this figure to prevent stupid warnings
    plt.close(fig)

# SMM: Checked and working 13/06/2018
# I've modified this so I think it gives us all the stuff we need, using the correct statistics (medians and quartiles)
def PlotEStarRStarWithinBasin(DataDirectory, FilenamePrefix, PlotDirectory, BasinID, minimum_traces = 50, Sc = 0.8, plot_mainstem_only = False, colour_by = "default"):
    """
    Makes a plot of E* against R* where the points are coloured by
    their distance from the outlet of the basin or k_sn or chi

    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): the file name prefix
        PlotDirectory (str): The directory into which the plots are saved
        BasinID (int): The basin to be plotted
        minimum_traces (int): the minimum number of traces before a data point is recorded
        Sc (float): the critical slope to be used in E* R* calculations
        plot_mainstem_onlt (bool): If true only plot the mainstem data
        colour_by (str): What the data points should be coloured by. Options are chi and distance, and anything else will be coloured by k_sn
        calculate_Es_Rs (bool): If false, reads E* R* from the hillslope file. If true calculates it directly

    Author: FJC, SMM
    """
    import math

    print("Plotting the E* R* curves for basin "+str(BasinID))
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

    # try to figure out the source key
    mainstem_source_key = BasinChannelData.source_key.iloc[0]
    print("The mainstem source key is: " +str(mainstem_source_key))

    # separate into main stem and trib data
    MainStemChannelData = BasinChannelData[BasinChannelData.source_key == mainstem_source_key]
    MainStemSegments = MainStemChannelData.segment_number.unique()

    # set up the figure
    Fig = CreateFigure()
    Ax = plt.subplot(111)

    #choose colormap
    ColourMap = cm.viridis


    CalcEsMed = []
    CalcEsLQ = []
    CalcEsUQ = []

    CalcRsMed = []
    CalcRsLQ = []
    CalcRsUQ = []

    AllDist = []
    AllChi = []
    AllKsn = []

    EsLowErr = []
    EsUpErr = []
    RsLowErr = []
    RsUpErr = []

    if plot_mainstem_only:
        print("I am only going to plot the main stem data.")

    total_data_points = 0

    for i in range (0, len(Segments)):
        SegmentHillslopeData = BasinHillslopeData[BasinHillslopeData.StreamID == Segments[i]]
        SegmentChannelData = BasinChannelData[BasinChannelData.segment_number == Segments[i]]

        # Get the number of traces
        N_traces = len(SegmentHillslopeData["i"].tolist())

        if N_traces > minimum_traces:

            # keep track of how many data points you have
            total_data_points = total_data_points+1

            # Calculate E* R*
            if plot_mainstem_only:
                if Segments[i] in MainStemSegments:
                    EStar = -2*SegmentHillslopeData.Cht*SegmentHillslopeData.Lh/Sc
                    RStar = SegmentHillslopeData.S/Sc

                    # Get distances and chi values
                    AllDist.append(SegmentChannelData.flow_distance.median()/1000)
                    AllChi.append(SegmentChannelData.chi.median())

                    # Now get the k_sn
                    AllKsn.append(SegmentChannelData.m_chi.unique()[0])

                    # Get the medians and the quartiles.
                    CalcEsMed.append(EStar.median())
                    CalcEsLQ.append(EStar.quantile(0.25))
                    CalcEsUQ.append(EStar.quantile(0.75))

                    EsLowErr.append(EStar.median()-EStar.quantile(0.25))
                    EsUpErr.append(EStar.quantile(0.75)-EStar.median())

                    CalcRsMed.append(RStar.median())
                    CalcRsLQ.append(RStar.quantile(0.25))
                    CalcRsUQ.append(RStar.quantile(0.75))

                    RsLowErr.append(RStar.median()-RStar.quantile(0.25))
                    RsUpErr.append(RStar.quantile(0.75)-RStar.median())

            else:
                EStar = -2*SegmentHillslopeData.Cht*SegmentHillslopeData.Lh/Sc
                RStar = SegmentHillslopeData.S/Sc

                # Get distances and chi values
                AllDist.append(SegmentChannelData.flow_distance.median()/1000)
                AllChi.append(SegmentChannelData.chi.median())

                # Now get the k_sn
                AllKsn.append(SegmentChannelData.m_chi.unique()[0])
                #print("Segment m_chi is: "+ str(SegmentChannelData.m_chi.unique()[0]))

                # Get the medians and the quartiles.
                CalcEsMed.append(EStar.median())
                CalcEsLQ.append(EStar.quantile(0.25))
                CalcEsUQ.append(EStar.quantile(0.75))

                EsLowErr.append(EStar.median()-EStar.quantile(0.25))
                EsUpErr.append(EStar.quantile(0.75)-EStar.median())

                CalcRsMed.append(RStar.median())
                CalcRsLQ.append(RStar.quantile(0.25))
                CalcRsUQ.append(RStar.quantile(0.75))

                RsLowErr.append(RStar.median()-RStar.quantile(0.25))
                RsUpErr.append(RStar.quantile(0.75)-RStar.median())

                #print("Rs: "+str(RStar.quantile(0.25))+","+str(RStar.median())+","+str(RStar.quantile(0.75)))


    plt.loglog()
    PlotEStarRStarTheoretical()

    #print(len(CalcEsMed))
    #print(len(EsLowErr))
    #print(len(EsUpErr))

    Ax.errorbar(CalcEsMed,CalcRsMed,xerr=[EsLowErr, EsUpErr], yerr=[RsLowErr,RsUpErr],fmt='.', ecolor='k',markersize=2,mec='k',mfc='k',zorder = 10, linewidth = 1, alpha = 0.5)

    if colour_by == "chi":
        Ax.scatter(CalcEsMed,CalcRsMed,c=AllChi,s=10, edgecolors='k', lw=0.1,cmap=ColourMap,zorder = 20)
    elif colour_by == "distance":
        Ax.scatter(CalcEsMed,CalcRsMed,c=AllDist,s=10, edgecolors='k', lw=0.1,cmap=ColourMap,zorder = 20)
    else:
        Ax.scatter(CalcEsMed,CalcRsMed,c=AllKsn,s=10, edgecolors='k', lw=0.1,cmap=ColourMap,zorder = 20)

    # Ax.set_yscale('log')
    plt.xlabel('$E*$')
    plt.ylabel('$R*$')


    # Make the plot
    plt.subplots_adjust(left=0.18,right=0.85, bottom=0.2, top=0.9)
    CAx = Fig.add_axes([0.86,0.2,0.02,0.7])
    m = cm.ScalarMappable(cmap=ColourMap)

    if colour_by == "chi":
        m.set_array(AllChi)
        plt.colorbar(m, cax=CAx,orientation='vertical', label='$\chi$ (m)')
        figappendstr = "_EStar_RStar_chi.png"
    elif colour_by == "distance":
        m.set_array(AllDist)
        plt.colorbar(m, cax=CAx,orientation='vertical', label='Distance from outlet (km)')
        figappendstr = "_EStar_RStar_dist.png"
    else:
        m.set_array(AllKsn)
        plt.colorbar(m, cax=CAx,orientation='vertical', label='$k_{sn}$')
        figappendstr = "_EStar_RStar_ksn.png"

    plt.text(0.6, 0.025, 'Basin: '+str(BasinID), horizontalalignment='left', verticalalignment='bottom', transform=Ax.transAxes, fontsize = 12)

    print("total_data_points: "+ str(total_data_points))

    if total_data_points==0:
        plt.text(0.05, 0.95, 'Basin'+str(BasinID)+", no data points.", horizontalalignment='left', verticalalignment='top', transform=Ax.transAxes)

    plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(BasinID).zfill(3) + figappendstr, dpi=300)
    plt.clf()
    plt.close(Fig)


# This is only functional for Mendocino so not general
# I don't think it would take too much effort, however, to look for the uplift file and just not plot uplift
# if the file is missing. However that is a task for another day (SMM, 13/06/2018)
def PlotHillslopeDataWithBasins(DataDirectory,FilenamePrefix,PlotDirectory):
    """
    Function to make plots of hillslope data vs basin id.
    At the moment this is hard coded for the MTJ because I need to add in the
    uplift data...sorry!

    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): the file name prefix
        PlotDirectory (str): The directory into which the plots are saved
        BasinID (int): The basin to be plotted

    Author: FJC
    """

    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)
    #print BasinChannelData

    # load the hillslopes data
    HillslopeData = ReadHillslopeData(DataDirectory, FilenamePrefix)

    basin_dict = MapBasinKeysToJunctions(DataDirectory,FilenamePrefix)
    basin_keys = list(basin_dict.keys())

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

    for key, jn in basin_dict.items():
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
    fig, ax = plt.subplots(nrows = 7, ncols=1, sharex=True, figsize=(6,12), facecolor='white')
    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=0)

    # get the colours
    cmap = cm.Dark2
    colors = LSDP.colours.list_of_hex_colours(7, cmap)

    # plot the hillslope length
    ax[0].errorbar(basin_keys,median_Lh,yerr=[Lh_lower_err, Lh_upper_err],fmt='o', ecolor='0.5',markersize=6,mec='k',mfc=colors[0])
    ax[0].set_ylabel('$L_H$')

    #plot the cht
    ax[1].errorbar(basin_keys,median_cht,yerr=[cht_lower_err, cht_upper_err],fmt='o', ecolor='0.5',markersize=6,mfc=colors[1],mec='k')
    ax[1].set_ylabel('$C_{HT}$')

    #plot the E*
    ax[2].errorbar(basin_keys,median_Estar,yerr=[Estar_lower_err, Estar_upper_err],fmt='o', ecolor='0.5',markersize=6,mfc=colors[2],mec='k')
    ax[2].set_ylabel('$E*$')

    #plot the R*
    ax[3].errorbar(basin_keys,median_Rstar,yerr=[Rstar_lower_err, Rstar_upper_err],fmt='o', ecolor='0.5',markersize=6,mfc=colors[3],mec='k')
    ax[3].set_ylabel('$R*$')

    #plot the Mchi
    ax[4].errorbar(basin_keys,median_mchi,yerr=[mchi_lower_err, mchi_upper_err],fmt='o', ecolor='0.5',markersize=6,mfc=colors[5],mec='k')
    ax[4].set_ylabel('$k_{sn}$')

    # read the uplift data in
    # read in the csv
    uplift_df = pd.read_csv(DataDirectory+'MTJ_basin_uplift.csv')
    dd_df = pd.read_csv(DataDirectory+FilenamePrefix+'_basin_dd.csv')

    # get the drainage density
    drainage_density = dd_df['drainage_density']*1000000
    ax[5].scatter(basin_keys, drainage_density, c=colors[6], edgecolors='k', s=30)
    ax[5].set_ylim(np.min(drainage_density)-1000, np.max(drainage_density)+1000)
    ax[5].set_ylabel('$D_d$ (m/km$^2$)')

    # get the data
    uplift_rate_old = uplift_df['Uplift_rate_old']
    uplift_rate_new = uplift_df['Uplift_rate_new']
    #ax[6].scatter(basin_keys, uplift_rate_new, c='None', edgecolors='k', label = '0 - 72 ka')
    ax[6].plot(basin_keys, uplift_rate_new, c='k', ls='--', label = '96 - 305 ka')
    ax[6].set_ylabel('Uplift rate (mm/yr)')
    #ax[6].legend(loc='upper right')

    # erosion rate
    be_erosion = uplift_df['Erosion_rate_Be']
    al_erosion = uplift_df['Erosion_rate_Al']
    al_min_erosion = uplift_df['Al_min']
    # # ax[7].scatter(basin_keys, be_erosion, c='k', label='Be')
    # ax[7].errorbar(basin_keys, be_erosion, xerr=None, yerr=uplift_df['Be_error'], fmt='o', ecolor='0.5',markersize=6,mec='k', mfc='k', label='Be')
    # # ax[7].scatter(basin_keys, al_erosion, c='None', edgecolors='k', marker='D', label='Al')
    # ax[7].errorbar(basin_keys, al_erosion, xerr=None, yerr=uplift_df['Al_error'], fmt='D', ecolor='0.5', markersize=6, mec='k', mfc='white', label='Al')
    # ax[7].scatter(basin_keys, al_min_erosion, c='None', edgecolors='k', marker='^', label='Al (min)')
    # ax[7].set_ylabel('Erosion rate (mm/yr)')
    # ax[7].legend(loc='upper right')

    # set the axes labels
    ax[6].set_xlabel('Basin ID')
    plt.xticks(np.arange(min(basin_keys), max(basin_keys)+1, 1), rotation=45, fontsize=8)
    plt.tight_layout()
    #plt.subplots_adjust(bottom=0.1)

    #save output
    plt.savefig(PlotDirectory+FilenamePrefix +"_basin_hillslope_data.pdf", dpi=300)
    plt.clf()

    output_list = [('basin_keys', basin_keys),
                   ('uplift_rate_old', uplift_rate_old),
                   ('uplift_rate_new', uplift_rate_new),
                   ('Erosion_rate_Be', be_erosion),
                   ('Be_error', uplift_df['Be_error']),
                   ('Erosion_rate_Al', al_erosion),
                   ('Al_error', uplift_df['Al_error']),
                   ('Erosion_rate_Al_min', al_min_erosion),
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

# SMM: Checked 13/06/2018 but not working since I do not know there the "_basin_hillslope_data.csv" comes from and don't have it.
def PlotKsnAgainstRStar(DataDirectory, FilenamePrefix, PlotDirectory):
    """
    Function to plot median Ksn against R* for a series of basins

    Author: FJC
    """

    # SMM: What generates this file?? I don't have it.
    input_csv = PlotDirectory+FilenamePrefix+'_basin_hillslope_data.csv'
    df = pd.read_csv(input_csv)

    # linregress
    slope, intercept, r_value, p_value, std_err = stats.linregress(df['mchi_median'],df['Rstar_median'])
    print(slope, intercept, r_value, p_value)
    x = np.linspace(0, 200, 100)
    new_y = slope*x + intercept

    # set up the figure
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=(5,5))

    ax.scatter(df['mchi_median'], df['Rstar_median'], c=df['basin_keys'], s=50, edgecolors='k', zorder=100, cmap=cm.viridis)
    ax.errorbar(df['mchi_median'], df['Rstar_median'], xerr=[df['mchi_lower_err'], df['mchi_upper_err']], yerr=[df['Rstar_lower_err'], df['Rstar_upper_err']], fmt='o', ecolor='0.5',markersize=1,mfc='white',mec='k')
    # ax.text(0.55, 0.1, '$y = $'+str(np.round(slope,4))+'$x + $'+str(np.round(intercept,2))+'\n$R^2 = $'+str(np.round(r_value,2))+'\n$p = $'+str(p_value), fontsize=9, color='black', transform=ax.transAxes)
    ax.plot(x, new_y, c='0.5', ls='--')
    ax.set_xlim(0,100)

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


# This seems to do the same as the PlotEStarRStarWithinBasin function!!!
# However it is not working since I don't have the _basin_hillslope_data.csv' file
def PlotEStarRStarBasins(DataDirectory, FilenamePrefix, PlotDirectory, Sc = 0.8):
    """
    Function to make an E*R* plot for a series of drainage basins.
    Changing so that we calculate E* and R* in the python script following
    Martin's example, so that we can test sensitivity to Sc.

    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): the file name prefix
        PlotDirectory (str): The directory into which the plots are saved
        Sc (float): The critical slope to be used in the analysis

    Author: FJC
    """

    # SMM: It is not clear where this file comes from
    input_csv = PlotDirectory+FilenamePrefix+'_basin_hillslope_data.csv'
    df = pd.read_csv(input_csv)

    # set up the figure
    fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=(5,5))
    PlotEStarRStarTheoretical()

    #choose colormap
    ColourMap = cm.viridis

    # get the basins
    basins = df['basin_keys'].unique()
    NoBasins = len(basins)
    MinMChi = df.mchi_median.min()
    MaxMChi = df.mchi_median.max()

    for basin_key in basins:
        Data = CalculateEStarRStar(DataDirectory,FilenamePrefix,basin_key, Sc=Sc)

        # colour code by basin number
        #colour = float(basin_key)/float(NoBasins)
        colour = df.mchi_median[df.basin_keys == basin_key].values[0]
        EStarMedian = Data.EStar.median()
        RStarMedian = Data.RStar.median()
        EStar_lower_err = np.percentile(Data.EStar.as_matrix(), 25)
        EStar_upper_err = np.percentile(Data.EStar.as_matrix(), 75)
        RStar_lower_err = np.percentile(Data.RStar.as_matrix(), 25)
        RStar_upper_err = np.percentile(Data.RStar.as_matrix(), 75)

        cNorm  = colors.Normalize(vmin=MinMChi, vmax=MaxMChi)
        plt.cm.ScalarMappable(norm=cNorm, cmap=ColourMap)

        # plot the rstar vs estar
        sc = ax.scatter(EStarMedian,RStarMedian,c=colour,s=50, edgecolors='k', zorder=100, norm=cNorm)
        ax.errorbar(EStarMedian,RStarMedian,xerr=[[EStarMedian-EStar_lower_err],[EStar_upper_err-EStarMedian]], yerr=[[RStarMedian-RStar_lower_err],[RStar_upper_err-RStarMedian]],fmt='o', zorder=1, ecolor='0.5',markersize=1,mfc='white',mec='k')

    # Finalise the figure
    plt.xlabel('$E^*={{-2\:C_{HT}\:L_H}/{S_C}}$')
    plt.ylabel('$R^*=S/S_C$')
    plt.xlim(0.1,20)
    plt.ylim(0.05,1)

    # add colour bar
    cbar = plt.colorbar(sc,cmap=ColourMap)
    colorbarlabel='Basin ID'
    cbar.set_label(colorbarlabel, fontsize=10)
    # tick_locator = ticker.MaxNLocator(nbins=5)
    # cbar.locator = tick_locator
    # cbar.update_ticks()

    #save output
    plt.savefig(PlotDirectory+FilenamePrefix +"_estar_vs_rstar{}.png".format(Sc), dpi=300)
    plt.clf()

def PlotEStarRStarSubPlots(DataDirectory, FilenamePrefix, PlotDirectory, Sc = 0.8):
    """
    Make a composite plot of E* R* and R* Ksn

    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): the file name prefix
        PlotDirectory (str): The directory into which the plots are saved
        Sc (float): The critical slope to be used in the analysis

    FJC
    """
    input_csv = PlotDirectory+FilenamePrefix+'_basin_hillslope_data.csv'
    df = pd.read_csv(input_csv)

    # set up the figure
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10,5))

    #choose colormap
    ColourMap = cm.RdYlBu

    # get the basins
    basins = df['basin_keys'].unique()
    NoBasins = len(basins)
    print(basins)

    sc = ax[0].scatter(df.Estar_median,df.Rstar_median,c=basins,s=50, edgecolors='k', zorder=100, cmap=ColourMap)
    ax[0].errorbar(df.Estar_median,df.Rstar_median,xerr=[df['Estar_lower_err'], df['Estar_upper_err']], yerr=[df['Rstar_lower_err'], df['Rstar_upper_err']],fmt='o', zorder=1, ecolor='0.5',markersize=1,mfc='white',mec='k')

    sc = ax[1].scatter(df.mchi_median,df.Rstar_median,c=basins,s=50, edgecolors='k', zorder=100, cmap=ColourMap)
    ax[1].errorbar(df.mchi_median,df.Rstar_median,xerr=[df['mchi_lower_err'], df['mchi_upper_err']], yerr=[df['Rstar_lower_err'], df['Rstar_upper_err']],fmt='o', zorder=1, ecolor='0.5',markersize=1,mfc='white',mec='k')

    # plot the theoretical relationships for each one
    # Calculate analytical relationship for estar rstar
    EStar_model = np.logspace(-1,3,1000)
    RStar_model = CalculateRStar(EStar_model)

    # Plot with open figure
    ax[0].plot(EStar_model,RStar_model,c='0.5',ls='--')

    # calculate linear fit for Rstar ksn
    slope, intercept, r_value, p_value, std_err = stats.linregress(df.mchi_median, df.Rstar_median)
    print(slope, intercept, r_value, p_value)
    x = df.mchi_median.values
    print(x)
    new_y = slope*x + intercept
    ax[1].plot(x, new_y, c='0.5', ls='--')

    # get the difference between the linear fit and the real R* for each basin and
    # print to csv for plotting
    residual = df.Rstar_median.values - new_y
    print(residual)
    df['rstar_ksn_residual'] = residual
    OutputFilename = PlotDirectory+FilenamePrefix+'_basin_hillslope_data_residuals.csv'
    df.to_csv(OutputFilename, index=False)

    # Finalise the figure
    ax[0].set_xlabel('$E^*={{-2\:C_{HT}\:L_H}/{S_C}}$')
    ax[0].set_ylabel('$R^*=S/S_C$')
    ax[0].set_xlim(0.1,25)
    ax[0].set_ylim(0.2,1)

    ax[1].set_xlim(10,90)
    ax[1].set_ylim(0.2,1)
    ax[1].set_xlabel('$k_{sn}$')

    # add colour bar
    m = cm.ScalarMappable(cmap=ColourMap)
    m.set_array(basins)
    cbar = plt.colorbar(m, ax=ax.ravel().tolist())
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.ax.invert_yaxis()
    cbar.set_label('Basin ID')

    print("Made the E*R* plots")

    #save output
    plt.savefig(PlotDirectory+FilenamePrefix +"_estar_rstar_subplots.png", dpi=300)
    plt.clf()

def PlotDataAgainstErosionRate(DataDirectory, FilenamePrefix, PlotDirectory):
    """
    Make plots of the data against erosion rate. This only works if you have
    used the function PlotHillslopeDataWithBasins to generate the correct csv
    file first. I wrote it for the MTJ analysis

    FJC 11/10/18
    """
    df = pd.read_csv(PlotDirectory+FilenamePrefix+'_basin_hillslope_data.csv')

    # set up the figure
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(12,5))
    #ax = ax.ravel()
    # make a big subplot to allow sharing of axis labels
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

    #choose colormap
    ColourMap = cm.viridis

    # get the basins
    basins = df['basin_keys'].unique()
    NoBasins = len(basins)
    print(basins)

    norm = colors.Normalize(vmin=basins.min(), vmax=basins.max())

    #sc = ax[0].scatter(df.Erosion_rate_Be,df.Rstar_median,c=basins,s=50, edgecolors='k', zorder=100)
    sc = ax[0].scatter(df.Erosion_rate_Al,df.Rstar_median,c=basins,s=50, edgecolors='k', zorder=100, norm=norm)
    sc = ax[0].scatter(df.Erosion_rate_Al_min,df.Rstar_median,c=basins,s=50, edgecolors='k', zorder=100, norm=norm)
    ax[0].errorbar(df.Erosion_rate_Al,df.Rstar_median,xerr=df['Al_error'], yerr=[df['Rstar_lower_err'], df['Rstar_upper_err']],fmt='o', zorder=1, ecolor='0.5',markersize=1,mfc='white',mec='k')
    ax[0].errorbar(df.Erosion_rate_Al_min, df.Rstar_median, yerr=[df['Rstar_lower_err'], df['Rstar_upper_err']],fmt='o', zorder=1, ecolor='0.5',markersize=1,mfc='white',mec='k')
    #ax[0].errorbar(df.Estar_median,df.Rstar_median,xerr=[df['Estar_lower_err'], df['Estar_upper_err']], yerr=[df['Rstar_lower_err'], df['Rstar_upper_err']],fmt='o', zorder=1, ecolor='0.5',markersize=1,mfc='white',mec='k')

    sc = ax[1].scatter(df.Erosion_rate_Al,df.Estar_median,c=basins,s=50, edgecolors='k', zorder=100, norm=norm)
    sc = ax[1].scatter(df.Erosion_rate_Al_min,df.Estar_median,c=basins,s=50, edgecolors='k', zorder=100, norm=norm)
    ax[1].errorbar(df.Erosion_rate_Al,df.Estar_median,xerr=df['Al_error'], yerr=[df['Estar_lower_err'], df['Estar_upper_err']],fmt='o', zorder=1, ecolor='0.5',markersize=1,mfc='white',mec='k')
    ax[1].errorbar(df.Erosion_rate_Al_min, df.Estar_median, yerr=[df['Estar_lower_err'], df['Estar_upper_err']],fmt='o', zorder=1, ecolor='0.5',markersize=1,mfc='white',mec='k')

    new_df = df[np.isnan(df.Erosion_rate_Al) == False]
    newbasins = new_df.basin_keys.unique()
    print(new_df.Erosion_rate_Al)

    sc = ax[2].scatter(new_df.Estar_median,new_df.Rstar_median,c=newbasins,s=50, edgecolors='k', zorder=100, norm=norm)
    #sc = ax[2].scatter(df.Erosion_rate_Al_min,df.cht_median,c=basins,s=50, edgecolors='k', zorder=100, norm=norm)
    #ax[2].errorbar(df.Erosion_rate_Al,df.cht_median,xerr=df['Al_error'], yerr=[df['cht_lower_err'], df['cht_upper_err']],fmt='o', zorder=1, ecolor='0.5',markersize=1,mfc='white',mec='k')
    ax[2].errorbar(new_df.Estar_median, new_df.Rstar_median, xerr=[new_df['Estar_lower_err'], new_df['Estar_upper_err']], yerr=[new_df['Rstar_lower_err'], new_df['Rstar_upper_err']],fmt='o', zorder=1, ecolor='0.5',markersize=1,mfc='white',mec='k')

    new_df = df[np.isnan(df.Erosion_rate_Al_min) == False]
    newbasins = new_df.basin_keys.unique()

    sc = ax[2].scatter(new_df.Estar_median,new_df.Rstar_median,c=newbasins,s=50, edgecolors='k', zorder=100, norm=norm)
    #sc = ax[2].scatter(df.Erosion_rate_Al_min,df.cht_median,c=basins,s=50, edgecolors='k', zorder=100, norm=norm)
    #ax[2].errorbar(df.Erosion_rate_Al,df.cht_median,xerr=df['Al_error'], yerr=[df['cht_lower_err'], df['cht_upper_err']],fmt='o', zorder=1, ecolor='0.5',markersize=1,mfc='white',mec='k')
    ax[2].errorbar(new_df.Estar_median, new_df.Rstar_median, xerr=[new_df['Estar_lower_err'], new_df['Estar_upper_err']], yerr=[new_df['Rstar_lower_err'], new_df['Rstar_upper_err']],fmt='o', zorder=1, ecolor='0.5',markersize=1,mfc='white',mec='k')

    # Finalise the figure
    ax[0].set_xlabel('Erosion rate (mm/yr)')
    ax[0].set_ylabel('$R^*=S/S_C$')

    ax[1].set_xlabel('Erosion rate (mm/yr)')
    ax[1].set_ylabel('$E^*={{-2\:C_{HT}\:L_H}/{S_C}}$')

    ax[2].set_xlabel('$E^*={{-2\:C_{HT}\:L_H}/{S_C}}$')
    ax[2].set_ylabel('$R^*=S/S_C$')
    #ax[2].set_ylim(0.005,0.0225)

    # add colour bar
    m = cm.ScalarMappable(cmap=ColourMap)
    m.set_array(basins)
    cax = fig.add_axes([0.91,0.1,0.02,0.8])
    cbar = fig.colorbar(m, cax=cax)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.set_label('Basin ID')
    #plt.tight_layout()
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.3)

    #save output
    plt.savefig(PlotDirectory+FilenamePrefix +"_hs_data_erosion_rate.png", dpi=300)
    plt.clf()

def Make3DHillslopePlot(DataDirectory, FilenamePrefix, PlotDirectory):
    """
    Function to make a 3d plot of E, R* and E* for a series of basins
    """
    df = pd.read_csv(PlotDirectory+FilenamePrefix+'_basin_hillslope_data.csv')

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    #choose colormap
    ColourMap = cm.viridis

    # get the basins
    basins = df['basin_keys'].unique()
    NoBasins = len(basins)
    print(basins)

    norm = colors.Normalize(vmin=basins.min(), vmax=basins.max())

    # plot the errorbars. this is annoying in 3d
    for i in range(len(df.Erosion_rate_Al)):
        # x error
        ax.plot([df.Erosion_rate_Al[i] + df.Al_error[i], df.Erosion_rate_Al[i] - df.Al_error[i]], [df.Estar_median[i], df.Estar_median[i]], [df.Rstar_median[i], df.Rstar_median[i]], c='0.4', zorder=-1)
        # y error
        ax.plot([df.Erosion_rate_Al[i], df.Erosion_rate_Al[i]], [df.Estar_median[i]+df.Estar_upper_err[i], df.Estar_median[i]-df.Estar_lower_err[i]], [df.Rstar_median[i], df.Rstar_median[i]], c='0.4', zorder=-1)
        # z error
        ax.plot([df.Erosion_rate_Al[i], df.Erosion_rate_Al[i]], [df.Estar_median[i], df.Estar_median[i]], [df.Rstar_median[i]+df.Rstar_upper_err[i], df.Rstar_median[i]-df.Rstar_lower_err[i]], c='0.4', zorder=-1)

    for i in range(len(df.Erosion_rate_Al_min)):
        # y error
        ax.plot([df.Erosion_rate_Al_min[i], df.Erosion_rate_Al_min[i]], [df.Estar_median[i]+df.Estar_upper_err[i], df.Estar_median[i]-df.Estar_lower_err[i]], [df.Rstar_median[i], df.Rstar_median[i]], c='0.4', zorder=-1)
        # z error
        ax.plot([df.Erosion_rate_Al_min[i], df.Erosion_rate_Al_min[i]], [df.Estar_median[i], df.Estar_median[i]], [df.Rstar_median[i]+df.Rstar_upper_err[i], df.Rstar_median[i]-df.Rstar_lower_err[i]], c='0.4', zorder=-1)

    # plot the data
    ax.scatter(df.Erosion_rate_Al, df.Estar_median, df.Rstar_median, c=basins, alpha=1, edgecolors='k', s=50, zorder=1, norm=norm)
    ax.scatter(df.Erosion_rate_Al_min, df.Estar_median, df.Rstar_median, c=basins, alpha=1, edgecolors='k', s=50, zorder=1, norm=norm)

    yflat = np.full_like(df.Estar_median, max(ax.get_ylim()))
    zflat = np.full_like(df.Rstar_median, min(ax.get_zlim()))

    new_df = df[np.isnan(df.Erosion_rate_Al) == False]
    newbasins = new_df.basin_keys.unique()
    xflat = np.full_like(new_df.Erosion_rate_Al, min(ax.get_xlim()))
    ax.scatter(xflat, new_df.Estar_median, new_df.Rstar_median,c=newbasins, alpha=0.2, edgecolors='k', s=50, zorder=1, norm=norm)

    new_df = df[np.isnan(df.Erosion_rate_Al_min) == False]
    newbasins = new_df.basin_keys.unique()
    xflat = np.full_like(new_df.Erosion_rate_Al_min, min(ax.get_xlim()))
    ax.scatter(xflat, new_df.Estar_median, new_df.Rstar_median,c=newbasins, alpha=0.2, edgecolors='k', s=50, zorder=-2, norm=norm)
    #ax.scatter(x2flat, df.Estar_median, df.Rstar_median,c=basins, alpha=0.5, edgecolors='k', s=50, zorder=1, norm=norm)
    ax.scatter(df.Erosion_rate_Al, yflat, df.Rstar_median,c=basins, alpha=0.2, edgecolors='k', s=50, zorder=-2, norm=norm)
    ax.scatter(df.Erosion_rate_Al, df.Estar_median, zflat,c=basins, alpha=0.2, edgecolors='k', s=50, zorder=-2, norm=norm)
    ax.scatter(df.Erosion_rate_Al_min, yflat, df.Rstar_median,c=basins, alpha=0.2, edgecolors='k', s=50, zorder=-2, norm=norm)
    ax.scatter(df.Erosion_rate_Al_min, df.Estar_median, zflat,c=basins, alpha=0.2, edgecolors='k', s=50, zorder=-2, norm=norm)


   # ax.plot([fx[i]+xerror[i], fx[i]-xerror[i]], [fy[i], fy[i]], [fz[i], fz[i]], marker="_")
   #  ax.plot([fx[i], fx[i]], [fy[i]+yerror[i], fy[i]-yerror[i]], [fz[i], fz[i]], marker="_")
   #  ax.plot([fx[i], fx[i]], [fy[i], fy[i]], [fz[i]+zerror[i], fz[i]-zerror[i]], marker="_")


    # ax.errorbar(df.Erosion_rate_Al, df.Rstar_median, df.Estar_median, xerr=df['Al_error'], yerr=[df['Rstar_lower_err'], df['Rstar_upper_err']], zerr=[df['Estar_lower_err'], df['Estar_upper_err']], fmt='o', zorder=1, ecolor='0.5',markersize=1,mfc='white',mec='k')

    # add colour bar
    m = cm.ScalarMappable(cmap=ColourMap)
    m.set_array(basins)
    cax = fig.add_axes([0.85,0.15,0.02,0.65])
    cbar = fig.colorbar(m, cax=cax)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.set_label('Basin ID')

    ax.set_xlabel('Erosion rate (mm/yr)')
    ax.set_ylabel('$E*$')
    ax.set_zlabel('$R*$')
    zed = [tick.label.set_fontsize(8) for tick in ax.xaxis.get_major_ticks()]
    zed = [tick.label.set_fontsize(8) for tick in ax.yaxis.get_major_ticks()]
    zed = [tick.label.set_fontsize(8) for tick in ax.zaxis.get_major_ticks()]
    plt.subplots_adjust(left=0.05, right=0.8, bottom=0.1, top=0.9)

    # make the grid lines dashed
    ax.xaxis._axinfo["grid"]['linestyle'] = ":"
    ax.yaxis._axinfo["grid"]['linestyle'] = ":"
    ax.zaxis._axinfo["grid"]['linestyle'] = ":"

    #save output
    plt.savefig(PlotDirectory+FilenamePrefix +"_hs_data_3d.png", dpi=300)
    plt.clf()


def PlotHillslopeTraces(DataDirectory, FilenamePrefix, PlotDirectory, CustomExtent=[-9999],FigSizeFormat="epsl"):
    """
    Function to plot a hillshade image with hilltops, hillslope traces and the channel network superimposed.

    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): the file name prefix
        PlotDirectory (str): The directory into which the plots are saved
        CustomExtent (list): If [-9999] then just use extent of raster. Otherwise a four lement list with extents of the area you want to plot
        FigSizeFormat (str): The format of the figure you want. Try your favourite journal. It may or may not be there.

    MDH
    """


    # Save the figure
    ImageName = PlotDirectory+FilenamePrefix+"_traces.png"
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
    MF.add_point_data(HilltopPoints,alpha=1,zorder=100,unicolor=[0.8,0,0],manual_size=2)

    # add channel heads
    ChannelHeadsDF = pd.read_csv(DataDirectory+FilenamePrefix+"_Wsources.csv")
    ChannelHeadPoints = LSDP.LSDMap_PointData(ChannelHeadsDF, data_type = "pandas", PANDEX = True)
    MF.add_point_data(ChannelHeadPoints,zorder=100,unicolor="blue",manual_size=8)

    # add channels
    ChannelDF = Helper.ReadChiDataMapCSV(DataDirectory,FilenamePrefix)
    ChannelPoints = LSDP.LSDMap_PointData(ChannelDF, data_type = "pandas", PANDEX = True)
    MF.add_point_data(ChannelPoints,show_colourbar="False", scale_points=True, max_point_size = 2.5, min_point_size = 0.5, column_for_scaling='drainage_area',zorder=90)

    # add hillslope traces
    ThinningFactor=1
    #HillslopeTracesShp = DataDirectory+FilenamePrefix+"_hillslope_traces.shp"
    #if os.path.exists(HillslopeTracesShp) == False:
    WriteHillslopeTracesShp(DataDirectory,FilenamePrefix,ThinningFactor,CustomExtent)

    MF.add_line_data(DataDirectory+FilenamePrefix+"_hillslope_traces.shp",zorder=80,alpha=0.9,linewidth=0.8)

    #finalise and save figure
    MF.SetRCParams(label_size=8)
    MF.save_fig(fig_width_inches = FigWidth_Inches, FigFileName = ImageName, FigFormat="png", Fig_dpi = 300)

def PlotEStarRStarProgression(Sc=0.71):
    """
    Plots the progression of hillslopes along Bolinas in Estar Rstar space. Maybe Fiona just copied this over since I'm not sure from where it will read data (SMM)

    Args:
        Sc (float): The critical slope

    MDH, September 2017

    """

    from scipy.stats import gaussian_kde

    # setup the figure
    Fig = CreateFigure(AspectRatio=1.2)

    #choose colormap
    ColourMap = cm.viridis

    #Plot analytical relationship
    plt.loglog()
    PlotEStarRStarTheoretical()

    #Store median values to plot the track through E* R* space
    EStarMedian = np.zeros(NoBasins)
    RStarMedian = np.zeros(NoBasins)

    # Setup extent for data density calcs
    ESmin = np.log10(0.1)
    ESmax = np.log10(100.)
    RSmin = np.log10(0.05)
    RSmax = np.log10(1.5)

    # setup grid for density calcs
    ESgrid = np.logspace(ESmin,ESmax,(ESmax-ESmin)*100.)
    RSgrid = np.logspace(RSmin,RSmax,(RSmax-RSmin)*100.)

    #loop through the basins
    for Basin in range(0,NoBasins):
    #for Basin in range(0,1):

        # Get the hillslope data for the basin
        Data = CalculateEStarRStar(Basin)

        # Get the convex hull
        #Points = np.column_stack((Data.EStar,Data.RStar))
        #Hull = ConvexHull(Points)

        # calculate the 2D density of the data given
        #Counts,Xbins,Ybins=np.histogram2d(Data.EStar,Data.RStar,bins=100)
        #Counts = Counts.T
        #X,Y = np.meshgrid(Xbins,Ybins)
        #plt.pcolormesh(X,Y,Counts)

        # calculate gaussian kernel density
        Values = np.vstack([np.log10(Data.EStar), np.log10(Data.RStar)])
        Density = gaussian_kde(Values)

        ES,RS = np.meshgrid(np.log10(ESgrid),np.log10(RSgrid))
        Positions = np.vstack([ES.ravel(), RS.ravel()])

        # colour code by basin number
        colour = float(Basin)/float(NoBasins)

        Density = np.reshape(Density(Positions).T, ES.shape)
        Density /= np.max(Density)
        #plt.pcolormesh(10**ES,10**RS,Density,cmap=cm.Reds)

        plt.contour(10**ES,10**RS,Density,[0.2,],colors=[ColourMap(colour),],linewidths=1.,alpha=0.5)
        #plt.plot(Data.EStar,Data.RStar,'k.',ms=2,zorder=32)

        # make the contour plot
        #plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),
        #    ybins.min(),ybins.max()],linewidths=3,colors='black',
        #    linestyles='solid')

        # colour code by basin number
        #colour = float(Basin)/float(NoBasins)

        # Get median EStar RStar
        EStarMedian[Basin] = Data.EStar.median()
        RStarMedian[Basin] = Data.RStar.median()
        plt.plot(Data.EStar.median(),Data.RStar.median(),'o',ms=5,color=ColourMap(colour), zorder=32)

        # Plot the Hull
        #if Basin % 4 == 0:


        #    Ind = np.append(Hull.vertices, Hull.vertices[0])
        #    plt.plot(Points[Ind,0], Points[Ind,1], '-', color=ColourMap(colour), lw=1,alpha=0.5)

    #plot the path
    #plt.plot(EStarMedian,RStarMedian,'k-')

    # Finalise the figure
    plt.xlabel('$E^*={{-2\:C_{HT}\:L_H}/{S_C}}$')
    plt.ylabel('$R^*=S/S_C$')
    plt.xlim(0.1,100)
    plt.ylim(0.05,1.5)

    # add colour bar
    m = cm.ScalarMappable(cmap=ColourMap)
    m.set_array(np.arange(0,NoBasins))
    cbar = plt.colorbar(m)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.set_label('Basin No.')

    #save the figure
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(PlotDirectory+FilenamePrefix + "_EStarRStarProgression.png", dpi=300)
    plt.close(Fig)

# SMM Checked, revised and working 13/06/2018
# This is now working and produces lots of nice profile plots.
# It could be revised so that it shows the main stem only
def PlotChiProfileHillslopeData(DataDirectory, FilenamePrefix, PlotDirectory, Basins = [], PlotKsn = False, Sc = 0.71, mainstem_only = False, minimum_traces = 50,
                               common_max_Es = -99, common_max_Ksn = -99):
    """
    This plots the data by basin showing the E*, R* and either the chi profile or the K_sn data as a function of chi

    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): the file name prefix
        PlotDirectory (str): The directory into which the plots are saved
        Basins (list): The basins to be plotted
        PlotKsn (bool): If true, the profile plot will be Ksn instead of elevation
        Sc (float): The critical slope
        mainstem_only (bool): If true, only plot the data from the main stem
        minimum_traces (int): The minimum number of traces required to plot the hillslope data
        common_max_Es (float): If this is positive, use as the maximum Es for all plots
        common_max_Ksn (float): If this is positive, use as the maximum Ksn for all plots

    Author: MDH

    """

    print("Hi there. Let me print some basin by basin plots for you.")
    if PlotKsn:
        print("You are plotting chi-k_sn rather than chi-elevation")
    else:
        print("You are plotting chi-elevation rather than chi-k_sn")

    #Load hillslope metrics data
    HillslopesDF = ReadHillslopeData(DataDirectory, FilenamePrefix)

    # Read in the raw channel data
    ChannelsDF = ReadChannelData(DataDirectory, FilenamePrefix)

    # Basins list and keys
    BasinsDict = np.loadtxt(DataDirectory+FilenamePrefix+'_junctions.list',dtype=int)

    # loop through basins
    for key in Basins:

        Basin = BasinsDict[key]

        # print basin to screen
        print(key, Basin)

        # isolate basin data
        BasinChannelData = ChannelsDF[ChannelsDF.basin_key == key]
        MinimumChi = BasinChannelData.chi.min()
        MaximumMChi = BasinChannelData.m_chi.max()
        MinKsn = BasinChannelData.m_chi.min()
        MaxKsn = BasinChannelData.m_chi.max()

        # how many segments are we dealing with?
        Segments = BasinChannelData.segment_number.unique()


        # try to figure out the source key
        mainstem_source_key = BasinChannelData.source_key.iloc[0]
        print("The mainstem source key is: " +str(mainstem_source_key))

        # separate into main stem and trib data
        MainStemChannelData = BasinChannelData[BasinChannelData.source_key == mainstem_source_key]
        MainStemSegments = MainStemChannelData.segment_number.unique()

        # setup the figure
        Fig = CreateFigure(FigSizeFormat="EPSL",AspectRatio=1)
        ax1 = Fig.add_axes([0.1,0.1,0.8,0.35])
        ax2 = Fig.add_axes([0.1,0.45,0.8,0.32])
        ax3 = Fig.add_axes([0.1,0.65,0.8,0.32])

        #choose colormap
        ColourMap = cm.viridis

        # create new dataframe for plotting
        PlotDF = pd.DataFrame(columns=['Chi','Ksn','EStarMedian','EStarLower',
                     'EStarUpper','RStarMedian','RStarLower','RStarUpper','NTraces'])

        # Get the data columns for plotting
        for i, Segment in np.ndenumerate(Segments):


            if mainstem_only:
                if Segments[i] in MainStemSegments:
                    # get metrics to plot
                    if PlotKsn:
                        KKsn = BasinChannelData.m_chi[BasinChannelData.segment_number == Segment]

                    Ksn = BasinChannelData.m_chi[BasinChannelData.segment_number == Segment].unique()[0]
                    Chi = BasinChannelData.chi[BasinChannelData.segment_number == Segment]
                    Elevation = BasinChannelData.elevation[BasinChannelData.segment_number == Segment]

                    #print("Sizes are:")
                    #print("Ksn: "+str(Ksn.size))
                    #print("Chi: "+str(Chi.size))
                    #print("Elevation: "+str(Elevation.size))

                    #normalise chi by outlet chi
                    Chi = Chi-MinimumChi

                    # plot the chi data
                    Colour = (Ksn-MinKsn)/(MaxKsn-MinKsn)

                    PlotMaxKsn = int(math.ceil(MaxKsn / 10.0)) * 10

                    if PlotKsn:
                        ax1.scatter(Chi,KKsn,marker='o', edgecolors='none', lw=0.5, c=[1.0,0.0,0.0], s=20, zorder=20)
                    else:
                        ax1.plot(Chi,Elevation,'-', lw=1.5,c=ColourMap(Colour), zorder=10)

                    # get hillslope data
                    SegmentHillslopeData = HillslopesDF[HillslopesDF.StreamID == Segment]
                    NTraces = len(SegmentHillslopeData["i"].tolist())

                    if NTraces<minimum_traces:
                        continue

                    # Calculate E* R*
                    EStar = -2*SegmentHillslopeData.Cht*SegmentHillslopeData.Lh/Sc
                    RStar = SegmentHillslopeData.S/Sc

                    EStarMedian = EStar.median()
                    EStarLower = EStar.quantile(0.25)
                    EStarUpper = EStar.quantile(0.75)

                    RStarMedian = RStar.median()
                    RStarLower = RStar.quantile(0.25)
                    RStarUpper = RStar.quantile(0.75)

                    # add to plot dataframe
                    PlotDF.loc[i]  = [Chi.median(),Ksn,EStarMedian,EStarLower, EStarUpper, RStarMedian, RStarLower, RStarUpper,NTraces]

            else:
                # get metrics to plot
                if PlotKsn:
                    KKsn = BasinChannelData.m_chi[BasinChannelData.segment_number == Segment]

                Ksn = BasinChannelData.m_chi[BasinChannelData.segment_number == Segment].unique()[0]
                Chi = BasinChannelData.chi[BasinChannelData.segment_number == Segment]
                Elevation = BasinChannelData.elevation[BasinChannelData.segment_number == Segment]

                #print("Sizes are:")
                #print("Ksn: "+str(Ksn.size))
                #print("Chi: "+str(Chi.size))
                #print("Elevation: "+str(Elevation.size))

                #normalise chi by outlet chi
                Chi = Chi-MinimumChi

                # plot the chi data
                Colour = (Ksn-MinKsn)/(MaxKsn-MinKsn)

                PlotMaxKsn = int(math.ceil(MaxKsn / 10.0)) * 10

                if PlotKsn:
                    ax1.scatter(Chi,KKsn,marker='o', edgecolors='none',lw=0.5, c=[1.0,0.0,0.0], s=20, zorder=20)
                else:
                    ax1.plot(Chi,Elevation,'-', lw=1.5,c=ColourMap(Colour), zorder=10)

                # get hillslope data
                SegmentHillslopeData = HillslopesDF[HillslopesDF.StreamID == Segment]
                NTraces = len(SegmentHillslopeData["i"].tolist())

                if NTraces<minimum_traces:
                    continue

                # Calculate E* R*
                EStar = -2*SegmentHillslopeData.Cht*SegmentHillslopeData.Lh/Sc
                RStar = SegmentHillslopeData.S/Sc

                EStarMedian = EStar.median()
                EStarLower = EStar.quantile(0.25)
                EStarUpper = EStar.quantile(0.75)

                RStarMedian = RStar.median()
                RStarLower = RStar.quantile(0.25)
                RStarUpper = RStar.quantile(0.75)

                # add to plot dataframe
                PlotDF.loc[i]  = [Chi.median(),Ksn,EStarMedian,EStarLower, EStarUpper, RStarMedian, RStarLower, RStarUpper,NTraces]

        # reset indices
        PlotDF = PlotDF.reset_index(drop=True)

        # Zip errors for plotting
        Es_max_err = PlotDF.EStarUpper.values-PlotDF.EStarMedian
        Es_min_err = PlotDF.EStarMedian-PlotDF.EStarLower.values
        Es_errors = np.array(zip(Es_min_err, Es_max_err)).T
        Rs_max_err = PlotDF.RStarUpper.values-PlotDF.RStarMedian
        Rs_min_err = PlotDF.RStarMedian-PlotDF.RStarLower.values
        Rs_errors = np.array(zip(Rs_min_err, Rs_max_err)).T

        #Get colours for plotting from Chi

        #plot ksn vs EStar and Rstar, colouring by Chi
        for i, row in PlotDF.iterrows():
            ax2.plot([row.Chi,row.Chi],[row.EStarLower, row.EStarUpper],'-',c=[0.5,0.9,0.7],lw=2)
            ax2.scatter(row.Chi, row.EStarMedian, marker='o', edgecolors='k',lw=0.5, c=[0.5,0.9,0.7], s=20, zorder=200)
            ax3.plot([row.Chi,row.Chi],[row.RStarLower, row.RStarUpper],'-',c=[0.5,0.7,0.9],lw=2)
            ax3.scatter(row.Chi, row.RStarMedian, marker='o', edgecolors='k',lw=0.5, c=[0.5,0.7,0.9], s=20, zorder=200)

        # Finalise the figure
        if PlotKsn:
            ax1.set_ylabel(r"$k_{sn}$")
        else:
            ax1.set_ylabel('Elevation (m)')

        ax1.set_xlabel(r"$\chi$ (m)")
        ax2.set_ylabel('Dimensionless $C_{\mathit{HT}}$')
        ax3.set_ylabel('Dimensionless Relief $(S/S_C)$')


        #add colourbar if you have a profile plot
        if not PlotKsn:
            CAx = Fig.add_axes([0.6,0.17,0.25,0.02])
            m = cm.ScalarMappable(cmap=ColourMap)
            m.set_array(PlotDF.Ksn)
            plt.colorbar(m, cax=CAx,orientation='horizontal')
            plt.xlabel('$k_{sn}$',fontsize=8)
            CAx.tick_params(axis='both', labelsize=8)

        # turn off ax2 overlap and x axis for superimposed plots
        ax1.patch.set_facecolor('none')
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.yaxis.set_ticks_position('left')
        ax1.xaxis.set_ticks_position('bottom')
        ax2.patch.set_facecolor('none')
        ax2.spines['left'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        ax2.xaxis.set_visible(False)
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")
        ax3.patch.set_facecolor('none')
        ax3.spines['right'].set_visible(False)
        ax3.spines['top'].set_visible(False)
        ax3.spines['bottom'].set_visible(False)
        ax3.xaxis.set_visible(False)
        ax3.yaxis.set_label_position("left")
        ax3.yaxis.set_ticks_position('left')

        # fix axis limits
        if PlotKsn:
            if common_max_Ksn > 0:
                ax1.set_ylim(0,common_max_Ksn)
            else:
                ax1.set_ylim(0,PlotMaxKsn)

        if common_max_Es > 0:
            ax2.set_ylim(0,common_max_Es)
        else:
            ax2.set_ylim(0,PlotDF.EStarUpper.max())

        ax3.set_ylim(0,1)

        #save output
        plt.suptitle('Basin ID ' + str(key) + " (" + str(Basin) + ")")

        if mainstem_only:
            if PlotKsn:
                plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(key).zfill(3) + "_ChiProfileEsRs_Ksn_ms.png", dpi=300)
            else:
                plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(key).zfill(3) + "_ChiProfileEsRs_ms.png", dpi=300)
        else:
            if PlotKsn:
                plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(key).zfill(3) + "_ChiProfileEsRs_Ksn.png", dpi=300)
            else:
                plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(key).zfill(3) + "_ChiProfileEsRs.png", dpi=300)
        plt.clf()
        plt.close()

# This has been taken from one of Martin's scripts
# Tested and working as of 13-6-2018 (SMM)
def PlotCatchmentKsnEsRs(DataDirectory, FilenamePrefix,PlotDirectory, Basins = [], Sc = 0.71, mainstem_only = False, minimum_traces = 20):
    """
    This prints plots of k_sn vs E* and R* for each basin. It colours points by the chi coordinate/

    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): the file name prefix
        PlotDirectory (str): The directory into which the plots are saved
        Basins (int list): A list of the basin numbers to plot
        Sc (float): The critical slope
        mainstem_only (bool): If true, only plot the data from the main stem
        minimum_traces (int): The minimum number of traces required to plot the hillslope data

    Author:
        MDH
        SMM (modified 13-06-2018)

    """

    #Load hillslope metrics data
    HillslopesDF = ReadHillslopeData(DataDirectory, FilenamePrefix)

    # Read in the raw channel data
    ChannelsDF = ReadChannelData(DataDirectory, FilenamePrefix)

    # Basins list and keys
    BasinsDict = np.loadtxt(DataDirectory+FilenamePrefix+'_junctions.list',dtype=int)

    # loop through basins
    #for key, Basin in np.ndenumerate(Basins):
    for key in Basins:
        # print basin to screen
        Basin = BasinsDict[key]
        print(key, Basin)

        # isolate basin data
        BasinChannelData = ChannelsDF[ChannelsDF.basin_key == key]
        MinimumChi = BasinChannelData.chi.min()

        # how many segments are we dealing with?
        Segments = BasinChannelData.segment_number.unique()

        # setup the figure
        Fig = CreateFigure(FigSizeFormat="EPSL")
        ax1 = Fig.add_axes([0.1,0.1,0.8,0.5])
        ax2 = Fig.add_axes([0.1,0.45,0.8,0.5])

        #choose colormap
        ColourMap = cm.viridis

        # create new dataframe for plotting
        PlotDF = pd.DataFrame(columns=['Chi','Ksn','EStarMedian','EStarLower',
                    'EStarUpper','RStarMedian','RStarLower','RStarUpper','NTraces'])

        # try to figure out the source key
        mainstem_source_key = BasinChannelData.source_key.iloc[0]
        print("The mainstem source key is: " +str(mainstem_source_key))

        # separate into main stem and trib data
        MainStemChannelData = BasinChannelData[BasinChannelData.source_key == mainstem_source_key]
        MainStemSegments = MainStemChannelData.segment_number.unique()

        # Get the data columns for plotting
        for i, Segment in np.ndenumerate(Segments):

            # A rather stupid way to ensure only mainstem but I don't have time to make this pythonic (SMM)
            if mainstem_only:
                if Segments[i] in MainStemSegments:
                    # get metrics to plot
                    Ksn = BasinChannelData.m_chi[BasinChannelData.segment_number == Segment].unique()[0]
                    Chi = BasinChannelData.chi[BasinChannelData.segment_number == Segment].median()

                    #normalise chi by outlet chi
                    Chi = Chi-MinimumChi

                    # get hillslope data
                    SegmentHillslopeData = HillslopesDF[HillslopesDF.StreamID == Segment]
                    NTraces = len(SegmentHillslopeData)

                    if NTraces<minimum_traces:
                        continue

                    # Calculate E* R*
                    EStar = -2*SegmentHillslopeData.Cht*SegmentHillslopeData.Lh/Sc
                    RStar = SegmentHillslopeData.S/Sc

                    EStarMedian = EStar.median()
                    EStarLower = EStar.quantile(0.25)
                    EStarUpper = EStar.quantile(0.75)

                    RStarMedian = RStar.median()
                    RStarLower = RStar.quantile(0.25)
                    RStarUpper = RStar.quantile(0.75)

                    # add to plot dataframe
                    PlotDF.loc[i]  = [Chi,Ksn,EStarMedian,EStarLower, EStarUpper, RStarMedian, RStarLower, RStarUpper,NTraces]

            else:
                # get metrics to plot
                Ksn = BasinChannelData.m_chi[BasinChannelData.segment_number == Segment].unique()[0]
                Chi = BasinChannelData.chi[BasinChannelData.segment_number == Segment].median()

                #normalise chi by outlet chi
                Chi = Chi-MinimumChi

                # get hillslope data
                SegmentHillslopeData = HillslopesDF[HillslopesDF.StreamID == Segment]
                NTraces = len(SegmentHillslopeData)

                if NTraces<minimum_traces:
                    continue

                # Calculate E* R*
                EStar = -2*SegmentHillslopeData.Cht*SegmentHillslopeData.Lh/Sc
                RStar = SegmentHillslopeData.S/Sc

                EStarMedian = EStar.median()
                EStarLower = EStar.quantile(0.25)
                EStarUpper = EStar.quantile(0.75)

                RStarMedian = RStar.median()
                RStarLower = RStar.quantile(0.25)
                RStarUpper = RStar.quantile(0.75)

                # add to plot dataframe
                PlotDF.loc[i]  = [Chi,Ksn,EStarMedian,EStarLower, EStarUpper, RStarMedian, RStarLower, RStarUpper,NTraces]

        # reset indices
        PlotDF = PlotDF.reset_index(drop=True)

        # Zip errors for plotting
        Es_max_err = PlotDF.EStarUpper.values-PlotDF.EStarMedian
        Es_min_err = PlotDF.EStarMedian-PlotDF.EStarLower.values
        Es_errors = np.array(zip(Es_min_err, Es_max_err)).T
        Rs_max_err = PlotDF.RStarUpper.values-PlotDF.RStarMedian
        Rs_min_err = PlotDF.RStarMedian-PlotDF.RStarLower.values
        Rs_errors = np.array(zip(Rs_min_err, Rs_max_err)).T

        #Get colours for plotting from Chi
        ChiArray = PlotDF.Chi.values.astype(float)
        MinChi = PlotDF.Chi.min()
        MaxChi = PlotDF.Chi.max()
        Colours = (ChiArray-MinChi)/(MaxChi-MinChi)

        #plot ksn vs EStar and Rstar, colouring by Chi
        for i, row in PlotDF.iterrows():
            ax1.plot([row.Ksn,row.Ksn],[row.EStarLower, row.EStarUpper],'-',c=ColourMap(Colours[i]))
            ax1.scatter(row.Ksn, row.EStarMedian, marker='o', edgecolors='k',lw=0.5, facecolors=ColourMap(Colours[i]), s=15, zorder=200)
            ax2.plot([row.Ksn,row.Ksn],[row.RStarLower, row.RStarUpper],'-',c=ColourMap(Colours[i]),lw=2)
            ax2.scatter(row.Ksn, row.RStarMedian, marker='o', edgecolors='k',lw=0.5, facecolors=ColourMap(Colours[i]), s=15, zorder=200)

        # Finalise the figure
        ax1.set_xlabel(r"$K_{sn}$ (m$^{0.62}$)")
        ax1.set_ylabel('Dimensionless $C_{\mathit{HT}}$')
        ax2.set_ylabel('Dimensionless Relief $(S/S_C)$')

        #add colourbar
        CAx = Fig.add_axes([0.02,0.9,0.2,0.02])
        m = cm.ScalarMappable(cmap=ColourMap)
        m.set_array(PlotDF.Chi)
        plt.colorbar(m, cax=CAx,orientation='horizontal')
        plt.xlabel('${\chi}$ (m)',fontsize=8)
        CAx.tick_params(axis='both', labelsize=8)

        # turn off ax2 overlap and x axis for superimposed plots
        ax1.patch.set_facecolor('none')
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.yaxis.set_ticks_position('left')
        ax1.xaxis.set_ticks_position('bottom')
        ax2.patch.set_facecolor('none')
        ax2.spines['left'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        ax2.xaxis.set_visible(False)
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")

        # fix axis limits
        ax1.set_ylim(0,25)
        ax2.set_ylim(0,1)

        #save output
        plt.suptitle('Basin ID ' + str(key) + " (" + str(Basin) + ")")

        if mainstem_only:
            plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(key).zfill(3) + "_KsnEsRs_ms.png", dpi=300)
        else:
            plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(key).zfill(3) + "_KsnEsRs.png", dpi=300)

        plt.clf()
        plt.close()

# This has been taken from one of Martin's scripts
# Tested and working as of 13-6-2018 (SMM)
def PlotStackedEsRsFxnChi(DataDirectory, FilenamePrefix,PlotDirectory, Basins = [], Sc = 0.71, mainstem_only = False):
    """
    This plots the E* and R* data as a function of where they are in chi space

    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): the file name prefix
        PlotDirectory (str): The directory into which the plots are saved
        Basins (int list): A list of the basin numbers to plot
        Sc (float): The critical slope
        mainstem_only (bool): If true, only plot the data from the main stem


    Author: SMM

    Date: 15-Jun-2018

    """

    #Load hillslope metrics data
    HillslopesDF = ReadHillslopeData(DataDirectory, FilenamePrefix)

    # Read in the raw channel data
    ChannelsDF = ReadChannelData(DataDirectory, FilenamePrefix)

    # Basins list and keys
    BasinsDict = np.loadtxt(DataDirectory+FilenamePrefix+'_junctions.list',dtype=int)

    # Create a dictionary for storing the plotting data
    PlotDataDict = {}

    # loop through basins
    #for key, Basin in np.ndenumerate(Basins):
    for key in Basins:
        # print basin to screen
        Basin = BasinsDict[key]
        print(key, Basin)

        # isolate basin data
        BasinChannelData = ChannelsDF[ChannelsDF.basin_key == key]
        MinimumChi = BasinChannelData.chi.min()

        # how many segments are we dealing with?
        Segments = BasinChannelData.segment_number.unique()



        # create new dataframe for plotting
        PlotDF = pd.DataFrame(columns=['Chi','Ksn','EStarMedian','EStarLower',
                    'EStarUpper','RStarMedian','RStarLower','RStarUpper','NTraces'])


        # try to figure out the source key
        mainstem_source_key = BasinChannelData.source_key.iloc[0]
        print("The mainstem source key is: " +str(mainstem_source_key))

        # separate into main stem and trib data
        MainStemChannelData = BasinChannelData[BasinChannelData.source_key == mainstem_source_key]
        MainStemSegments = MainStemChannelData.segment_number.unique()

        # Get the data columns for plotting
        for i, Segment in np.ndenumerate(Segments):

            if mainstem_only:
                print("I am not gonna do a thing.")
            else:


                # get hillslope data
                SegmentHillslopeData = HillslopesDF[HillslopesDF.StreamID == Segment]
                NTraces = len(SegmentHillslopeData)

                if NTraces>0:

                    # get metrics to plot
                    Ksn = BasinChannelData.m_chi[BasinChannelData.segment_number == Segment].unique()[0]
                    Chi = BasinChannelData.chi[BasinChannelData.segment_number == Segment].median()

                    #normalise chi by outlet chi
                    Chi = Chi-MinimumChi

                    # Calculate E* R*
                    EStar = -2*SegmentHillslopeData.Cht*SegmentHillslopeData.Lh/Sc
                    RStar = SegmentHillslopeData.S/Sc

                    EStarMedian = EStar.median()
                    EStarLower = EStar.quantile(0.25)
                    EStarUpper = EStar.quantile(0.75)

                    RStarMedian = RStar.median()
                    RStarLower = RStar.quantile(0.25)
                    RStarUpper = RStar.quantile(0.75)

                    # add to plot dataframe
                    PlotDF.loc[i]  = [Chi,Ksn,EStarMedian,EStarLower, EStarUpper, RStarMedian, RStarLower, RStarUpper,NTraces]

            # reset indices
            PlotDF = PlotDF.reset_index(drop=True)

        # Zip errors for plotting
        Es_max_err = PlotDF.EStarUpper.values-PlotDF.EStarMedian
        Es_min_err = PlotDF.EStarMedian-PlotDF.EStarLower.values
        Es_errors = np.array(zip(Es_min_err, Es_max_err)).T
        Rs_max_err = PlotDF.RStarUpper.values-PlotDF.RStarMedian
        Rs_min_err = PlotDF.RStarMedian-PlotDF.RStarLower.values
        Rs_errors = np.array(zip(Rs_min_err, Rs_max_err)).T

        #Get colours for plotting from Chi
        ChiArray = PlotDF.Chi.values.astype(float)
        MinChi = PlotDF.Chi.min()
        MaxChi = PlotDF.Chi.max()
        Colours = (ChiArray-MinChi)/(MaxChi-MinChi)

        # order the data by chi
        Sorted_PlotDF = PlotDF.sort_values(by=['Chi'])

        # Add to a dict of
        PlotDataDict[key] = Sorted_PlotDF

    # Now we make the plot

    #choose colormap
    ColourMap = cm.viridis

    Fig = CreateFigure(FigSizeFormat="EPSL", AspectRatio = 0.6)
    ax1 = Fig.add_axes([0.1,0.1,0.8,0.8])

    # Add an offset counter
    offset = 0
    offsetter = 0.8


    for key in PlotDataDict:
        ThisPlotDF = PlotDataDict[key]
        Chi_data = ThisPlotDF.as_matrix(columns = ["Chi"])[:,0]
        Es_data = ThisPlotDF.as_matrix(columns = ["EStarMedian"])[:,0]
        Rs_data = ThisPlotDF.as_matrix(columns = ["RStarMedian"])[:,0]

        chi = Chi_data.astype(float)
        Es = Es_data.astype(float)
        Rs = Rs_data.astype(float)
        chi = chi.tolist()
        Es = Es.tolist()
        Rs = Rs.tolist()
        #filler = np.asarray(chi)
        #filler.fill(offset)
        #Es = Es_data+offset

        # This fucking shit is required because of a fucking crazy error I about utypes I can only fix this way
        c = []
        e = []
        f = []
        r = []
        for i in range(0,len(chi)):
            c.append(float(chi[i]))
            e.append(float(Es[i]))
            r.append(float(Rs[i]))
            f.append(float(offset))

        C = np.asarray(c)
        E = np.asarray(e)
        F = np.asarray(f)
        R = np.asarray(r)
        EE = E+offset
        RR = R+offset

        max_Es =EE.max()
        max_Rs = RR.max()

        #plt.plot(C,EE, color="k", alpha = 0.8)
        #plt.fill_between(C,F,EE, facecolor='red', interpolate=True, alpha = 0.5)
        plt.plot(C,RR, color="k", alpha = 0.8)
        plt.fill_between(C,F,RR, facecolor='red', interpolate=True, alpha = 0.5)

        offset = offset+offsetter

    # Finalise the figure
    ax1.set_xlabel(r"$\chi$ (m)")
    #ax1.set_ylabel('Dimensionless $C_{\mathit{HT}}$')
    ax1.set_ylabel('Dimensionless Relief')

    # turn off ax2 overlap and x axis for superimposed plots
    ax1.patch.set_facecolor('none')
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')

    # fix axis limits
    #ax1.set_ylim(0,max_Es)
    ax1.set_ylim(0,max_Rs)
    ax1.set_xlim(0,60)

    plt.savefig(PlotDirectory+FilenamePrefix + "_Stack.png", dpi=300)

    plt.clf()
    plt.close()


def GetClusteredDataPlotDict(DataDirectory, FilenamePrefix, Sc = 0.71, mainstem_only = False, BasinsCluster = []):
    """
    This function reads the hillslope and channel data and returns a data dict that can be used for plotting

    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): the file name prefix
        Sc (float): The critical slope
        mainstem_only (bool): If true, only plot the data from the main stem
        BasinsCluster (list of int lists): This is a list of lists that has the basin numbers for clustering

    Author: SMM

    Returns:
        BasinCluster (list): A list of integer lists.
    """


    # Basins list and keys
    BasinsDict = np.loadtxt(DataDirectory+FilenamePrefix+'_junctions.list',dtype=int)
    print("The basins dict is: ")
    print(BasinsDict)

    Sorted = range(0,len(BasinsDict))
    #Sorted =  range(0,6)

    if len(BasinsCluster) == 0:
        print("You haven't supplied me with a list of basins so I am going to create a series of sequential lists")
        print("These will be of different sizes")
        seq = Sorted

        # Make a sequantial cluster
        combined_list = []
        for i in range(6,8):
            this_listlist = chunkIt(seq, i)
            for a_list in this_listlist:
                combined_list.append(a_list)

        print("Combined list is:")
        print(combined_list)
        BasinsCluster = combined_list


    #Load hillslope metrics data
    HillslopesDF = ReadHillslopeData(DataDirectory, FilenamePrefix)

    # Read in the raw channel data
    ChannelsDF = ReadChannelData(DataDirectory, FilenamePrefix)


    # Create a dictionary for storing the plotting data
    PlotDataDict = {}
    for cluster_index,Basins in enumerate(BasinsCluster):
        # Each basin cluster has several basins in it
        print("This cluster has the following basins:")
        print(Basins)

        # create new dataframe for plotting
        PlotDF = pd.DataFrame(columns=['Chi','Ksn','EStarMedian','EStarLower',
                        'EStarUpper','RStarMedian','RStarLower','RStarUpper','NTraces'])

        # loop through basins
        #for key, Basin in np.ndenumerate(Basins):
        for key in Basins:
            # print basin to screen
            Basin = BasinsDict[key]
            #print(key, Basin)

            # isolate basin data
            BasinChannelData = ChannelsDF[ChannelsDF.basin_key == key]
            MinimumChi = BasinChannelData.chi.min()

            # how many segments are we dealing with?
            Segments = BasinChannelData.segment_number.unique()

            # try to figure out the source key
            mainstem_source_key = BasinChannelData.source_key.iloc[0]
            #print("The mainstem source key is: " +str(mainstem_source_key))

            # separate into main stem and trib data
            MainStemChannelData = BasinChannelData[BasinChannelData.source_key == mainstem_source_key]
            MainStemSegments = MainStemChannelData.segment_number.unique()

            # Get the data columns for plotting
            for i, Segment in np.ndenumerate(Segments):

                if mainstem_only:
                    if Segment in MainStemSegments:
                        # get hillslope data
                        SegmentHillslopeData = HillslopesDF[HillslopesDF.StreamID == Segment]
                        NTraces = len(SegmentHillslopeData)

                        if NTraces>20:

                            # get metrics to plot
                            Ksn = BasinChannelData.m_chi[BasinChannelData.segment_number == Segment].unique()[0]
                            Chi = BasinChannelData.chi[BasinChannelData.segment_number == Segment].median()
                            BasinKey = BasinChannelData.basin_key[BasinChannelData.segment_number == Segment].median()

                            #normalise chi by outlet chi
                            Chi = Chi-MinimumChi

                            # Calculate E* R*
                            EStar = -2*SegmentHillslopeData.Cht*SegmentHillslopeData.Lh/Sc
                            RStar = SegmentHillslopeData.S/Sc

                            EStarMedian = EStar.median()
                            EStarLower = EStar.quantile(0.25)
                            EStarUpper = EStar.quantile(0.75)

                            RStarMedian = RStar.median()
                            RStarLower = RStar.quantile(0.25)
                            RStarUpper = RStar.quantile(0.75)

                            # add to plot dataframe
                            PlotDF.loc[i]  = [Chi,Ksn,EStarMedian,EStarLower, EStarUpper, RStarMedian, RStarLower, RStarUpper,NTraces]
                else:

                    # get hillslope data
                    SegmentHillslopeData = HillslopesDF[HillslopesDF.StreamID == Segment]
                    NTraces = len(SegmentHillslopeData)

                    if NTraces>20:

                        # get metrics to plot
                        Ksn = BasinChannelData.m_chi[BasinChannelData.segment_number == Segment].unique()[0]
                        Chi = BasinChannelData.chi[BasinChannelData.segment_number == Segment].median()

                        #normalise chi by outlet chi
                        Chi = Chi-MinimumChi

                        # Calculate E* R*
                        EStar = -2*SegmentHillslopeData.Cht*SegmentHillslopeData.Lh/Sc
                        RStar = SegmentHillslopeData.S/Sc

                        EStarMedian = EStar.median()
                        EStarLower = EStar.quantile(0.25)
                        EStarUpper = EStar.quantile(0.75)

                        RStarMedian = RStar.median()
                        RStarLower = RStar.quantile(0.25)
                        RStarUpper = RStar.quantile(0.75)

                        # add to plot dataframe
                        PlotDF.loc[i]  = [Chi,Ksn,EStarMedian,EStarLower, EStarUpper, RStarMedian, RStarLower, RStarUpper,NTraces]

        # reset indices
        PlotDF = PlotDF.reset_index(drop=True)

        print("The cluster index is: "+str(cluster_index))
        PlotDataDict[cluster_index] = PlotDF

    return BasinsCluster,PlotDataDict

# Tested and working as of 15-6-2018 (SMM)
def PlotClusteredEsRsFxnChi(DataDirectory, FilenamePrefix,PlotDirectory, Sc = 0.71, mainstem_only = False, BasinsCluster = []):
    """
    This plots the E* and R* data as a function of where they are in chi space

    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): the file name prefix
        PlotDirectory (str): The directory into which the plots are saved
        Sc (float): The critical slope
        mainstem_only (bool): If true, only plot the data from the main stem
        BasinsCluster (list of int lists): This is a list of lists that has the basin numbers for clustering

    Author: SMM

    Date: 15-Jun-2018
    """

    BasinsCluster,PlotDataDict = GetClusteredDataPlotDict(DataDirectory, FilenamePrefix, Sc, mainstem_only, BasinsCluster)

    for key in PlotDataDict:

        # Each basin cluster has several basins in it
        Basins = BasinsCluster[key]
        print("This cluster has the following basins:")
        print(Basins)

        # For labelling
        #StrBasins = str(Basins)
        clusterstr = ",".join(str(i) for i in Basins)

        # setup the figure
        Fig1 = CreateFigure(FigSizeFormat="EPSL")
        ax1 = Fig1.add_axes([0.1,0.1,0.8,0.7])

        # setup the figure
        Fig2 = CreateFigure(FigSizeFormat="EPSL")
        ax2 = Fig2.add_axes([0.1,0.1,0.8,0.7])

        #choose colormap
        ColourMap = cm.viridis

        # create new dataframe for plotting
        PlotDF = PlotDataDict[key]

        # Get the colourmap
        KsnArray = PlotDF.Ksn.values.astype(float)
        MinKsn = PlotDF.Ksn.min()
        MaxKsn = PlotDF.Ksn.max()
        #MinKsn = 0
        #MaxKsn = 20
        Colours1 = (KsnArray-MinKsn)/(MaxKsn-MinKsn)


         #plot ksn vs EStar and Rstar, colouring by Chi
        for i, row in PlotDF.iterrows():
            ax1.plot([row.Chi,row.Chi],[row.EStarLower, row.EStarUpper],'-',c=ColourMap(Colours1[i]))
            ax1.scatter(row.Chi, row.EStarMedian, marker='o', edgecolors='k',lw=0.5, facecolors=ColourMap(Colours1[i]), s=15, zorder=200)

            ax2.plot([row.Chi,row.Chi],[row.RStarLower, row.RStarUpper],'-',c=ColourMap(Colours1[i]))
            ax2.scatter(row.Chi, row.RStarMedian, marker='o', edgecolors='k',lw=0.5, facecolors=ColourMap(Colours1[i]), s=15, zorder=200)

        # Finalise the figure
        ax1.set_xlabel(r"$\chi$ (m)")
        ax1.set_ylabel('Dimensionless $C_{\mathit{HT}}$')

        # Finalise the figure
        ax2.set_xlabel(r"$\chi$ (m)")
        ax2.set_ylabel('Dimensionless relief')

        #add colourbar
        CAx = Fig1.add_axes([0.02,0.9,0.2,0.02])
        m = cm.ScalarMappable(cmap=ColourMap)
        m.set_array(PlotDF.Ksn)
        plt.colorbar(m, cax=CAx,orientation='horizontal')
        CAx.set_xlabel('${k_{sn}}$ (m)',fontsize=8)
        CAx.tick_params(axis='both', labelsize=8)

        #add colourbar
        CAx2 = Fig2.add_axes([0.02,0.9,0.2,0.02])
        m2 = cm.ScalarMappable(cmap=ColourMap)
        m2.set_array(PlotDF.Ksn)
        plt.colorbar(m2, cax=CAx2,orientation='horizontal')
        CAx2.set_xlabel('${k_{sn}}$ (m)',fontsize=8)
        CAx2.tick_params(axis='both', labelsize=8)

        # turn off ax2 overlap and x axis for superimposed plots
        ax1.patch.set_facecolor('none')
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.yaxis.set_ticks_position('left')
        ax1.xaxis.set_ticks_position('bottom')

        # turn off ax2 overlap and x axis for superimposed plots
        ax2.patch.set_facecolor('none')
        ax2.spines['right'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax2.yaxis.set_ticks_position('left')
        ax2.xaxis.set_ticks_position('bottom')


        # fix axis limits
        ax1.set_ylim(0,25)
        ax2.set_ylim(0,1)

        #save output
        Fig1.suptitle('Cluster number is ' + str(key)+ "\nBasins are: "+ clusterstr)
        Fig2.suptitle('Cluster number is ' + str(key)+ "\nBasins are: "+ clusterstr)

        if mainstem_only:
            Fig1.savefig(PlotDirectory+"Es_cluster_ms_"+FilenamePrefix + "_" + str(key).zfill(2) + ".png", dpi=300)
            Fig2.savefig(PlotDirectory+"Rs_cluster_ms_"+FilenamePrefix + "_" + str(key).zfill(2) + ".png", dpi=300)
        else:
            Fig1.savefig(PlotDirectory+"Es_cluster_"+FilenamePrefix + "_" + str(key).zfill(2) + ".png", dpi=300)
            Fig2.savefig(PlotDirectory+"Rs_cluster_"+FilenamePrefix + "_" + str(key).zfill(2) + ".png", dpi=300)

        # Clean up
        Fig1.clf()
        plt.close(Fig1)
        Fig2.clf()
        plt.close(Fig2)

# Working on this 15-06-2018
def PlotClusteredEsRs(DataDirectory, FilenamePrefix,PlotDirectory, Sc = 0.71, mainstem_only = False, BasinsCluster = [], colour_by = "chi"):
    """
    This plots the E* and R* data coloured by chi and ksn in clusters

    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): the file name prefix
        PlotDirectory (str): The directory into which the plots are saved
        Sc (float): The critical slope
        mainstem_only (bool): If true, only plot the data from the main stem
        BasinsCluster (list of int lists): This is a list of lists that has the basin numbers for clustering

    Author: SMM

    Date: 15-Jun-2018
    """

    BasinsCluster,PlotDataDict = GetClusteredDataPlotDict(DataDirectory, FilenamePrefix, Sc, mainstem_only, BasinsCluster)

    for key in PlotDataDict:

        # Each basin cluster has several basins in it
        Basins = BasinsCluster[key]
        print("This cluster has the following basins:")
        print(Basins)

        # For labelling
        #StrBasins = str(Basins)
        clusterstr = ",".join(str(i) for i in Basins)

        # setup the figure
        Fig1 = CreateFigure(FigSizeFormat="EPSL")
        ax1 = Fig1.add_axes([0.1,0.1,0.8,0.7])

        plt.loglog()
        PlotEStarRStarTheoretical()

        #choose colormap
        ColourMap = cm.viridis

        # create new dataframe for plotting
        PlotDF = PlotDataDict[key]

        # Get the colourmap
        KsnArray = PlotDF.Ksn.values.astype(float)
        MinKsn = PlotDF.Ksn.min()
        MaxKsn = PlotDF.Ksn.max()
        #MinKsn = 0
        #MaxKsn = 20
        Colours1 = (KsnArray-MinKsn)/(MaxKsn-MinKsn)

        ChiArray = PlotDF.Chi.values.astype(float)
        MinChi = PlotDF.Chi.min()
        MaxChi = PlotDF.Chi.max()
        #MinKsn = 0
        #MaxKsn = 20
        Colours2 = (ChiArray-MinChi)/(MaxChi-MinChi)


        EsArray = PlotDF.EStarMedian.values.astype(float)
        EsLArray = PlotDF.EStarLower.values.astype(float)
        EsLowErr = EsArray-EsLArray
        EsUArray = PlotDF.EStarUpper.values.astype(float)
        EsUpErr = EsUArray-EsArray

        RsArray = PlotDF.RStarMedian.values.astype(float)
        RsLArray = PlotDF.RStarLower.values.astype(float)
        RsLowErr = RsArray-RsLArray
        RsUArray = PlotDF.RStarUpper.values.astype(float)
        RsUpErr = RsUArray-RsArray

        #plot ksn vs EStar and Rstar, colouring by Chi
        ax1.errorbar(EsArray,RsArray,xerr=[EsLowErr, EsUpErr], yerr=[RsLowErr,RsUpErr],fmt='.', ecolor='k',markersize=2,mec='k',mfc='k',zorder = 10, linewidth = 1, alpha = 0.5)
        #ax2.errorbar(EsArray,RsArray,xerr=[EsLowErr, EsUpErr], yerr=[RsLowErr,RsUpErr],fmt='.', ecolor='k',markersize=2,mec='k',mfc='k',zorder = 10, linewidth = 1, alpha = 0.5)



        if colour_by == "chi":
            ax1.scatter(EsArray,RsArray,c=ChiArray,s=10, edgecolors='k', lw=0.1,cmap=ColourMap,zorder = 20)
        else:
            ax1.scatter(EsArray,RsArray,c=KsnArray,s=10, edgecolors='k', lw=0.1,cmap=ColourMap,zorder = 20)

        #add colourbar
        CAx = Fig1.add_axes([0.02,0.9,0.2,0.02])
        if colour_by == "chi":
            m = cm.ScalarMappable(cmap=ColourMap)
            m.set_array(PlotDF.Chi)
            plt.colorbar(m, cax=CAx,orientation='horizontal')
            CAx.set_xlabel('$\chi$ (m)',fontsize=8)
        else:
            m = cm.ScalarMappable(cmap=ColourMap)
            m.set_array(PlotDF.Ksn)
            plt.colorbar(m, cax=CAx,orientation='horizontal')
            CAx.set_xlabel('${k_{sn}}$ (m)',fontsize=8)

        CAx.tick_params(axis='both', labelsize=8)

        # Finalise the figure
        ax1.set_xlabel(r"Dimensionless $C_{\mathit{HT}}$")
        ax1.set_ylabel('Dimensionless Relief')

        # Finalise the figure
        #ax2.set_xlabel(r"$\chi$ (m)")
        #ax2.set_ylabel('Dimensionless relief')



        #add colourbar
        #CAx2 = Fig2.add_axes([0.02,0.9,0.2,0.02])
        #m2 = cm.ScalarMappable(cmap=ColourMap)
        #m2.set_array(PlotDF.Ksn)
        #plt.colorbar(m2, cax=CAx2,orientation='horizontal')
        #CAx2.set_xlabel('${k_{sn}}$ (m)',fontsize=8)
        #CAx2.tick_params(axis='both', labelsize=8)

        # turn off ax2 overlap and x axis for superimposed plots
        ax1.patch.set_facecolor('none')
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.yaxis.set_ticks_position('left')
        ax1.xaxis.set_ticks_position('bottom')

        # turn off ax2 overlap and x axis for superimposed plots
        #ax2.patch.set_facecolor('none')
        #ax2.spines['right'].set_visible(False)
        #ax2.spines['top'].set_visible(False)
        #ax2.yaxis.set_ticks_position('left')
        #ax2.xaxis.set_ticks_position('bottom')


        # fix axis limits
        #ax1.set_ylim(0,25)
        #ax2.set_ylim(0,1)

        #save output
        Fig1.suptitle('Cluster number is ' + str(key)+ "\nBasins are: "+ clusterstr)
        #Fig2.suptitle('Cluster number is ' + str(key)+ "\nBasins are: "+ clusterstr)

        if mainstem_only:
            Fig1.savefig(PlotDirectory+"Es_Rs_cluster_ms_"+FilenamePrefix + "_" + str(key).zfill(2) + ".png", dpi=300)
            #Fig2.savefig(PlotDirectory+"Rs_cluster_ms_"+FilenamePrefix + "_" + str(key).zfill(2) + ".png", dpi=300)
        else:
            Fig1.savefig(PlotDirectory+"Es_Rs_cluster_"+FilenamePrefix + "_" + str(key).zfill(2) + ".png", dpi=300)
            #Fig2.savefig(PlotDirectory+"Rs_cluster_"+FilenamePrefix + "_" + str(key).zfill(2) + ".png", dpi=300)

        # Clean up
        Fig1.clf()
        plt.close(Fig1)
        #Fig2.clf()
        #plt.close(Fig2)


def chunkIt(seq, num):
    """
    This comes from https://stackoverflow.com/questions/2130016/splitting-a-list-into-n-parts-of-approximately-equal-length
    I will use it to create a bunch of lists for sequential clustering

    Args:
        seq: The initial list for chunking
        num: The number of items in each chunk

    Return:
        A chunked list with roughly equal numbers of elements

    Author: Max Shawabkeh


    """
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out
