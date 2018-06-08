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
        ThinningFactor (int): This determines how many of the traces are discarded. Higher numbers mean more traces are discarded
        CustomExtent (list): if this is [-9999] the extent is grabbed from the raster. Otherwise you give it a 4 element list giving the extent of the area of interest. 


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
    Following Grieve et al. 2016 How Long is a Hillslope?

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
def PlotChiElevationSegments(DataDirectory, FilenamePrefix, PlotDirectory, BasinID):
    """


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
    """


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
    """


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
    """


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

def PlotHillslopeDataVsDistance(DataDirectory, FilenamePrefix, PlotDirectory, BasinID):
    """
    This function makes some composite plots of the hillslope data vs
    distance upstream from the outlet

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
        if SegmentHillslopeData.size > 50: # remove traces with < 50 pixels
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

def PlotEStarRStarWithinBasin(DataDirectory, FilenamePrefix, PlotDirectory, BasinID):
    """
    Makes a plot of E* against R* where the points are coloured by
    their distance from the outlet of the basin

    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): the file name prefix
        PlotDirectory (str): The directory into which the plots are saved
        BasinID (int): The basin to be plotted
        
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
        if SegmentHillslopeData.size > 50:
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
    uplift_rate = uplift_df['Uplift_rate']
    ax[6].plot(basin_keys, uplift_rate, c='k', ls='--')
    ax[6].set_ylabel('Uplift rate (mm/yr)')

    # set the axes labels
    ax[6].set_xlabel('Basin ID')
    plt.xticks(np.arange(min(basin_keys), max(basin_keys)+1, 1), rotation=45, fontsize=8)
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

def PlotEStarRStarBasins(DataDirectory, FilenamePrefix, PlotDirectory, Sc = 0.8):
    """
    Function to make an E*R* plot for a series of drainage basins.
    Changing so that we calculate E* and R* in the python script following
    Martin's example, so that we can test sensitivity to Sc.
    Author: FJC
    """
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

    FJC
    """
    input_csv = PlotDirectory+FilenamePrefix+'_basin_hillslope_data.csv'
    df = pd.read_csv(input_csv)

    # set up the figure
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10,5))

    #choose colormap
    ColourMap = cm.viridis

    # get the basins
    basins = df['basin_keys'].unique()
    NoBasins = len(basins)

    sc = ax[0].scatter(df.Estar_median,df.Rstar_median,c=df.basin_keys,s=50, edgecolors='k', zorder=100)
    ax[0].errorbar(df.Estar_median,df.Rstar_median,xerr=[df['Estar_lower_err'], df['Estar_upper_err']], yerr=[df['Rstar_lower_err'], df['Rstar_upper_err']],fmt='o', zorder=1, ecolor='0.5',markersize=1,mfc='white',mec='k')

    sc = ax[1].scatter(df.mchi_median,df.Rstar_median,c=df.basin_keys,s=50, edgecolors='k', zorder=100)
    ax[1].errorbar(df.mchi_median,df.Rstar_median,xerr=[df['mchi_lower_err'], df['mchi_upper_err']], yerr=[df['Rstar_lower_err'], df['Rstar_upper_err']],fmt='o', zorder=1, ecolor='0.5',markersize=1,mfc='white',mec='k')

    # plot the theoretical relationships for each one
    # Calculate analytical relationship for estar rstar
    EStar_model = np.logspace(-1,3,1000)
    RStar_model = CalculateRStar(EStar_model)

    # Plot with open figure
    ax[0].plot(EStar_model,RStar_model,c='0.5',ls='--')

    # calculate linear fit for Rstar ksn
    slope, intercept, r_value, p_value, std_err = stats.linregress(df.mchi_median, df.Rstar_median)
    print (slope, intercept, r_value, p_value)
    x = np.linspace(0, 200, 100)
    new_y = slope*x + intercept
    ax[1].plot(x, new_y, c='0.5', ls='--')

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
    cbar.set_label('Basin ID')

    #save output
    plt.savefig(PlotDirectory+FilenamePrefix +"_estar_rstar_subplots.png", dpi=300)
    plt.clf()


def PlotHillslopeTraces(DataDirectory, FilenamePrefix, PlotDirectory, CustomExtent=[-9999],FigSizeFormat="epsl"):
    """
    Function to plot a hillshade image with hilltops, hillslope traces and the channel network superimposed.


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
    Plots the progression of hillslopes along Bolinas in Estar Rstar space

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
    
def PlotChiProfileHillslopeData(DataDirectory, FilenamePrefix, PlotDirectory, Basins = [], PlotKsn = False, Sc = 0.71):
    """
    This plots the data by basin showing the E*, R* and either the chi profile or the K_sn data as a function of chi
    
    Args:
        DataDirectory (str): the data directory
        FilenamePrefix (str): the file name prefix
        PlotDirectory (str): The directory into which the plots are saved
        Basins (list): The basins to be plotted
        PlotKsn (bool): If true, the profile plot will be Ksn instead of elevation
        Sc (float): The critical slope
    
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
    Basins = np.loadtxt(DataDirectory+FilenamePrefix+'_junctions.list',dtype=int)
    
    # loop through basins
    for key, Basin in np.ndenumerate(Basins):

        # print basin to screen
        print(key[0], Basin)
        
        # isolate basin data
        BasinChannelData = ChannelsDF[ChannelsDF.basin_key == key]
        MinimumChi = BasinChannelData.chi.min()
        MaximumMChi = BasinChannelData.m_chi.max()
        MinKsn = BasinChannelData.m_chi.min()
        MaxKsn = BasinChannelData.m_chi.max()
        
        # how many segments are we dealing with?    
        Segments = BasinChannelData.segment_number.unique()
        
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
                ax1.scatter(Chi,KKsn,marker='o', edgecolors='k',lw=0.5, c=[0.8,0.8,0.8], s=20, zorder=20)
            else:    
                ax1.plot(Chi,Elevation,'-', lw=1.5,c=ColourMap(Colour), zorder=10)

            # get hillslope data
            SegmentHillslopeData = HillslopesDF[HillslopesDF.StreamID == Segment]
            NTraces = len(SegmentHillslopeData)

            if NTraces<20:
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
            ax1.set_ylim(0,PlotMaxKsn)
        
        ax2.set_ylim(0,PlotDF.EStarUpper.max())
        ax3.set_ylim(0,1)

        #save output
        plt.suptitle('Basin ID ' + str(key[0]) + " (" + str(Basin) + ")")
        
        if PlotKsn:
            plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(key[0]) + "_ChiProfileEsRsKsn.png", dpi=300)
        else:    
            plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(key[0]) + "_ChiProfileEsRs.png", dpi=300)
        plt.clf()
        plt.close()    
    
