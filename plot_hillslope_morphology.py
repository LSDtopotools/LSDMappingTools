# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#import modules
from __future__ import print_function
from geopandas import GeoDataFrame
from shapely.geometry import LineString, shape, Point, MultiPolygon, Polygon
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
    rcParams['font.sans-serif'] = ['Liberation Sans']
    rcParams['font.size'] = 8
    rcParams['text.usetex'] = False
        
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

def ProcessSegmentedData(DataDirectory, FilenamePrefix):
    """
    This function reads channel and hillslope data and organises it by basins and 
    segments, taking median values for each segment and 16/84 percentiles as range estimates
    
    Args:
        DataDirectory: the data directory
        FilenamePrefix: the file name prefix

    Author: MDH
    """

    #load the hillslopes data
    HillslopeData = ReadHillslopeData(DataDirectory, FilenamePrefix)
    
    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)
    
    # Get a list of unique basins to loop through
    Basins = ChannelData.basin_key.unique()
    
    #loop through the basins
    for Basin in Basins:
        
        # isolate basin data
        BasinChannelData = ChannelData[ChannelData.basin_key == Basin]
        BasinJunctions = HillslopeData.BasinID.unique()
        BasinHillslopeData = HillslopeData[HillslopeData.BasinID == BasinJunctions[BasinID]]
        
        # how many segments are we dealing with?    
        Segments = BasinChannelData.segment_number.unique()
    

    
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

    
def MapBasinChannelHillslopes(BasinID):
    
    """
    Makes a plot of the basin, this must be the first script to run since the basin sets the plot extent
    
    MDH, September 2017
    
    """
    
    #import module for converting hillslope data to lat long
    from pyproj import Proj, transform
    
    # Open basins shapefile, convert to latlong, and select basin by ID
    PolygonDict, CRS = rmt.ReadShapeFile(ShapeFile)
    
    # get a sorted list of keys for selecting the basin    
    Keys = np.sort(np.array(PolygonDict.keys()))
    
    # select the polygon of the basin we're interested in
    BasinPoly = dict((key,value) for key, value in PolygonDict.iteritems() if key == Keys[BasinID])
    
    # create a shapefile of the selected basin
    BasinShapeFile = DataDirectory+FilenamePrefix+"_basin_"+str(BasinID)+".shp"
    HillshadeRaster = Directory+FilenamePrefix+"_hs_resample_latlong.bil"
    
    # define the output coordinate system
    rmt.WriteShapeFile(BasinPoly,CRS,BasinShapeFile)
    
    # Create a map figure with basin shapefile as extent
    Fig, Ax, Map = rmt.CreateMapFigure(BasinShapeFile)
    
    # plot the hillshade as base layer
    rmt.PlotRaster(HillshadeRaster,Map)
    rmt.PlotShapefile(BasinShapeFile,Map,Ax,'k','w')
    
    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)
    
    # isolate basin data
    BasinChannelData = ChannelData[ChannelData.basin_key == BasinID]
    lon = BasinChannelData.longitude.as_matrix()
    lat = BasinChannelData.latitude.as_matrix()
    x,y = Map(lon,lat)
    
    # channels marked by chi
    Chi = BasinChannelData.chi.as_matrix()
    Chi -= np.min(Chi)
    MarkerSize = 2.-(Chi/np.max(Chi))*2.
    plt.scatter(x,y,marker='.',color='b',s=MarkerSize)
    
    #load the hillslopes data
    HillslopeData = ReadHillslopeData(DataDirectory, FilenamePrefix)
    BasinJunctions = HillslopeData.BasinID.unique()
    BasinHillslopeData = HillslopeData[HillslopeData.BasinID == Keys[BasinID]]
    X = BasinHillslopeData.X.as_matrix()
    Y = BasinHillslopeData.Y.as_matrix()
    
    #Convert to Lat Long
    Output_CRS = Proj({'init': "epsg:4326"})
    lon, lat = transform(CRS, Output_CRS, X, Y)
    x,y = Map(lon,lat)
    
    plt.scatter(x,y,marker='.',color='r',s=2)
    
    plt.savefig(PlotDirectory+FilenamePrefix+"_basin_plot.png", dpi=300)
    plt.close()
    
    
def PlotChiElevationSegments(BasinID):
    
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
    plt.close()
    
def PlotLongProfileSegments(BasinID):
    
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
    plt.close()
    
def PlotChiElevationMChi(BasinID):
    
    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)

    if (BasinChannelData.count == 0):
        print("No Channel Data for Basin ID " + str(BasinID))

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
    plt.close()
    
def PlotLongProfileMChi(BasinID):
    
    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)

    # isolate basin data
    BasinChannelData = ChannelData[ChannelData.basin_key == BasinID]
    
    if (BasinChannelData.count == 0):
        print("No Channel Data for Basin ID " + str(BasinID))
        
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
    plt.close()
    
def PlotLongProfileMChiCht(BasinID):
    
    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)

    #load the hillslopes data
    HillslopeData = ReadHillslopeData(DataDirectory, FilenamePrefix)
    
    # isolate basin data
    BasinChannelData = ChannelData[ChannelData.basin_key == BasinID]

    BasinJunctions = HillslopeData.BasinID.unique()
    
    BasinHillslopeData = HillslopeData[HillslopeData.BasinID == BasinJunctions[BasinID]]
        
    MinimumDistance = BasinChannelData.flow_distance.min()
    MaximumMChi = BasinChannelData.m_chi.max()
    
    # how many segments are we dealing with?    
    Segments = BasinChannelData.segment_number.unique()
    
    # setup the figure
    Fig = CreateFigure(FigSizeFormat="JGR")
    Ax = plt.subplot(111)
    
    #choose colormap
    ColourMap = cm.viridis
    
    #empty lists for hilltop data
    ChtMedian=np.zeros(len(Segments))
    Cht25=np.zeros(len(Segments))
    Cht75=np.zeros(len(Segments))
    Distances=np.zeros(len(Segments))
    MinDistances=np.zeros(len(Segments))
    MaxDistances=np.zeros(len(Segments))
    NTraces=np.zeros(len(Segments))
    
    # Get the data columns for plotting
    for i in range(0, len(Segments)):
        
        #get data arrays
        Dist = ChannelData.flow_distance[ChannelData.segment_number == Segments[i]]
        Elevation = ChannelData.elevation[ChannelData.segment_number == Segments[i]]
        SegmentedElevation = ChannelData.segmented_elevation[ChannelData.segment_number == Segments[i]]
        MChi = ChannelData.m_chi[ChannelData.segment_number == Segments[i]].unique()[0]
        
        # get hillslope data
        SegmentHillslopeData = BasinHillslopeData[BasinHillslopeData.StreamID == Segments[i]]
        ChtMedian[i] = SegmentHillslopeData.Cht.quantile(0.5)
        Cht25[i] = SegmentHillslopeData.Cht.quantile(0.25)
        Cht75[i] = SegmentHillslopeData.Cht.quantile(0.75)
        NTraces[i] = SegmentHillslopeData.size
            
        #normalise distance by outlet distance
        Dist = Dist-MinimumDistance
        Distances[i] = Dist.mean()
        MinDistances[i] = Dist.min()
        MaxDistances[i] = Dist.max()
        
        #plot, colouring segments
        Colour = MChi/MaximumMChi
        plt.plot(Dist/1000,Elevation,'k--',dashes=(2,2), lw=0.5,zorder=10)
        plt.plot(Dist/1000, SegmentedElevation, '-', lw=2, c=ColourMap(Colour),zorder=9)
    
    # Finalise the figure
    plt.xlabel('Distance (km)')
    plt.ylabel('Elevation (m)')
    
    
    #plot the hillslope data
    Ax2 = plt.twinx()
    Ax.set_zorder(Ax2.get_zorder()+1) # put ax in front of ax2 
    Ax.patch.set_visible(False) # hide the 'canvas' 

    print("plotting errorbar")
    plt.errorbar(Distances/1000.,ChtMedian,yerr=np.asarray([Cht75,Cht25]),fmt='s',color=[0.6,0.6,0.6],ms=4,elinewidth=1,zorder=-32)
    plt.ylabel("$C_{HT}$ (m$^{-1}$)")
    
    #add colourbar
    CAx = Fig.add_axes([0.64,0.32,0.2,0.02])
    m = cm.ScalarMappable(cmap=ColourMap)
    m.set_array(ChannelData.m_chi)
    plt.colorbar(m, cax=CAx,orientation='horizontal')
    plt.xlabel('$M_{\chi}$ m$^{0.64}$',fontsize=8)
    CAx.tick_params(axis='both', labelsize=8)
    
    #save output
    plt.suptitle('Basin ID ' + str(BasinID))
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(BasinID) + "_LongProfMChiCht.png", dpi=300)
    plt.close()

def PlotMChiCht(BasinID):
    
    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)
    
    #load the hillslopes data
    HillslopeData = ReadHillslopeData(DataDirectory, FilenamePrefix)
        
    # isolate basin data
    BasinChannelData = ChannelData[ChannelData.basin_key == BasinID]
    BasinJunctions = np.sort(HillslopeData.BasinID.unique())
    BasinHillslopeData = HillslopeData[HillslopeData.BasinID == BasinJunctions[BasinID]]
    
    # segments in the hillslope data
    #Segments = BasinHillslopeData.StreamID.unique()
    Segments = BasinChannelData.segment_number.unique()
    
    # setup the figure
    Fig = CreateFigure()
    plt.subplot(111)
    
    #choose colormap
    ColourMap = cm.viridis
    
    # For each segment get the MChi value and collect dimensionless hillslope data
    # record the number of traces in each segment inbto a new dataframe
    Data = pd.DataFrame(columns=['SegmentNo','MChi','FlowLength','SegmentLength','ChtMedian','ChtLower','ChtUpper','NTraces'])
    
    for i in range(0, len(Segments)):
        
        #Get segment hillslope data
        SegmentHillslopeData = BasinHillslopeData[BasinHillslopeData.StreamID == float(Segments[i])]
        
        #Get segment channel data and calculate flow length
        SegmentChannelData = BasinChannelData[BasinChannelData.segment_number == Segments[i]]
        
        #channels
        MChi = SegmentChannelData.m_chi.unique()[0]
        TempFL = SegmentChannelData.flow_distance
        FlowLength = np.median(TempFL)
        SegmentLength = np.max(TempFL)-np.min(TempFL)
        
        #hillslopes
        ChtMedian = SegmentHillslopeData.Cht.quantile(0.5)
        ChtLower = SegmentHillslopeData.Cht.quantile(0.25)
        ChtUpper = SegmentHillslopeData.Cht.quantile(0.75)
        NTraces = SegmentHillslopeData.size
        
        #add to data frame
        Data.loc[i] = [Segments[i],MChi,FlowLength,SegmentLength,ChtMedian,ChtLower,ChtUpper,NTraces]
  
    # remove rows with no data (i.e. no hillslope traces)
    Data = Data.dropna()
    Data = Data[Data.NTraces > 50]
    
    # colour code by flow length
    MinFlowLength = Data.FlowLength.min()
    Data.FlowLength = Data.FlowLength-MinFlowLength
    MaxFlowLength = Data.FlowLength.max()
    colours = (Data.FlowLength/MaxFlowLength)
    
    # Error bars with colours but faded (alpha)
    for i, row in Data.iterrows(): 
        ChtErr = np.array([[row.ChtLower],[row.ChtUpper]])
        plt.plot([row.MChi,row.MChi],ChtErr,'-', lw=1.5, color=ColourMap(colours[i]), alpha=0.5,zorder=9)
        plt.plot(row.MChi,row.ChtMedian,'s',color=ColourMap(colours[i]),zorder=32)
        
    
    # Finalise the figure
    plt.xlabel('$M_{\chi}$ m$^{0.64}$')
    plt.ylabel('$C_{HT}$ m$^{-1}$')
    plt.suptitle('Basin ID ' + str(BasinID))
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    # add colour bar
    m = cm.ScalarMappable(cmap=ColourMap)
    m.set_array(FlowLength)
    cbar = plt.colorbar(m)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.set_label('Flow Length (m)')
    
    #save output
    plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(BasinID) + "_MChi_Cht.png", dpi=300)
    plt.close()    
    
def PlotMChiEstar(BasinID,Sc=1.):
    
    # load the channel data
    ChannelData = ReadChannelData(DataDirectory, FilenamePrefix)
    
    #load the hillslopes data
    HillslopeData = ReadHillslopeData(DataDirectory, FilenamePrefix)
    
    # isolate basin data
    BasinChannelData = ChannelData[ChannelData.basin_key == BasinID]
    BasinJunctions = HillslopeData.BasinID.unique()
    BasinHillslopeData = HillslopeData[HillslopeData.BasinID == BasinJunctions[BasinID]]
    MinimumDistance = BasinChannelData.flow_distance.min()
    MaximumMChi = BasinChannelData.m_chi.max()
    
    # how many segments are we dealing with?    
    Segments = BasinChannelData.segment_number.unique()
    
    # setup the figure
    Fig = CreateFigure(AspectRatio=1.)
    Ax = plt.subplot(111)
    
    #choose colormap
    ColourMap = cm.viridis
    
    # Calculate E_Star
    BasinHillslopeData.E_Star = 2.*BasinHillslopeData.Cht*BasinHillslopeData.Lh/Sc
    
    # For each segment get the MChi value and collect dimensionless hillslope data
    # record the number of traces in each segment
    MChi = []
    EstarMean = []
    EstarStD = []
    EstarStE = []
    NTraces = []
    for i in range(0, len(Segments)):
        SegmentHillslopeData = BasinHillslopeData[BasinHillslopeData.StreamID == Segments[i]]
        if SegmentHillslopeData.size != 0:
            MChi.append(ChannelData.m_chi[ChannelData.segment_number == Segments[i]].unique()[0])
            EstarMean.append(SegmentHillslopeData.E_Star.mean())
            EstarStD.append(SegmentHillslopeData.E_Star.std())
            EstarStE.append(SegmentHillslopeData.E_Star.std()/np.sqrt(float(SegmentHillslopeData.size)))
            NTraces.append(SegmentHillslopeData.size)
    
    #make the plot and colour code by number of hillslope traces
    NTraces = np.asarray(NTraces,dtype=float)
    [plt.errorbar(MChi[i],EstarMean[i],yerr=EstarStD[i],fmt='s',elinewidth=1.5,ecolor=ColourMap(NTraces[i]/np.max(NTraces)),mfc=ColourMap(NTraces[i]/np.max(NTraces))) for i in range(0,len(MChi))]
    
    # Finalise the figure
    plt.xlabel('$M_{\chi}$ m$^{0.64}$')
    plt.ylabel('$E*$')
    plt.text(-0.2,-0.3,'Basin ID ' + str(BasinID),transform = Ax.transAxes,color=[0.35,0.35,0.35])
    plt.tight_layout()
    
    # add colour bar
    m = cm.ScalarMappable(cmap=ColourMap)
    m.set_array(FlowLength)
    cbar = plt.colorbar(m)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.set_label('Flow Length (m)')
    
    #save output
    plt.savefig(PlotDirectory+FilenamePrefix + "_" + str(BasinID) + "_MChiEstar.png", dpi=300)
    plt.close()

def CalculateRStar(EStar):
    """
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

def CalculateEStarRStar(Basin,Sc=0.71):

    """

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
    Data = Data[Data.NTraces > 50]
    
    return Data

def PlotEStarRStar(Basin, Sc=0.71):
    """
    MDH
    """
    
    Data = CalculateEStarRStar(Basin)
    
    # setup the figure
    Fig = CreateFigure(AspectRatio=1.2)
        
    #choose colormap
    ColourMap = cm.viridis

    #Plot analytical relationship
    PlotEStarRStarTheoretical()
    
    # colour code by flow length
    MinFlowLength = Data.FlowLength.min()
    Data.FlowLength = Data.FlowLength-MinFlowLength
    MaxFlowLength = Data.FlowLength.max()
    colours = (Data.FlowLength/MaxFlowLength)
    
    #plot the data
    plt.loglog()
    
    # Error bars with colours but faded (alpha)
    for i, row in Data.iterrows(): 
        EStarErr = np.array([[row.EStarLower],[row.EStarUpper]])
        RStarErr = np.array([[row.RStarLower],[row.RStarUpper]])
        plt.plot([row.EStar,row.EStar],RStarErr,'-', lw=1, color=ColourMap(colours[i]), alpha=0.5,zorder=9)
        plt.plot(EStarErr,[row.RStar,row.RStar],'-', lw=1, color=ColourMap(colours[i]), alpha=0.5,zorder=9)
        plt.plot(row.EStar,row.RStar,'o',ms=4,color=ColourMap(colours[i]),zorder=32)

    # Finalise the figure
    plt.xlabel('$E^*={{-2\:C_{HT}\:L_H}/{S_C}}$')
    plt.ylabel('$R^*=S/S_C$')
    plt.xlim(0.1,1000)
    plt.ylim(0.01,1.5)
        
    # add colour bar
    m = cm.ScalarMappable(cmap=ColourMap)
    m.set_array(Data.FlowLength)
    cbar = plt.colorbar(m)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.set_label('Distance to Outlet (m)')
    
    plt.suptitle("Basin "+str(Basin)+" Dimensionless Hillslope Morphology")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(PlotDirectory+FilenamePrefix + "_" + "%02d" % Basin + "_EStarRStar.png", dpi=300)
    plt.close(Fig)
    
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

def PlotOrthogonalResiduals(ModelX, ModelY, DataX, DataY):

    """
    """
    
    # setup the figure
    Fig = CreateFigure(AspectRatio=1.2)
    
    # Get residuals
    Residuals, OrthoX, OrthoY = OrthogonalResiduals(ModelX, ModelY, DataX, DataY)
    
    # plot model and data
    plt.loglog()
    plt.axis('equal')
    plt.plot(ModelX,ModelY,'k-', lw=1)
    plt.plot(DataX,DataY,'k.', ms=2)
    
    # plot orthogonals
    for i in range(0,len(DataX)):
        plt.plot([DataX[i],OrthoX[i]],[DataY[i],OrthoY[i]],'-',color=[0.5,0.5,0.5])
    
    plt.savefig(PlotDirectory+FilenamePrefix + "_ESRSOrthoResiduals.png", dpi=300)
    
def OrthogonalResiduals(ModelX, ModelY, DataX, DataY):

    """
    Orthogonal residuals between data and model using dot product approach
    
    MDH, September 2017
    
    """
    
    
    # Start search at beginning of model vector
    StartSearch = 0
    EndSearch = len(ModelX)

    #holder for residuals and points on the model line
    Residuals = np.zeros(len(DataX))
    XOrtho = np.zeros(len(DataX)) 
    YOrtho = np.zeros(len(DataX)) 
    
    # loop through the data and do the orthogonal thing    
    for i in range(0,len(DataX)):
        
        # pull out X and Y
        X = DataX[i]
        Y = DataY[i]
        
        # fractional distance along the model line (0-1 if orthogonal)        
        t = 0
        
        # loop across model to search for orthogonal        
        for j in range(StartSearch,EndSearch-1):
            
            # Get model vector dx and dy
            dX12 = ModelX[j+1]-ModelX[j]
            dY12 = ModelY[j+1]-ModelY[j]
            
            # Get Data to Model Vector dx and dy
            dX13 = X-ModelX[j]
            dY13 = Y-ModelY[j]
            
            # Calculate the Dot Product
            DotProduct = dX12*dX13 + dY12*dY13
            
            # Calculate the fraction of distance along the line
            lastt = t;
            t = DotProduct/(dX12*dX12 + dY12*dY12)
            
            # Check we haven't passed the solution
            if ((lastt < 0. and t > 1.) or (lastt > 1. and t < 0)):
                
                #solution is a model node                
                XOrtho[i] = ModelX[j]
                YOrtho[i] = ModelY[j]
            
            # If we haven't reached the solution then continue along the model
            elif (t<0. or t>1.):
                continue
            
            # Otherwise find the solution interpolated along the model vector
            else:
                
                #Find point along line
			XOrtho[i] = ModelX[j] + t*dX12
			YOrtho[i] = ModelY[j] + t*dY12
   
            # Now calculate the residual and track if it is positive or negative
            dX = XOrtho[i]-X
            dY = YOrtho[i]-Y
            Residuals[i] = np.sqrt(dX**2. + dY**2.)*(dY/np.abs(dY))
            break

    return Residuals, XOrtho, YOrtho
    
def EStarRStarResiduals(Sc=0.71):

    """
    MDH
    
    """

    # setup the figure
    Fig = CreateFigure(FigSizeFormat="EPSL",AspectRatio=4.)
    Ax = plt.subplot(111)
    
    # Calculate analytical relationship
    EStar = np.logspace(-1,3,1000)
    RStar = CalculateRStar(EStar)
    
    plt.plot([-6,41],[0,0],'k--',lw=0.5,zorder=1)
    plt.fill([-6,41,41,-6],[0,0,0.5,0.5],color=[1.0,0.8,0.8],zorder=0)
    plt.fill([-6,41,41,-6],[0,0,-0.5,-0.5],color=[0.8,0.9,1.0],zorder=0)
    plt.text(-5,0.35,"Growing")
    plt.text(-5,-0.45,"Decaying")
    
    #loop through the basins
    for Basin in range(0,NoBasins):
    #for Basin in range(0,1):

        # Get the hillslope data for the basin        
        Data = CalculateEStarRStar(Basin)
        
        # Calculate log transformed orthogonal residuals
        Residuals, Xortho, Yortho = OrthogonalResiduals(np.log10(EStar), np.log10(RStar), np.log10(Data.EStar.as_matrix()), np.log10(Data.RStar.as_matrix()))

        # plot violin of residuals
        violin_parts = plt.violinplot(-Residuals, [Basin,],showmeans=False, showmedians=False, showextrema=False)
        
        # set the colour of the distribution
        for fc in violin_parts['bodies']:
            fc.set_facecolor([0.2,0.2,0.2])
        
        # plot the median
        Median = -np.median(Residuals)
        plt.plot([Basin-0.2,Basin+0.2],[Median,Median],'k-',lw=1)
            
        
        
    plt.xlabel("Basin Number")
    plt.ylabel("Orthogonal Residuals\n($log_{10}\:E* R*$)")
    plt.xlim(-6,41)
    plt.ylim(-0.5,0.5)
    plt.xticks(np.arange(0,NoBasins,2))
    Ax.tick_params(axis='x',which='minor',bottom='on')
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(PlotDirectory+FilenamePrefix + "_ESRSOrthoResidualsViolin.png", dpi=300)
    
def PlotHillslopeHistograms(Basin):
    """
    Plots histograms of Hillslope Length, Hillslope Gradient, Relief and Hilltop Curvature
    
    MDH, September 2017
    
    """

    # load the hillslopes data and isolate segments in basin
    HillslopeData = ReadHillslopeData(DataDirectory, FilenamePrefix)
    Basins = np.sort(HillslopeData.BasinID.unique())
    BasinHillslopeData = HillslopeData[HillslopeData.BasinID == Basins[Basin]]
    
    # hillslope length
    CreateFigure(FigSizeFormat="EPSL", AspectRatio=1.)
    
    plt.subplot(221)
    freq, bin_edges = np.histogram(BasinHillslopeData.Lh,bins=np.arange(0,400,20.))
    freq_norm = freq.astype(np.float)/float(np.max(freq))
    plt.bar(bin_edges[:-1],freq_norm,width=20,align='edge',edgecolor='k',linewidth=0.5,color='steelblue')
    plt.xlabel("Hillslope Length (m)")
    plt.ylabel("Normalised Frequency")
    plt.xlim(0,400)
    
    # hillslope gradient
    ax2 = plt.subplot(222)
    freq, bin_edges = np.histogram(BasinHillslopeData.S,bins=np.arange(0,0.8,0.04))
    freq_norm = freq.astype(np.float)/float(np.max(freq))
    plt.bar(bin_edges[:-1],freq_norm,width=0.04,align='edge',edgecolor='k',linewidth=0.5,color='thistle')
    plt.xlabel("Mean Slope (m/m)")
    plt.xlim(0,0.8)
    ax2.yaxis.set_ticklabels([])
    
    # hilltop curvature
    plt.subplot(223)
    freq, bin_edges = np.histogram(BasinHillslopeData.R,bins=np.arange(0,200.,10.))
    freq_norm = freq.astype(np.float)/float(np.max(freq))
    plt.bar(bin_edges[:-1],freq_norm,width=10.,align='edge',edgecolor='k',linewidth=0.5,color='sandybrown')
    plt.xlabel("Hillslope Relief (m)")
    plt.ylabel("Normalised Frequency")
    plt.xlim(0,200)
    
    
    # hillslope relief
    ax4 = plt.subplot(224)
    freq, bin_edges = np.histogram(BasinHillslopeData.Cht,bins=np.arange(-0.1,0,0.005))
    freq_norm = freq.astype(np.float)/float(np.max(freq))
    plt.bar(bin_edges[:-1],freq_norm,width=0.005,align='edge',edgecolor='k',linewidth=0.5,color='salmon')
    plt.xlabel("Hilltop Curvature (m$^{-1}$)")
    plt.xlim(0,-0.1)
    ax4.yaxis.set_ticklabels([])
    
    plt.suptitle("Basin "+str(Basin)+" Hillslopes")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(PlotDirectory+FilenamePrefix + "_" + "%02d" % Basin + "_HillslopeHists.png", dpi=300)
    plt.close()

def JoyPlot(HillslopeData,Column,XLabel,Colour,Outfile,BinMin,BinSpacing,BinMax):
    
    CreateFigure(AspectRatio=0.5,FigSizeFormat="small")
    Ax = plt.subplot(111)
    
    Basins = np.sort(HillslopeData.BasinID.unique())
    
    for Basin in range(0,NoBasins):
        #Get hillslopes for basin
        BasinHillslopeData = HillslopeData[HillslopeData.BasinID == Basins[Basin]]

        #create the PDF
        freq, BinEdges = np.histogram(BasinHillslopeData[Column],bins=np.arange(BinMin,BinMax+BinSpacing,BinSpacing))
        BinMidpoints = BinEdges+BinSpacing*0.5
        freq_norm = freq.astype(np.float)/float(np.max(freq))

        #plot, offset by Basin #
        plt.plot(BinMidpoints[:-1],freq_norm-Basin,'k-',linewidth=1)
        plt.fill_between(BinMidpoints[:-1],freq_norm-Basin,-Basin,color=Colour)
    
    if np.abs(BinMin) < np.abs(BinMax):
        plt.xlim(BinMin,BinMax)
    else:
        plt.xlim(BinMax,BinMin)
        BinSpacing *= -1
        
    plt.xlabel(XLabel)
    plt.text(-BinSpacing,0,"North-West",rotation=90,verticalalignment='top')
    plt.text(-BinSpacing,-(NoBasins-1),"South-East",rotation=90,verticalalignment='bottom')
    
    #only display bottom axis
    Ax.spines['right'].set_visible(False)
    Ax.spines['top'].set_visible(False)
    Ax.spines['left'].set_visible(False)
    Ax.yaxis.set_visible(False)
    
    plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.98])
    plt.savefig(PlotDirectory+Outfile, dpi=300)
    plt.clf()
    
def PlotHillslopeJoyPlots():
    
    """
    Make a joy plot of slope for all basins
    
    MDH
    
    """    
    
    # load the hillslopes data and isolate segments in basin
    HillslopeData = ReadHillslopeData(DataDirectory, FilenamePrefix)
    
    #Hillslope Length
    BinMin = 0; BinSpacing = 20.; BinMax = 400
    JoyPlot(HillslopeData,"Lh","Hillslope Length (m)","steelblue","Length_Joy.png",BinMin,BinSpacing,BinMax)
    
    #Slope
    BinMin = 0; BinSpacing = 0.04; BinMax = 0.8
    JoyPlot(HillslopeData,"S","Slope (m/m)","thistle","Slope_Joy.png",BinMin,BinSpacing,BinMax)
    
    #Hillslope Relief
    BinMin = 0; BinSpacing = 10.; BinMax = 200.
    JoyPlot(HillslopeData,"R","Hillslope Relief (m)","sandybrown","Relief_Joy.png",BinMin,BinSpacing,BinMax)
    
    #Hilltop Curv
    BinMin = -0.1; BinSpacing = 0.005; BinMax = 0.
    JoyPlot(HillslopeData,"Cht","Hilltop Curvature (m$^{-1}$)","Salmon","HilltopCurv_Joy.png",BinMin,BinSpacing,BinMax)

def DetermineSc():
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
    
def StitchImages2Animation(InputFile, OutputFile):
    
    """
    This is reliant on mencoder

    MDH
    
    """    
    # import the subprocess module for launching external commands
    import subprocess, os
    
    # create the file list
    f = open("filelist.txt","w")
    for Basin in range(0,NoBasins):
        f.write(FilenamePrefix + "_" + "%02d" % Basin + InputFile + "\n")
    
    os.chdir(PlotDirectory)
    args = ["mencoder", "mf://@filelist.txt", "-mf", "fps=2:type=png", "-ovc", "copy", "-oac", "copy", "-o", OutputFile]
    subprocess.call(args)
        
if __name__ == "__main__":
    #do something
    print("Running...")

    # linux
    if platform == "linux" or platform == "linux2":    
        Directory = "/home/mhurst/bolinas_paper/"
        DataDirectory = Directory+"data/"+"hillslope_data/"
        PlotDirectory = Directory+"plots/"
    
    # windows
    elif platform == "win32":
        Directory = "C:\Users\\Martin Hurst\\Dropbox\\Glasgow\\Projects\\Bolinas\\analysis\\bolinas\\"
        DataDirectory = Directory+"data\\"+"hillslope_data\\"
        PlotDirectory = Directory+"plots\\"
    
    else:
        print("Platform not recognised")
        
    #file names and extensions
    RasterExtension = "bil"
    FilenamePrefix = "bolinas"
    #ShapeFile = Directory+FilenamePrefix+"_AllBasins.shp"
    #NoBasins = 41
    
    #ModelX = np.arange(0.,100.)
    #ModelY = np.arange(0.,100.)
    #DataX = np.arange(10.,90.)
    #DataY = DataX + 5.*np.random.randn(len(DataX))
    #PlotOrthogonalResiduals(ModelX,ModelY,DataX,DataY)
    #EStarRStarResiduals()
    #PlotEStarRStarProgression()
    #PlotEStarRStarResiduals()
    #SaveHillslopeDataByBasin(DataDirectory,FilenamePrefix)    
    #SaveChannelDataByBasin(DataDirectory,FilenamePrefix)
    
    #print("Basin is ", end='')
    #Transience = np.zeros(NoBasins)
    #for Basin in range(0,NoBasins):
        #print(Basin, end=',')
        #stdout.flush()
        #PlotMChiCht(Basin)
        #Transience[Basin] = PlotEStarRStar(Basin)
        #PlotHillslopeHistograms(Basin)
    
    #print(Transience)
    
    #CreateFigure()
    #plt.plot(np.arange(0,NoBasins),Transience,'ko')
    #plt.plot([0,NoBasins],[0,0],'k--')
    #plt.fill([0,NoBasins,NoBasins,0],[0,0,1,1],color=[1.0,0.9,0.9])
    #plt.fill([0,NoBasins,NoBasins,0],[0,0,-1,-1],color=[0.9,0.9,1.0])
    #plt.xlabel("Basin No")
    #plt.ylabel("Percentage of points above Steady-state")
    #plt.tight_layout()
    #plt.savefig(PlotDirectory+"Transience.png",dpi=300)
    #PlotChiElevationMChi(0)
    #PlotLongProfileMChiCht(0)
    #PlotMChiCht(9)
    #PlotHillshadeBasins(Rotation=45.)
    #PlotEStarRStar(12)
    #PlotHillslopeHistograms(0)
    
    #MapBasinChannelHillslopes(0)
    
    #make movies
    #InputFile = "_HillslopeHists.png"
    #OutputFile = PlotDirectory+"HillslopeHists.avi"
    #StitchImages2Animation(InputFile, OutputFile)
    
    #PlotHillslopeJoyPlots()
    
    #DetermineSc()
    
    WriteHillslopeTracesShp(DataDirectory,FilenamePrefix)
    
    print("Done")
