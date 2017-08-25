"""This class contains calculation for sorting data in order to make the plotting scripts lighter
@author Boris Gailleton"""
import numpy as np
import pandas as bamboo_bears
from pyproj import Proj, transform

def source_num_of_longest_river(arr):
    """
    This function will return the source number of the longuest channel by counting the number of occurence of each sources number.
    It may be an issue if your channel have similar size, I will create a fancier function that calculate the actual lenght of each channel if required.
    """
    counount = 0
    max_count = -999
    max_count_index = -999

    for i in range(1,arr.shape[0]):
        if(arr[i] == arr[i-1]):
            counount += 1
        else:
            if (counount > max_count):
                max_count_index = arr[i-1]
                max_count = counount
            counount = 0

    return max_count_index

def crop_point_data_to_base_raster(raster_name, raster_directory, csv_file, EPSG_code = 0):
    """
    This function create a new csv file cropped to the base raster. It can lower the processing time if your point data is on a significantly larger area than the base raster.
    """
    print("ok let me load your dataset and your hdr file")
    # Read the file
    df = bamboo_bears.read_csv(csv_file)
    # Read and sort the csv_info
    with open(raster_directory+raster_name+".hdr","r") as hdr_file:
        print("I got these")
        for line in hdr_file:
            if(line[0:8] == "map info"):
                info = line[12:-2]
                info = info.split(",")
                x_min = float(info[3])
                y_max = float(info[4])
                x_res = float(info[5])
                y_res = float(info[6])
                utm_zone = int(info[7])
                utm_hemisphere = info[8]
            else:
                if(line[0:7] == "samples"):
                    num_col = line.replace(" ","").split("=")[1]
                    print("there are " + str(num_col) + " columns")
                    num_col = int(num_col)
                else:
                    if(line[0:5] == "lines"):
                        num_lines = line.replace(" ","").split("=")[1]
                        print("there are " + str(num_lines) + " lines")
                        num_lines = int(num_lines)

    # Now I calculate the size of the dem
    x_max = x_min + x_res*num_col
    y_min = y_max - y_res*num_lines
    # Conversion UTM to lat/long
    inProj = Proj(init='epsg:'+str(EPSG_code))
    outProj = Proj(init='epsg:4326')
    long_min,lat_min = transform(inProj,outProj,x_min,y_min)
    long_max,lat_max = transform(inProj,outProj,x_max,y_max)

    # data sorting
    df = df[df.longitude<long_max]
    df = df[df.latitude<lat_max]
    df = df[df.latitude>lat_min]
    df = df[df.longitude>long_min]
    df.to_csv(csv_file[:-4]+"_"+raster_name+"_filtered.csv", index = False)

    #return the name of the new csv file
    return csv_file[:-4]+"_"+raster_name+"_filtered.csv"

def select_data_from_one_col(csv_file, col_name = "None",values =[], wanted = True):

    """
    This function selects specific data from a csv file using pandas (the package, not the majestic animal).
    It requires the name+path of csv, the name of the column, a python list containing the wanted values and the wanted parameter (True to select the data, False to select the opposite)
    """
    # Ignore this
    word = "excluding"
    if(wanted):
        word = "including"
    ###########
    print("I am going to sort your data " + col_name +" by " + word + " the following value(s):")
    print(values)
    # import the file
    df = bamboo_bears.read_csv(csv_file,sep = ",")
    # Select the data
    if(wanted):
        retArray = df[df[col_name].isin(values)]
    else:
        if(wanted == False):
            retArray = df[~df[col_name].isin(values)]

    retArray.to_csv(csv_file[:-4]+"_SDFOC.csv", index = False)
    print("Your new file is ready and is there: " + str(csv_file[:-4]+"_SDFOC.csv"))
    return csv_file[:-4]+"_SDFOC.csv"

def quick_load_bil_raster(raster_name, raster_path):
    """
    This function return a numpy 2Darray from a bil-ENVI raster. It uses numpy.fromfile(), quite efficient.
    """
    # Read the number of raw/col to recast the array
    with open(raster_path+raster_name+".hdr","r") as hdr_file:
        print("I am opening your raster info")
        for line in hdr_file:
            #testing the data type
            if(line[0:9] == "data type"):
                info_dtype = line[-2:]
                info_dtype = int(info_dtype)
            else:
                if(line[0:8] == "map info"):
                    info = line[12:-2]
                    info = info.split(",")
                    x_min = float(info[3])
                    y_max = float(info[4])
                    x_res = int(info[5])
                    y_res = int(info[6])
                    utm_zone = int(info[7])
                    utm_hemisphere = info[8]
                else:
                    if(line[0:7] == "samples"):
                        num_col = line.replace(" ","").split("=")[1]
                        print("there are " + str(num_col) + " columns")
                        num_col = int(num_col)
                    else:
                        if(line[0:5] == "lines"):
                            num_lines = line.replace(" ","").split("=")[1]
                            print("there are " + str(num_lines) + " lines")
                            num_lines = int(num_lines)

        #The type of data representation:
        #1 = Byte: 8-bit unsigned integer
        #2 = Integer: 16-bit signed integer
        #3 = Long: 32-bit signed integer
        #4 = Floating-point: 32-bit single-precision
        #5 = Double-precision: 64-bit double-precision floating-point
        #6 = Complex: Real-imaginary pair of single-precision floating-point
        #9 = Double-precision complex: Real-imaginary pair of double precision floating-point
        #12 = Unsigned integer: 16-bit
        #13 = Unsigned long integer: 32-bit
        #14 = 64-bit long integer (signed)
        #15 = 64-bit unsigned long integer (unsigned)

    if(info_dtype == 1):
        data_type = np.dtype('uint8')
    else:
        if(info_dtype == 2):
            data_type = np.dtype('int16')
        else:
            if(info_dtype == 3):
                data_type = np.dtype('int32')
            else:
                if(info_dtype == 4):
                    data_type = np.dtype('Float32')
                else:
                    if(info_dtype == 5):
                        data_type = np.dtype('Float64')
                    else:
                        if(info_dtype == 12):
                            data_type = np.dtype('uint16')
    print("your data type is "  + str(data_type))
    # Alright now loading and converting the data
    print("I am now ingesting your raster")
    raster_data = np.fromfile(raster_path+raster_name+".bil",data_type).reshape(num_lines,num_col)
    print("I nailed it")
    x_max = x_min + x_res*num_col
    y_min = y_max - y_res*num_lines
    print("I am returning the raster array and info")
    return raster_data,x_min,x_max,y_min,y_max,x_res,y_res

def get_basin_middle(df,wdir = "", key = "basin_key", keyDA = "drainage area"):
    """
    Create a dataframe containing the lat/long of the center of a basin and the key of this one. Usefull to visualize basin key when QGIS is not working...
    Default values adapted for output from chi_mapping tool
    Author: Boris
    """

    dfbasin = []

    lst_basins = df[key].unique()
    nb_basin = lst_basins.shape[0]
    for bas in lst_basins:
        print("I am processing basin %s/%s " %(bas, nb_basin))
        temp = []
        long_max = df["longitude"][df[key] == bas].max()
        long_min = df["longitude"][df[key] == bas].min()
        lat_max = df["latitude"][df[key] == bas].max()
        lat_min = df["latitude"][df[key] == bas].min()
        DA = df[keyDA][df[key] == bas].max()
        BK = bas
        temp = [(lat_min+lat_max)/2,(long_min+long_max)/2,BK,DA]
        dfbasin.append(temp)

    header = ["latitude", "longitude", "basin_key", "drainage area"]
    print("I am now creating a dataframe")
    dfbasin = bamboo_bears.DataFrame(dfbasin,columns = header)
    print("I am saving a csv file f you want to reload directly this data as %sBASIN_LOCA.csv" %(wdir))
    dfbasin.to_csv(wdir+"BASIN_LOCA.csv", index = False)
    print(dfbasin)
    print("done with the basin location")
    return dfbasin

def pandas_to_unique_df(idf,columns = 'source_ID', otype = 'list'):
    """
        Take a dataframe, return a list or dictionary of dataframes from the uniques values of one of its colums.
        Can be useful to plot basins_river more easily, I already needed it in some cases, I don't really remember why!
        But now it is coded anyway.
        @args:
            idf (pandas dataframe): The input dataframe
            columns (str): the name of the column containing the unique values
            type (str): Do you want a list or a dict (the key would be the unique value)

        @return:
            A list of dataframes
    """
    if otype == 'list':
        odf = [] # the Output list
        for i in idf.unique(): # looping through the unique values of the dataframe
            odf.append(idf[idf[columns] == i]) # appending a single dataframe

    else:
        odf = {}
        for i in idf.unique():
            odf[i] = idf[idf[columns] == i]
    return odf










    #
