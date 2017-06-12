"""
A set of functions to do some simple statistical analyses
Created on Thu Jun 8th 2017

    Author: SMM
"""
from __future__ import absolute_import, division, print_function, unicode_literals


import numpy as np
import pandas as pd

# This function comes from
# https://github.com/joferkington/oost_paper_code/blob/master/utilities.py
def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh

def add_outlier_column_to_PD(df, column = "none", threshold = "none"):

    """
    Takes a pandas dataframe and returns the same with added boolean columns (True if outlier).
    Uses the function is_outlier to detect the outliers. Can also take a list of dataframes
    Args:
        df (Pandas dataframe): The dataframe or a list of dataframe
        column (list or string): name of the column(s) you want to outlier-check
        threshold  (list or float): list of threshold for each columns (must be same size as column)
    returns
        Pandas.DataFrame
    """
    print("DEBUG")

    if(isinstance(df,list) == False and isinstance(df,dict) == False):
        lst_df = [df]

    else:
        lst_df = df

    if(isinstance(column,str) and column =="none"):
        print("you need to give me the name of at least a column, or a list ([])")
        quit()
    if(isinstance(threshold,str) and threshold =="none"):
        print("you need to give me the name of at least a column, or a list ([])")
        quit()

    for instance in lst_df:

        if(isinstance(column,str)):
            column = [column]

        if(isinstance(threshold,float) or isinstance(threshold,int)):
            threshold = [threshold]

        if(len(threshold) != len(column)):
            print("You need to assign one threshold per columns name")

        for i in range(len(column)):
            # print(lst_df[instance])
            is_outliers = is_outlier(lst_df[instance][column[i]],threshold[i])
            coln =column[i]+"_outlier"
            lst_df[instance][coln] = pd.Series(is_outliers,index = lst_df[instance].index)

    if len(lst_df) == 1:
        return lst_df[0]
    else:
        return lst_df

def binning_PD(df, column = "", values = [], log = False):
    """
    takes a dataframe (Pandas) and return a list of dataframes binned by one columns.
    Args:
        df: The pandas dataframe
        column (str): name of the column that hold the data
        values (list): _ list of the upper values of each binning, another binning category will incorporate everything superior to the last value
                       _ you also can give the value "auto_power_10" to this. it will automatically bin the data each 10**n until the max
                       _ "unique" will bin the df for each values of the column (Basin/ source key for example)
        log (bool): if you want to compare values to log of column
    return:
        dictionnary of pandas dataframe, the key being the upper value
    """
    # check the function parameters
    if(column == ""):
        print("You need to give a valid column name")
    if(isinstance(values, list) and len(values) < 2):
        print("You need at least two values to bin the dataframe")
    if(isinstance(values,str) and values == "auto_power_10"):
        print("I am automatically choosing the binning values each 10**n, thanks for trusting me")
        max_val = df[column].max()
        min_value = df[column].min()

        po = 0
        values = []

        while(max_val>10**po):
            if(min_value<10**po):
                values.append(10**po)
            po +=1

        del values[-1] # delete the last value to keep last bin > to last value
        print("Your binning values are: ")
        print(values)
    else:
        if(isinstance(values,str) and values == "unique"):
            print("I am automatically choosing the binning values for each unique values, thanks for trusting me")
            values = df[column].unique()
            print("Your binning values are: ")
            print(values)
    cumul_lines = 0# check if all the values are inside bins


    # log the data if required
    if(log):
        return_DF = [df[np.log10(df[column])<values[0]]]
        cumul_lines += return_DF[0].shape[0]
        for i in range(1,len(values)):
            tempdf = df[np.log10(df[column])<values[i]]
            tempdf = tempdf[tempdf[column]>values[i-1]]
            return_DF.append(tempdf)
            cumul_lines += return_DF[i].shape[0]
        tempdf= df[np.log10(df[column])>=values[-1]]
        cumul_lines += tempdf.shape[0]
        return_DF.append(tempdf)
    else:
        return_DF = [df[df[column]<values[0]]]
        cumul_lines += return_DF[0].shape[0]
        for i in range(1,len(values)):
            tempdf = df[df[column]<values[i]]
            tempdf = tempdf[tempdf[column]>values[i-1]]
            return_DF.append(tempdf)
            cumul_lines += return_DF[i].shape[0]

        cumul_lines += df[df[column]>=values[-1]].shape[0]
        return_DF.append(df[df[column]>=values[-1]]) # last overweight bin

        print("DEBUG: " +str(cumul_lines) +" lines detected over " +str(df.shape[0]))
    # compile the results in a dictionnary
    dict_return = {}
    for i in range(len(values)):
        dict_return[str(values[i])] = return_DF[i]

    dict_return['>'+str(values[-1])] = return_DF[-1]

    return dict_return

















#
