"""
A set of functions to do some simple statistical analyses
Created on Thu Jun 8th 2017

    Author: SMM
"""
from __future__ import absolute_import, division, print_function, unicode_literals


import numpy as np
import pandas as pd
import scipy.stats as ss

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

    # make sure that you don't get a divide by zero.
    # If MAD is 0, then there are no outliers
    if med_abs_deviation == 0:
        modified_z_score = diff * 0
    else:
        modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh


def get_MAD(points):
    """
    Returns the MAD of a ampled population:
    -----------
        points : An numobservations by numdimensions array of observations
        
    Returns:
    --------
        MAD.

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

   
    return med_abs_deviation


def get_outlier_from_KernelDensityStuff(df, column = "", binning = "", threshold = 6, method = "gaussian", sort_by = ""):

    from sklearn.neighbors.kde import KernelDensity as harry


    if(column ==""):
        print("I need a column")
        quit()





    out_df = pd.DataFrame(data = None, columns = df.columns)

    if(binning != ""):
        for potter in df[binning].unique():
            tdf = df[df[binning] == potter]
            if(sort_by!= ""):
                tdf = tdf.sort_values(sort_by)

            if(column == "deriv_ksn"):
                tdf["deriv_ksn"] = pd.Series(np.abs(tdf.ksn.diff()/tdf.chi.diff()),index = tdf.index)
                tdf["deriv_ksn"].iloc[0] = 0
                print(tdf.shape[0])

            dumbledore = np.copy(tdf[column].values.reshape((-1,1)))

            severus = harry(kernel = method).fit(dumbledore)

            snake = np.abs(severus.score_samples(dumbledore))

            McGonagal = []
            for gobelins in snake:
                if gobelins<threshold:
                    McGonagal.append(False)
                else:
                    McGonagal.append(True)

            aCat = tdf[McGonagal]


            out_df = pd.concat([out_df,aCat])


    return out_df





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

    # Check the DataType
    if(isinstance(df,list) == False and isinstance(df,dict) == False):
        lst_df = [df]

    else:
        lst_df = df
    # check the data validity
    if(isinstance(column,str) and column =="none"):
        print("you need to give me the name of at least a column, or a list ([])")
        quit()
    if(isinstance(threshold,str) and threshold =="none"):
        print("you need to give me the name of at least a column, or a list ([])")
        quit()



    # calculate the outliers
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

def dixon_test(data, left=True, right=True, q_dict = ""):
    """
    Keyword arguments:
        data = A ordered or unordered list of data points (int or float).
        left = Q-test of minimum value in the ordered list if True.
        right = Q-test of maximum value in the ordered list if True.
        q_dict = A dictionary of Q-values for a given confidence level,
            where the dict. keys are sample sizes N, and the associated values
            are the corresponding critical Q values. E.g.,
            {3: 0.97, 4: 0.829, 5: 0.71, 6: 0.625, ...}

    Returns a list of 2 values for the outliers, or None.
    E.g.,
       for [1,1,1] -> [None, None]
       for [5,1,1] -> [None, 5]
       for [5,1,5] -> [1, None]

    """
    assert(left or right), 'At least one of the variables, `left` or `right`, must be True.'
    assert(len(data) >= 3), 'At least 3 data points are required'
    assert(len(data) <= max(q_dict.keys())), 'Sample size too large'

    sdata = sorted(data)
    Q_mindiff, Q_maxdiff = (0,0), (0,0)

    if left:
        Q_min = (sdata[1] - sdata[0])
        try:
            Q_min /= (sdata[-1] - sdata[0])
        except ZeroDivisionError:
            pass
        Q_mindiff = (Q_min - q_dict[len(data)], sdata[0])

    if right:
        Q_max = abs((sdata[-2] - sdata[-1]))
        try:
            Q_max /= abs((sdata[0] - sdata[-1]))
        except ZeroDivisionError:
            pass
        Q_maxdiff = (Q_max - q_dict[len(data)], sdata[-1])

    if not Q_mindiff[0] > 0 and not Q_maxdiff[0] > 0:
        outliers = [None, None]

    elif Q_mindiff[0] == Q_maxdiff[0]:
        outliers = [Q_mindiff[1], Q_maxdiff[1]]

    elif Q_mindiff[0] > Q_maxdiff[0]:
        outliers = [Q_mindiff[1], None]

    else:
        outliers = [None, Q_maxdiff[1]]

    return outliers

def linregress_residuals(xdata,ydata):
    """
    This function performs a linear regression and then gets the residuals
    
    Args:
        xdata (array-like): The x data
        ydata (array-like): The y data
        
    Returns:
        residuals: the residuals of the regression
        slope: the slope of regression line
        intercept: intercept of the regression line
        rvalue: correlation coeffficient
        pvalue: two-sided p-value for a hypothesis test whose null hypothesis is that the slope is zero.
        stderr: standard error of the estimated gradient
    
    Author: SMM
    """
    
    from scipy import stats
    
    # Get the regression
    (m,b,r,pvalue,stderr)=stats.linregress(xdata,ydata)   
    #print("In regress1, m is: "+str(m))    
    # get the residuals
    residuals = np.copy(xdata)
    for idx,x in enumerate(xdata):
        yfit = m*x+b
        residuals[idx] = yfit-ydata[idx]
        
    return (residuals,m,b,r,pvalue,stderr)

def remove_outlying_residuals(xdata,ydata,residuals):
    """
    This function removes data with outying residuals
    
    Args:
        xdata (array-like): The x data
        ydata (array-like): The y data    
        residuals: The residuals

        
    Returns:
        new_x: data with outlier removed
        new_y: data with outliers removed
        is_outlier_vec: A vec of bools denoiting if they are outliers
        m: the slope of the regression
        b: the intercept of the regression
            
    Author: SMM
    """
    from scipy import stats

    # get the outliers
    is_outlier_vec = is_outlier(residuals, thresh=3.5)
    
    # now remove the outliers from the dataset
    new_x = []
    new_y = []
    for idx,x in enumerate(xdata):
        if not is_outlier_vec[idx]:
            new_x.append(xdata[idx])
            new_y.append(ydata[idx])
  
    NX = np.asarray(new_x)
    NY = np.asarray(new_y)
    
    # now get the new regression
    (m,b,r,pvalue,stderr)=stats.linregress(NX,NY)

    return (NX,NY, is_outlier_vec, m,b)
    
    
    


def extract_outliers_by_header(df, data_column_name = "diff", header_for_group = "source_key", threshold = 3.5):
    """
    Extract outliers from a dataframe, groupped by a specific column. 
    for example, as in default, extract the outliers in the diff column, groupped by source_key

    param:
        df (pandas dataframe): The dataframe containing at least the two columns in param
        data_column_name (str): Name of the column containing the data
        header_for_group (str): Name of the column for groupping

    return: a dataframe containing only the outliers

    Author: BG - 05/10/2017

    """

    ### Let's do it

    # Output dataframe initialized as an integer, you'll see why later 
    out_df = 0
    for i in df[header_for_group].unique():
        # I am grouping the dataset for each value of the header, for example each basins
        tdf = df[df[header_for_group]==i]
        tdfpositive = tdf[tdf["sign"] == 1]
        tdfnegative = tdf[tdf["sign"] == -1]

        # masking the outliers
        maskpositive = is_outlier(tdfpositive[data_column_name].values, thresh = threshold)
        masknegative = is_outlier(tdfnegative[data_column_name].values, thresh = threshold)
        # Applying the mask
        tdfpositive = tdfpositive[maskpositive]
        tdfnegative = tdfnegative[masknegative]
        tdf = pd.concat([tdfpositive,tdfnegative])
        # Feeding the out dataset
        
        if(isinstance(out_df, int)):
            out_df = tdf
        else:
            out_df = pd.concat([out_df,tdf])
    print("I selected %s outliers" %(out_df.shape[0]))

    return out_df


















#
