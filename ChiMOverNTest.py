# -*- coding: utf-8 -*-
"""
Created on Mon May 29 15:07:06 2017

@author: smudd
"""

import numpy as np
import LSDPlottingTools as LSDP
import matplotlib.pyplot as plt

def ChiMOverNTest(start_movern = 0.1, d_movern = 0.1, n_movern = 6):

    label_size = 10

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure,
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[25:100,10:95])


    
    DataDirectory = "T:\\analysis_for_papers\\movern_testing\\"
    movern_profile_file = "Irian_Jaya_PP_movern.csv"
    #Base_file = "Mega_clip"

    # load the m_over_n data file
    thisPointData = LSDP.LSDMap_PointData(DataDirectory+movern_profile_file)    

    end_movern = start_movern+d_movern*(n_movern-1)
    m_over_n_values = np.linspace(start_movern,end_movern,n_movern)
    
    print("m over n values are: ")
    print(m_over_n_values)
    
    mn_legends = []
    for mn in m_over_n_values:
        mn_legends.append("m_over_n = "+str(mn))
        
    print("The mn labels are: ")
    print(mn_legends)
        


    elevation = thisPointData.QueryData('elevation')
    elevation = [float(x) for x in elevation]
    #m_chi = thisPointData.QueryData('m_chi')
    #m_chi = [float(x) for x in m_chi]
    basin = thisPointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    source = thisPointData.QueryData('source_key')
    source = [int(x) for x in source]
    
    # loop through m over n
    
    for mn,idx in enumerate(m_over_n_values):
        
        # First we need to mask the data
    

if __name__ == "__main__":
    ChiMOverNTest()