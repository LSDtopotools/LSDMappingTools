# -*- coding: utf-8 -*-
"""
Created on Mon May 29 15:07:06 2017

@author: smudd
"""

import numpy as np
import LSDPlottingTools as LSDP
import matplotlib.pyplot as plt

def ChiMOverNTest(start_movern = 0.1, d_movern = 0.1, n_movern = 6):

    #DataDirectory = "T:\\analysis_for_papers\\movern_testing\\"
    DataDirectory = "C:\\VagrantBoxes\\LSDTopoTools\\Topographic_projects\\Meghalaya_chi_test\\"
    #movern_profile_file = "Irian_Jaya_PP_movern.csv"
    movern_profile_file = "Mega_divide_movern.csv"
    #Base_file = "Mega_clip"

    label_size = 10

    # Set up fonts for plots
    #rcParams['font.family'] = 'sans-serif'
    #rcParams['font.sans-serif'] = ['arial']
    #rcParams['font.size'] = label_size
    size_format = "default"
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
    basin = thisPointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    source = thisPointData.QueryData('source_key')
    source = [int(x) for x in source]
    
    # need to convert everything into arrays so we can mask different basins
    Elevation = np.asarray(elevation)
    Basin = np.asarray(basin)
    Source = np.asarray(source)
    
    print(Elevation)
    
    for idx,mn in enumerate(m_over_n_values):
        # first get the chi values for this m_over_n
        print ("mn is: " + str(mn))
        print("index is: "+str(idx))
        mn_legend = "m_over_n = "+str(mn)
        print("I am looking for the data element: "+mn_legend)
        this_chi = thisPointData.QueryData(mn_legend)
        
        # convert to a numpy array for masking
        Chi = np.asarray(this_chi)
        
        # some info about the chi and elevation values
        max_chi = np.amax(Chi)
        max_Elevation = np.amax(Elevation)
        min_Elevation = np.amin(Elevation)        

        z_axis_min = int(min_Elevation/10)*10
        z_axis_max = int(max_Elevation/10)*10+10

        chi_axis_max = int(max_chi/5)*5+5
        
        # Now mask the data. Initially we will do only basin 0
        basin_key = 0
        if basin_key == 0:      # We dont use this but I am putting conditional statement here so we can have proper indent
            # this gets the mask (for the chosen basin)        
            m = np.ma.masked_where(Basin!=basin_key, Basin)
            
            # this is the masked chi value
            maskX = np.ma.masked_where(np.ma.getmask(m), Chi)
            maskElevation = np.ma.masked_where(np.ma.getmask(m), Elevation)
            
            # now plot the data with a colourmap
            ax.scatter(maskX,maskElevation,s=2.0, c=maskElevation,cmap="cubehelix",edgecolors='none')

            # some formatting of the figure
            ax.spines['top'].set_linewidth(1)
            ax.spines['left'].set_linewidth(1)
            ax.spines['right'].set_linewidth(1)
            ax.spines['bottom'].set_linewidth(1)

            # make the lables
            ax.set_xlabel("$\chi$")
            ax.set_ylabel("Elevation (m)")

            # This affects all axes because we set share_all = True.
            ax.set_ylim(z_axis_min,z_axis_max)
            ax.set_xlim(0,chi_axis_max) 
            
            #save the plot
            newFilename = DataDirectory+"test.png"
                
            # This gets all the ticks, and pads them away from the axis so that the corners don't overlap    
            ax.tick_params(axis='both', width=1, pad = 2)
            for tick in ax.xaxis.get_major_ticks():
                tick.set_pad(2)

            FigFormat = "png"
            plt.savefig(newFilename,format=FigFormat,dpi=500)
            fig.clf()
         
    

if __name__ == "__main__":
    ChiMOverNTest()