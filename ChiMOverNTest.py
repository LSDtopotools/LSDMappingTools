# -*- coding: utf-8 -*-
"""
Created on Mon May 29 15:07:06 2017

@author: smudd
"""

import numpy as np
import LSDPlottingTools as LSDP
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd

def ChiMOverNTest(start_movern = 0.1, d_movern = 0.1, n_movern = 6):

    DataDirectory = "T:\\analysis_for_papers\\movern_testing\\"
    #DataDirectory = "C:\\VagrantBoxes\\LSDTopoTools\\Topographic_projects\\Meghalaya_chi_test\\"
    movern_profile_file = "Irian_Jaya_PP_movern.csv"
    movern_basin_stats_file = "Irian_Jaya_PP_movernstats_basinstats.csv"
    #movern_profile_file = "Mega_divide_movern.csv"
    #Base_file = "Mega_clip"
 

    end_movern = start_movern+d_movern*(n_movern-1)
    m_over_n_values = np.linspace(start_movern,end_movern,n_movern)    
    
    
    # get the maximum MLE of each basin
    pd_DF = pd.DataFrame.from_csv(DataDirectory+movern_basin_stats_file)
    shp = pd_DF.shape
    max_MLEs = []
    max_MLEs_index = []
    for i in range(0,shp[0]):
        #print("I is: "+str(i))
        a = pd_DF.loc[[i]]
        b = np.asarray(a)
        c = b[0,:]
        max_MLEs.append(max(c))
        max_MLEs_index.append(np.argmax(c))
    print("max_MLEs are: ")
    print(max_MLEs)
    m_over_n_of_max = []
    
    for idx in max_MLEs_index:
        m_over_n_of_max.append(m_over_n_values[idx])
           
    print("The m over n of these max are: ")    
    print(m_over_n_of_max)
    
    n_basins = len(max_MLEs)
        
    

    label_size = 12

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size
    size_format = "default"
    # make a figure,
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        #l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        #l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
        #l_pad = -35

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax = fig.add_subplot(gs[5:100,10:95])


    


    # load the m_over_n data file
    thisPointData = LSDP.LSDMap_PointData(DataDirectory+movern_profile_file)
    allBasinStatsData = LSDP.LSDMap_PointData(DataDirectory+movern_basin_stats_file)     

    
    print("m over n values are: ")
    print(m_over_n_values)
    
    mn_legends = []
    for mn in m_over_n_values:
        mn_legends.append("m_over_n = "+str(mn))
        
    print("The mn labels are: ")
    print(mn_legends)
        

    # get the data form the profiles
    elevation = thisPointData.QueryData('elevation')
    elevation = [float(x) for x in elevation]
    basin = thisPointData.QueryData('basin_key')
    basin = [int(x) for x in basin]
    source = thisPointData.QueryData('source_key')
    source = [int(x) for x in source]
    
    # get the basin keys
    allstats_basinkeys = allBasinStatsData.QueryData("basin_key")
    allstats_basinkeys = [int(x) for x in allstats_basinkeys]
    
    # need to convert everything into arrays so we can mask different basins
    Elevation = np.asarray(elevation)
    Basin = np.asarray(basin)
    Source = np.asarray(source)
    
    #print(Elevation)
    
    # Loop through m/n values aggregating data
    for idx,mn in enumerate(m_over_n_values):
        
        counter = str(idx).zfill(3)
        print("Counter is: "+counter)
        
        # first get the chi values for this m_over_n
        #print ("mn is: " + str(mn))
        #print("index is: "+str(idx))
        mn_legend = "m_over_n = "+str(mn)
        #print("I am looking for the data element: "+mn_legend)
        this_chi = thisPointData.QueryData(mn_legend)
        
        # get the MLE value for this m/n
        this_MLE = allBasinStatsData.QueryData(mn_legend)
        
        
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
        basin_key = 12
        if(basin_key == 12):                  # We dont use this but I am putting conditional statement here so we can have proper indent
        #for basin_key in range(0,n_basins-1):      
            
            # now we need to find out if this basin is in the allstats file, 
            # and if so what index it is
            #this_basin_index = -99
            #for ii,bk in enumerate(allstats_basinkeys):
            #    #print("index: "+str(ii)+" and basin_key is: "+str(bk))
            #    if (bk == basin_key):
            #        this_basin_index = idx
                    
            #if(this_basin_index != -99):        
            #    MLE = this_MLE[this_basin_index]
            #else:
            #    MLE = "NaN"
            MLE = this_MLE[basin_key] 
            #MLE_str = str(MLE)
            short_MLE = str("%03.02e" % round(MLE,2))
            print("The short MLE is: "+short_MLE)
            
            #print("The MLE of this basin for this m over n is: "+short_MLE)
                        
        
        
            # this gets the mask (for the chosen basin)        
            m = np.ma.masked_where(Basin!=basin_key, Basin)
            
            # this is the masked chi value
            maskX = np.ma.masked_where(np.ma.getmask(m), Chi)
            maskElevation = np.ma.masked_where(np.ma.getmask(m), Elevation)
            
            # now plot the data with a colourmap
            ax.scatter(maskX,maskElevation,s=2.5, c=maskElevation,cmap="terrain",edgecolors='none')

            # some formatting of the figure
            ax.spines['top'].set_linewidth(1)
            ax.spines['left'].set_linewidth(1)
            ax.spines['right'].set_linewidth(1)
            ax.spines['bottom'].set_linewidth(1)

            # make the lables
            ax.set_xlabel("$\chi$ (m)")
            ax.set_ylabel("Elevation (m)")

            # This affects all axes because we set share_all = True.
            #ax.set_ylim(z_axis_min,z_axis_max)
            #ax.set_ylim(z_axis_min,z_axis_max)
            #ax.set_xlim(0,chi_axis_max)
            #plt.title("Basin = " +mn_legend+", MLE = "+short_MLE)
            
            #newline = "\n"
            title_string = "Basin "+str(basin_key)+", $m/n$ = "+str(mn)
            title_string2 = "MLE = "+short_MLE
            ax.text(0.05, 0.95, title_string,
                    verticalalignment='top', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='black', fontsize=10)
            print("The basin index is: "+str(basin_key)+" and the max index is: "+str(max_MLEs_index[basin_key]))
            if( idx == max_MLEs_index[basin_key]):
                print("This m/n is: "+str(mn)+" and it is the maximum MLE")
                ax.text(0.05, 0.88, title_string2+", maximum MLE in basin.",
                    verticalalignment='top', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='red', fontsize=10)
            else:
                ax.text(0.05, 0.88, title_string2,
                    verticalalignment='top', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='black', fontsize=10)                
            
            #save the plot
            newFilename = DataDirectory+"Chi_profiles_basin_"+str(basin_key)+"_"+counter+".png"
                
            # This gets all the ticks, and pads them away from the axis so that the corners don't overlap    
            ax.tick_params(axis='both', width=1, pad = 2)
            for tick in ax.xaxis.get_major_ticks():
                tick.set_pad(2)

            FigFormat = "png"
            plt.savefig(newFilename,format=FigFormat,dpi=300)
            ax.cla()
            #plt.show()
         
    

if __name__ == "__main__":
    ChiMOverNTest()