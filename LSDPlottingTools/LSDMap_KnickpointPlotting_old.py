## LSDMap_KnickpointPlotting.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## These functions are tools for analysing and plotting knickpoint data
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## FJC
## 05/06/2017
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import LSDPlottingTools.LSDMap_PointTools as LSDMap_PD
import matplotlib.pyplot as plt
import time as clock
from matplotlib import rcParams
import matplotlib.cm as cm
from LSDMapFigure.PlottingRaster import MapFigure
from LSDMapFigure.PlottingRaster import BaseRaster
from LSDMapFigure import PlottingHelpers as Helper
from LSDPlottingTools import colours as lsdcolours
from LSDPlottingTools import statsutilities as SUT
from LSDPlottingTools import init_plotting_DV
import LSDPlottingTools as LSDP
import sys
import os
import pandas as pd
from scipy.stats import norm
import utm
pd.options.mode.chained_assignment = None



class KP_plotting(object):
    """
        This class is a development version of the knickpoint algorithm. 
        Its aim is to deal with knickpoint picking stuffs via python code, 
        in a development aim as I am significantly changing methods relatively often.
        B.G.
    """


    def __init__(self, fpath,fprefix, basin_key = [], source_key = [], min_length = 0):

        print("Let me first preprocess and check your files")

        # Loading the attributes
        self.fpath = fpath # the path of your file : /home/your/path/
        self.fprefix = fprefix # the common prefix of all your files

        # Loading the files

        print("Loading the knickpoint-related files")
        
        try:
            self.df_river = Helper.ReadMChiSegCSV(self.fpath, self.fprefix, type = "knickpoint")
            self.df_kp_raw = Helper.ReadKnickpointCSV(self.fpath, self.fprefix, ftype = "raw")
            self.df_kp = Helper.ReadKnickpointCSV(self.fpath, self.fprefix)
            self.df_SK = Helper.readSKKPstats(self.fpath, self.fprefix)
        except IOError:
            print("I didnae find your knickpoint related files make sure that:")
            print("- Your path is good and finishing by '/' e.g. /home/name/kazakhstan/ ")
            quit()


        print("Managing the data:")
        if(basin_key == []):
            print("All the basins are selected:")
            print(self.df_SK["basin_key"].unique().tolist())
        else:
            print("You selected the following basins:")
            print(basin_key)
            self.df_river = self.df_river[self.df_river["basin_key"].isin(basin_key)]
            self.df_kp_raw = self.df_kp_raw[self.df_kp_raw["basin_key"].isin(basin_key)]
            self.df_kp = self.df_kp[self.df_kp["basin_key"].isin(basin_key)]
            self.df_SK = self.df_SK[self.df_SK["basin_key"].isin(basin_key)]

        if(source_key == [] and min_length == 0):
            print("All the sources are selected:")
            print(self.df_SK["source_key"].unique().tolist())
        elif(min_length > 0):
            print("Let me remove the river smaller than " +str(min_length))
            self.df_SK = self.df_SK[self.df_SK["length"]>min_length]
            source_key = self.df_SK["source_key"].unique()
            self.df_river = self.df_river[self.df_river["source_key"].isin(source_key)]
            self.df_kp_raw = self.df_kp_raw[self.df_kp_raw["source_key"].isin(source_key)]
            self.df_kp = self.df_kp[self.df_kp["source_key"].isin(source_key)]
            self.df_SK = self.df_SK[self.df_SK["source_key"].isin(source_key)]
        else:
            print("You selected the following Sources: ")
            print(source_key)
            self.df_river = self.df_river[self.df_river["source_key"].isin(source_key)]
            self.df_kp_raw = self.df_kp_raw[self.df_kp_raw["source_key"].isin(source_key)]
            self.df_kp = self.df_kp[self.df_kp["source_key"].isin(source_key)]
            self.df_SK = self.df_SK[self.df_SK["source_key"].isin(source_key)]

        # TODO Adding some sorting functions based on the source key

        print("Done now")



    def DEBUG_print_ksn_filters(self):
        """
            This function is used to print one ksn profile per river to check the effect of the different filters on the dataset
            BG - 12/01/2018
        """
        plt.clf()
        print("I will now print ksn(chi) with the different filter")
        svdir = self.fpath+'river_plots/'
        if not os.path.isdir(svdir):
            os.makedirs(svdir)

        for SK in self.df_river["source_key"].unique():
            print("printing river: " +str(SK))

            # Selecting the river
            df = self.df_river[self.df_river["source_key"] == SK]

            fig = plt.figure(1, facecolor='white',figsize=(9,5))

            gs = plt.GridSpec(100,100,bottom=0.10,left=0.10,right=0.95,top=0.95)
            ax1 = fig.add_subplot(gs[0:100,0:100])

            ax1.scatter(df["chi"], df["m_chi"], c = "r", s = 1, marker = "o", label = "ksn")
            ax1.scatter(df["chi"], df["lumped_ksn"], c = "g", s = 1, marker = "s", label = "lumped ksn")
            ax1.scatter(df["chi"], df["TVD_ksn"], c = "k", s = 1, marker = "+", label = "TVD ksn")

            ax1.legend()

            ax1.set_xlabel(r'$ \chi$')
            ax1.set_ylabel(r'$ k_{sn}$')

            plt.savefig(svdir + self.fprefix + "_ksn_SK_" +str(SK)+".png", dpi = 300)
            plt.clf()

    def DEBUG_print_ksn_outliers(self):
        """
            This function is used to print one ksn profile per river to check the effect of the different filters on the dataset
            BG - 12/01/2018
        """
        plt.clf()
        print("I will now print ksn(chi) with outliers")
        svdir = self.fpath+'river_plots/'
        if not os.path.isdir(svdir):
            os.makedirs(svdir)

        for SK in self.df_river["source_key"].unique():
            print("printing river: " +str(SK))

            # Selecting the river
            df = self.df_river[self.df_river["source_key"] == SK]
            dfo = self.df_kp[self.df_kp["source_key"] == SK]

            fig = plt.figure(1, facecolor='white',figsize=(9,5))

            gs = plt.GridSpec(100,100,bottom=0.10,left=0.10,right=0.95,top=0.95)
            ax1 = fig.add_subplot(gs[0:100,0:100])

            ax1.scatter(df["chi"], df["m_chi"], c = "r", s = 1, marker = "o", label = "ksn", alpha = 0.15)
            # ax1.scatter(df["chi"], df["lumped_ksn"], c = "g", s = 1, marker = "s", label = "lumped ksn")
            # ax1.scatter(df["chi"], df["TVD_ksn_NC"], c = "purple", s = 2, marker = "x", label = "TVD ksn non corrected")
            ax1.scatter(df["chi"], df["TVD_ksn"], c = "k", s = 1, marker = "+", label = "TVD ksn")

            ax1.scatter(dfo["chi"][dfo["out"]==1], dfo["delta_ksn"][dfo["out"]==1], c = "purple" , marker = "s", s = 4)
            lim = ax1.get_ylim()
            ax1.vlines(dfo["chi"][dfo["out"]==1],-1000,1000, lw = 0.5, alpha = 0.5)
            ax1.set_ylim(lim)


            ax1.legend()

            ax1.set_xlabel(r'$ \chi$')
            ax1.set_ylabel(r'$ k_{sn}$')

            plt.savefig(svdir + self.fprefix + "_outksn_SK_" +str(SK)+".png", dpi = 300)
            plt.clf()

    def DEBUG_print_ksn_dksndchi(self):
        """
            This function is used to print one ksn profile per river to check the effect of the different filters on the dataset
            BG - 12/01/2018
        """
        plt.clf()
        print("I will now print ksn(chi) with outliers")
        svdir = self.fpath+'river_plots/'
        if not os.path.isdir(svdir):
            os.makedirs(svdir)

        for SK in self.df_river["source_key"].unique():
            print("printing river: " +str(SK))

            # Selecting the river
            df = self.df_river[self.df_river["source_key"] == SK]
            dfo = self.df_kp_raw[self.df_kp_raw["source_key"] == SK]

            fig = plt.figure(1, facecolor='white',figsize=(9,5))

            gs = plt.GridSpec(100,100,bottom=0.10,left=0.10,right=0.90,top=0.95)
            ax1 = fig.add_subplot(gs[0:100,0:100])
            ax2 = fig.add_subplot(gs[0:100,0:100], facecolor = "None")



            ax1.scatter(df["chi"], df["m_chi"], c = "gray", s = 1, marker = "o", label = "ksn")
            # ax1.scatter(df["chi"], df["lumped_ksn"], c = "g", s = 1, marker = "s", label = "lumped ksn")
            
            ax1.scatter(df["chi"], df["TVD_ksn"], c = "k", s = 1, marker = "+", label = "TVD ksn")
            ylim = ax1.get_ylim()
            ax2.scatter(dfo["chi"], dfo["delta_ksn"], c = "r", s = 3, marker = "s", label = r'$ \frac{d(TVD_ksn)}{d\chi}$')

            ax2.yaxis.set_label_position('right')
            ax2.yaxis.set_ticks_position('right')
            ax2.xaxis.set_visible(False)
            ax1.yaxis.set_label_position('left')
            ax1.yaxis.set_ticks_position('left')



            ax2.set_xlim(ax1.get_xlim())
            ax1.set_ylim(ylim)
            # ax1.scatter(dfo["chi"][dfo["out_MZS"]==1], dfo["delta_ksn"][dfo["out_MZS"]==1], c = "r" , marker = "s", s = 2)


            # ax1.legend()

            ax1.set_xlabel(r'$ \chi$')
            ax1.set_ylabel(r'$ k_{sn}$')
            ax2.set_ylabel(r'$ \frac{d(TVD_ksn)}{d\chi}$')

            plt.savefig(svdir + self.fprefix + "_ksn_rawkp_SK_" +str(SK)+".png", dpi = 300)
            # plt.show()
            plt.clf()



    def DEBUG_print_KDE(self):
        """
            This function is used to print one ksn profile per river to check the effect of the different filters on the dataset
            BG - 12/01/2018
        """
        plt.clf()
        print("I will now print ksn(chi) with the different KDE")
        svdir = self.fpath+'river_plots/'
        if not os.path.isdir(svdir):
            os.makedirs(svdir)

        for SK in self.df_kp_raw["source_key"].unique():
            print("printing river: " +str(SK))

            # Selecting the river
            df = self.df_kp_raw[self.df_kp_raw["source_key"] == SK]

            fig = plt.figure(1, facecolor='white',figsize=(9,5))

            gs = plt.GridSpec(100,100,bottom=0.10,left=0.10,right=0.95,top=0.95)
            ax1 = fig.add_subplot(gs[0:100,0:100])

            ax1.scatter(df["dksn/dchi"], df["KDE"], c = "k", s = 1, marker = "+", label = "ksn")

            ax1.set_xlabel(r'$ \frac{dk_{sn}}{\chi}$')
            ax1.set_ylabel(r'$ KDE_pdf $')

            plt.savefig(svdir + self.fprefix + "_KDE_SK_" +str(SK)+".png", dpi = 300)
            plt.clf()

    def DEBUG_print_ksn_comparison(self):
        """
            This function is used to print one ksn profile per river to check the effect of the different filters on the dataset
            BG - 12/01/2018
        """
        plt.clf()
        print("I will now print ksn(chi) with the different filter")
        svdir = self.fpath+'river_plots/'
        if not os.path.isdir(svdir):
            os.makedirs(svdir)

        for SK in self.df_river["source_key"].unique():
            print("printing river: " +str(SK))

            # Selecting the river
            df = self.df_river[self.df_river["source_key"] == SK]

            df_kp = self.df_kp[self.df_kp["source_key"] == SK]
            df_kp_raw = self.df_kp_raw[self.df_kp_raw["source_key"] == SK]
            fig = plt.figure(1, facecolor='white',figsize=(9,5))

            gs = plt.GridSpec(100,100,bottom=0.10,left=0.10,right=0.95,top=0.95)
            ax1 = fig.add_subplot(gs[0:100,0:100])

            ax1.scatter(df_kp["chi"][df_kp["out"] == 1], df_kp["delta_ksn"][df_kp["out"] == 1], c = "purple",s=5 , marker = "s", label = "kp final outliers")
            ax1.scatter(df["chi"], df["TVD_ksn"], c = "k", s = 2, marker = "o", label = "ksn", alpha = 0.5)
            ax1.scatter(df_kp_raw["chi"], df_kp_raw["delta_ksn"], c = "green" , marker = "x", label = "kp raw", alpha = 0.5)
            ax1.scatter(df_kp["chi"], df_kp["delta_ksn"], c = "red" , marker = "+", label = "kp final", alpha = 0.5)

            ax1.legend()

            ax1.set_xlabel(r'$ \chi$')
            ax1.set_ylabel(r'$ k_{sn}$')

            plt.savefig(svdir + self.fprefix + "_DEBUG_kp_SK_" +str(SK)+".png", dpi = 300)
            plt.clf()




































############################## Old versions to keep ################################




# Dev version obsolete now
class KP_dev(object):
    """
        This class is a development version of the knickpoint algorithm. 
        Its aim is to deal with knickpoint picking stuffs via python code, 
        in a development aim as I am significantly changing methods relatively often.
        B.G.
    """


    def __init__(self, fpath,fprefix, binning = "source_key", basins = [], sources = [], bd_KDE = 5, first_min_coeff = 0.001, load_last = False, save_last = True, remove_nodata = True,
        remove_river_threshold_node = 0, min_chi_dist_for_kp = 0.1, max_coeff = 1):
        """
            Initialization method, it creates a Knickpoint object and preprocess the stat before plotting.

        """
        
        print("Let me first preprocess and check your files")

        # Loading the attributes
        self.fpath = fpath # the path of your file : /home/your/path/
        self.fprefix = fprefix # the common prefix of all your files
        self.binning = binning # The binning method you want
        self.basins = basins
        self.sources = sources
        self.bd_KDE = bd_KDE
        self.dict_of_KDE_ksn = {}
        self.dict_of_KDE_dksn = {}
        self.first_min_coeff = first_min_coeff
        self.min_chi_dist_for_kp = min_chi_dist_for_kp
        self.max_coeff = max_coeff


        # Loading the file

        self.df = pd.read_csv(self.fpath+self.fprefix+"_KsnKn.csv")
        
        self.CNMC = pd.read_csv(self.fpath+self.fprefix+"_MChiSegmented_Ks.csv")
        if(remove_nodata):
            self.CNMC = self.CNMC[self.CNMC["m_chi"]>= 0]

        # TODO remove small river option

        print("I have loaded the files from the cpp code")

        if(len(basins)>0 or len(sources)>0):
            self.sort_my_df()

        if (remove_river_threshold_node>0):
            self.remove_small_river(remove_river_threshold_node)


        # sort the basins and/or source_keys

        # Correcting the first knickpoint of each rivers
        self.links = pd.read_csv(self.fpath+self.fprefix+"_SourcesLinks.csv")

        if(load_last):
            self.df = pd.read_pickle(self.fpath+self.fprefix+"_LASTDATASET.boris")
        else:
            self.add_knickpoints_from_source()
            print("I have preprocessed your selected sources/basin")


            # Derivative per river
            print("derivative processing")
            self.derivative_per_river()
            
            self.df = self.df.replace([np.inf, -np.inf], np.nan)
            self.df = self.df.dropna()


            # Get the KDE
            print("KDE calculation")
            self.KDE_from_scipy(bd_KDE = self.bd_KDE)

            # Get the outliers
            print("outliers selection")
            self.Select_the_outliers(first_min_coeff = self.first_min_coeff, max_coeff = self.max_coeff)

            # Now merging the knickpoints from composite break in slope
            print("merging the composite knickpoints") 
            self.merge_composite_knickpoints(chi_th = self.min_chi_dist_for_kp)
            #self.final_out = self.df[self.df["out_KDE_ksn"] >0]
            print("\n#\n##\n###")
            print("I have done my statistics and I have selected the main variation in Mx/k_sn")

        print("I have ingested and preprocessed your dataset, I am ready to plot")
        if(save_last):
            self.df.to_pickle(self.fpath+self.fprefix+"_LASTDATASET.boris")
            self.final_out.to_pickle(self.fpath+self.fprefix+"_LASTDATASET_OUT.boris")


    def remove_small_river(self, threshold):
        """
            Preprocessing alluviating function that remove the rivers under a certain umber of nodes
        """
        out_df = pd.DataFrame(data = None, columns = self.CNMC.columns)
        out_df_2 = pd.DataFrame(data = None, columns = self.df.columns)
        for source in self.CNMC["source_key"].unique():
            temp_df = self.CNMC[self.CNMC["source_key"] == source]
            temp_df2 = self.df[self.df["source_key"] == source]
            if(temp_df.shape[0]>threshold):
                out_df = pd.concat([out_df,temp_df])
                out_df_2 = pd.concat([out_df_2,temp_df2])

        # saving the new state of the global variables and resetting the index for more readability if someone wants to use the index base thing
        self.CNMC = out_df.copy()
        self.CNMC = self.CNMC.reset_index(drop = True)
        self.df = out_df_2.copy()
        self.df = self.df.reset_index(drop = True)




    def add_knickpoints_from_source(self):
        """
        Function to alleviate the initialization function
        """
        
        out_df = pd.DataFrame(data = None, columns = self.df.columns)
        for source in self.df["source_key"].unique():
            # getting the right dfs
            temp_df = self.df[self.df["source_key"] == source]
            temp_CNDF = self.CNMC[self.CNMC["source_key"] ==  source]
            if(temp_df.shape[0]>0 and temp_CNDF.shape[0]>0):
                # make sure that it is in the right order
                temp_inter = self.links[self.links["source_key"] == source]
                temp_df = temp_df.sort_values("chi")

                #creating the date for the first row           
                data = []
                singularity = temp_CNDF.sort_values("elevation")
                singularity = singularity.iloc[0]
                data = pd.DataFrame.from_dict({"X":singularity["X"], "Y": singularity["Y"], "latitude" : singularity["latitude"], "longitude" : singularity["longitude"], "elevation": singularity["elevation"],"flow_distance": singularity["flow_distance"], "chi": temp_inter["chi"],"drainage_area" : singularity["drainage_area"],"ksn": (singularity["m_chi"] - temp_inter["m_chi"]),"rksn":  (singularity["m_chi"] / temp_inter["m_chi"]),"rad" : -9999, "cumul_ksn" : -9999,  "cumul_rksn" : -9999,"cumul_rad" : -9999, "source_key" : source, "basin_key" : singularity["basin_key"]})
                # implementing the dataset
                # print(temp_df["chi"].iloc[0])
                # print("###############################################")
                # print(temp_df["chi"].iloc[1])
                # print("###############################################")
                # print("###############################################")
                temp_df = pd.concat([data, temp_df], ignore_index = True)
                temp_df = temp_df.sort_values("chi")
                # print(temp_df["chi"].iloc[0])
                # print("###############################################")
                # print(temp_df["chi"].iloc[1])
                # quit()

                #preparing the save
                out_df = pd.concat([out_df,temp_df], ignore_index = True)

        # saving the new state of the global variables and resetting the index for more readability if someone wants to use the index base thing
        self.df = out_df.copy()
        self.df = self.df.reset_index(drop = True)






    def sort_my_df(self):
        """
        Another function to alleviate the main one, it sorts the df, removing the unwanted sources and basins
        """

        ### TO DO ###
        # Add more exceptions and warning 


        # I am sorting the df by basins first and then by sources
        if (len(self.basins)>0 and len(self.sources)>0):
            
            print("\n \n \n WARNING: You gave me a list a basins_keys and sources_keys to sort your basin. If your sources are not in the basin you have choosen they won't be saved as well! WARNING \n")

        if(len(self.basins)>0):
            self.df = self.df[self.df["basin_key"].isin(self.basins)]
            self.CNMC = self.CNMC[self.CNMC["basin_key"].isin(self.basins)]
        if(len(self.sources)>0):
            self.df = self.df[self.df["source_key"].isin(self.sources)]
            self.CNMC = self.CNMC[self.CNMC["source_key"].isin(self.sources)]



    def derivative_per_river(self):
        """
            d(delta(m_chi))/d(chi).

        """

        out_df = pd.DataFrame(data = None, columns = self.df.columns)
        out_df["deriv_delta_ksn"] = pd.Series(data = None, index = out_df.index)
        out_df["deriv_ksn"] = pd.Series(data = None, index = out_df.index)

        for source in self.df["source_key"].unique():
            # getting the right dfs
            temp_df = self.df[self.df["source_key"] == source]
            # derivative
            t = (temp_df["ksn"]/temp_df["chi"].diff())
            

            t.iloc[0] = 0
            temp_df["deriv_ksn"] = pd.Series(data = t.abs(), index =temp_df.index)

            t = (temp_df["ksn"].diff()/temp_df["chi"].diff())
            t.iloc[0] = 0
            temp_df["deriv_delta_ksn"] = pd.Series(data = t.abs(), index =temp_df.index)
            # The first value is NaN, resetting to 0, it won't have any effect on the result and most of the plotting routines panic when NaN values are involved
            out_df = pd.concat([out_df,temp_df])


        # saving the new state of the global variables and resetting the index for more readability if someone wants to use the index base thing
        self.df = out_df.copy()
        self.df = self.df.reset_index(drop = True)


    def KDE_from_scipy(self, bd_KDE = 5):
        """
        Function to alleviate the initialization. Calculate the KDE
        """

        from scipy.stats import gaussian_kde as gKDE

        out_df = pd.DataFrame(data = None, columns = self.df.columns)
        out_df["KDE_ksn"] = pd.Series(data = None, index = self.df.index)
        out_df["KDE_delta_ksn"] = pd.Series(data = None, index = self.df.index)


        for source in self.df[self.binning].unique():

            try:
                temp_df = self.df[self.df[self.binning] == source]
                self.kernel_deriv_ksn = gKDE(temp_df["deriv_ksn"].values)
                self.kernel_deriv_ksn.set_bandwidth(bw_method = bd_KDE)
                #self.kernel_deriv_ksn.set_bandwidth(bw_method=self.kernel_deriv_ksn.factor )
                self.dict_of_KDE_ksn[source] = self.kernel_deriv_ksn

                temp_df["KDE_ksn"] = pd.Series(data = self.kernel_deriv_ksn(temp_df["chi"].values), index = temp_df.index)


                self.kernel_deriv_delta_ksn = gKDE(temp_df["deriv_delta_ksn"].values)
                self.kernel_deriv_delta_ksn.set_bandwidth(bw_method = bd_KDE)
                #self.kernel_deriv_delta_ksn.set_bandwidth(bw_method=self.kernel_deriv_delta_ksn.factor )
                self.dict_of_KDE_dksn[source] = self.kernel_deriv_delta_ksn

                temp_df["KDE_delta_ksn"] = pd.Series(data = self.kernel_deriv_ksn( temp_df["chi"].values), index = temp_df.index)

                out_df = pd.concat([out_df,temp_df])
            except ValueError:
                    print("I am ignoring source " +str(source) +", the river is probably too small")


        # saving the new state of the global variables and resetting the index for more readability if someone wants to use the index base thing
        self.df = out_df.copy()
        self.df = self.df.reset_index(drop = True)

    def Select_the_outliers(self, first_min_coeff = 0.01, max_coeff = 1):
        """
         Another function to alluviate the initialization one. this one add columns to the dataset to indentify the outliers
         I will implement the different methods!
         Author: B.G. 02/12/2017
        """
        out_df = pd.DataFrame(data = None, columns = self.df.columns)
        out_df["out_KDE_ksn"] = pd.Series(data = None, index = self.df.index)
        out_df["out_KDE_dksn"] = pd.Series(data = None, index = self.df.index)
        ignored_sources = []
        for source in self.df[self.binning].unique():
            

            if not np.isnan(source):
                this_coeff = first_min_coeff
                failure = True
                
                this_df = self.df[self.df[self.binning] == source]
                if(this_df.shape[0] >0):
                    this_df["out_KDE_ksn"] = pd.Series(data = np.zeros(this_df.shape[0]), index = this_df.index)
                    X = np.linspace(0,this_df["deriv_ksn"].max(), 1000)
                    Y = self.dict_of_KDE_ksn[source](X)
                    this_coeff = first_min_coeff*Y.max()
                    while(failure):
                        try:
                            th = np.min(X[Y<this_coeff])
                            this_df["out_KDE_ksn"][this_df["deriv_ksn"]>=th] = 1

                            # this_df["out_KDE_dksn"] = pd.Series(data = np.zeros(this_df.shape[0]), index = this_df.index)
                            # X = np.linspace(0,this_df["deriv_ksn"].max(), 1000)
                            # Y = self.dict_of_KDE_dksn[source](X)
                            # th = np.min(X[Y<first_min_coeff*Y.max()])
                            # this_df["out_KDE_dksn"][this_df["deriv_delta_ksn"]>=th] = 1

                            out_df = pd.concat([out_df, this_df])
                            failure = False

                        except ValueError:
                            this_coeff = Y.min()*1.1
                            if th>max_coeff:
                                print("ignoring source: %s" %(source))
                                ignored_sources.append(source)
                                failure = False
                else:
                    ignored_sources.append(source)


        
        print("I ignored " + str(len(ignored_sources))+ " sources out of " +str(self.df[self.binning].unique().shape[0]) + ", usually because they are knickpointless.")

        # saving the new state of the global variables and resetting the index for more readability if someone wants to use the index base thing
        self.df = out_df.copy()
        self.df = self.df.reset_index(drop = True)  


    def merge_composite_knickpoints(self, chi_th = 0.3):
        """
            alleviating function for initialization: merge the composite knickpoints that are detected in between two fuzzy segment:
            
            |oooo
            |
            |    ooooo
        ksn |         o
            |           o
            |           o
            |             ooooo
            |___________________
                    chi

                       ^



        """

        # looping trhough the sources
        out_df = pd.DataFrame(data = None, columns = self.df.columns)
        for source in self.df["source_key"].unique():

            # Sorting some stuffs
            temp_df = self.df[self.df["source_key"] == source]
            temp_df = temp_df[temp_df["out_KDE_ksn"]>0]
            temp_df = temp_df.sort_values("chi")
            if(temp_df.shape[0] > 1):

                # Delta chi for each knickpoints. If there is a serie of really close one, I merge them
                temp_df["dchi_kp_out"] = pd.Series(data = temp_df["chi"].diff(), index = temp_df.index)
                temp_df["group_id"] = pd.Series(data = np.zeros(temp_df.shape[0]), index = temp_df.index)
                temp_df["dchi_kp_out"].iloc[0] = 0
                incrementatorus = 1
                num_line = 0
                last_tested = True
                for index,row in temp_df.iterrows():
                    if(num_line>0):
                        if row["dchi_kp_out"] <= chi_th:
                            temp_df["group_id"].iloc[num_line-1] = incrementatorus
                            temp_df["group_id"].iloc[num_line] = incrementatorus
                            last_tested = True
                        elif(last_tested == True):
                            last_tested = False
                            incrementatorus += 1
                    num_line +=1

                for un in temp_df["group_id"].unique() :
                    if un != 0:
                        dalaf = temp_df[temp_df["group_id"] == un]
                        olaf = pd.DataFrame(data = [dalaf.mean().values], columns = dalaf.columns.values)
                        olaf["ksn"].iloc[0] = dalaf["ksn"].sum()
                        out_df = pd.concat([out_df,olaf])

                
                print temp_df["group_id"].unique()
                out_df = pd.concat([out_df,temp_df[temp_df["group_id"]==0]])
            elif (temp_df.shape[0] == 1):
                out_df = pd.concat([out_df,temp_df])

        self.final_out = out_df.copy()
        self.final_out = self.final_out.reset_index(drop = True)








    ################## Plotting routines for the object #########################     



    def plot_KDE(self, method = "deriv_ksn", group = "source_key"):
        """
            Statistical plot to calibrate the KDE per rivers
        """
        print("I am going to print one plot per rivers to display the statistics I used. you can calibrate your data using these plots")
                # check if a directory exists for the chi plots. If not then make it.
        svdir = self.fpath+'statplots/'
        if not os.path.isdir(svdir):
            os.makedirs(svdir)
        if method == "deriv_ksn":
            out_method = "out_KDE_ksn"

        elif method == "deriv_delta_ksn":
            out_method = "out_KDE_dksn"


        
        for source in self.df[group].unique():
            if not np.isnan(source):
                this_df = self.df[self.df[group] == source]
                if(this_df.shape[0]>0):
                    try:
                        fig = plt.figure(1, facecolor='white',figsize=(9,5))
                        gs = plt.GridSpec(100,100,bottom=0.10,left=0.10,right=0.95,top=0.95)
                        #ax1 = fig.add_subplot(gs[0:100,0:100])
                        ax2 = fig.add_subplot(gs[0:100,0:100])

                        #print this_df["KDE_ksn"]
                        #ax1.scatter(this_df["chi"], this_df["elevation"], s = 1, c = "b", lw = 0, label = "")
                        X = np.linspace(0,this_df[method].max(), 1000)
                        Y = self.dict_of_KDE_ksn[source](X)
                        #print Y
                        ax2.fill_between(X, 0, Y)
                        ax2.scatter(this_df[method], np.zeros(this_df.shape[0])+0.10*(Y.min()+Y.max()), s = 10, marker = "+", lw = 1, c = 'k' , label = "d(k_sn)/d(chi)")
                        ou = this_df[this_df[out_method]==1]
                        ax2.scatter(ou[method], np.zeros(ou.shape[0])+0.10*(Y.min()+Y.max()), s = 50, lw = 1, facecolors = "None", edgecolors = "r" , label = "outliers")
                        ax2.scatter(X, Y, s =5, lw = 0, c = "k", label = "KDE d(ksn)/d(chi)")
                        this_df = this_df.sort_values(method)
                        
                        ax2.set_ylim(Y.min(),Y.max())
                        ax2.set_xlabel(r'$\frac{\delta M_{\chi}}{\delta\chi}}$')
                        ax2.set_ylabel("PDF from KDE")
                        #ax2.set_xlim(0,1000)
                        #ax2.scatter(this_df["deriv_delta_ksn"], this_df["KDE_delta_ksn"], s = 100, marker = "x", lw = 1, c = "k",label = "KDE d2(ksn)/d(chi)")
                        ax2.legend()
                        plt.title(group + " #" +str(source)+", nodes: " + str(this_df.shape[0]))

                        plt.savefig(svdir+"KDE_plot_source_" + str(source) + "_out_" + method+ ".png", dpi = 300)
                        plt.clf()
                    except KeyError:
                        print("Ignoring %s #%s, something wrong happened in the KDE estimation"%(group,source))

    def plot_mchi_segments(self, method = "deriv_ksn", group = "source_key"):
        """
            Statistical plot to calibrate the KDE per rivers
        """
        print("I am going to plot the m_chi/ksn segments to check my knickpoints")
                # check if a directory exists for the chi plots. If not then make it.
        svdir = self.fpath+'river_plots/'
        if not os.path.isdir(svdir):
            os.makedirs(svdir)

        if method == "deriv_ksn":
            out_method = "out_KDE_ksn"

        elif method == "deriv_delta_ksn":
            out_method = "out_KDE_dksn"

        
        for source in self.df[group].unique():
            if not np.isnan(source):
                this_df = self.df[self.df[group] == source]
                this_MCdf = self.CNMC[self.CNMC[group] == source]
                fig = plt.figure(2, facecolor='white',figsize=(9,5))
                gs = plt.GridSpec(100,100,bottom=0.10,left=0.10,right=0.95,top=0.95)
                ax1 = fig.add_subplot(gs[0:100,0:100])
                ax2 = fig.add_subplot(gs[0:100,0:100], facecolor = "None")

                ax2.scatter(this_MCdf["chi"], this_MCdf["m_chi"], s = 1, lw = 1, c = this_MCdf["m_chi"], cmap = "RdBu_r")
                ou = this_df[this_df[out_method]==1]
                ax1.scatter(ou["chi"], ou["ksn"], s = 40, c = ou["source_key"], lw = 1,marker = "+", label = "outliers before merging", cmap = "jet")
                ax1.vlines(ou["chi"], ou["ksn"].min() , ou["ksn"].max(), linestyles  = "dashdot", lw = 0.5 )
                ouf = self.final_out[self.final_out[group] == source]
                ax1.scatter(ouf["chi"], ouf["ksn"], marker = "x", s = 50, lw = 1,c = "g", label = "outliers after merging" )


                ax2.set_xlabel(r'$\chi$')
                ax2.set_ylabel(r"$M_\chi$")
                #ax2.set_xlim(0,1000)
                #ax2.scatter(this_df["deriv_delta_ksn"], this_df["KDE_delta_ksn"], s = 100, marker = "x", lw = 1, c = "k",label = "KDE d2(ksn)/d(chi)")
                ax1.xaxis.set_visible(False)
                ax1.yaxis.set_ticks_position("right")
                ax1.yaxis.set_label_position("right")
                ax1.set_ylabel(r'$\Delta K_{sn}$')
                ax1.set_xlim(ax2.get_xlim())
                ax1.legend()

                plt.title(group + " #" +str(source))

                plt.savefig(svdir+"M_chi_plot_"+ str(group)+ "_" + str(source) + "_out_" + method+ ".png", dpi = 300)
                plt.clf()


    def plot_chi_profiles(self, method = "deriv_ksn", group = "source_key"):
        """
            Statistical plot to calibrate the KDE per rivers
        """
        print("I am going to print the chi profile per river/basins depending on what you asked")
                # check if a directory exists for the chi plots. If not then make it.
        svdir = self.fpath+'river_plots/'
        if not os.path.isdir(svdir):
            os.makedirs(svdir)

        if method == "deriv_ksn":
            out_method = "out_KDE_ksn"

        elif method == "deriv_delta_ksn":
            out_method = "out_KDE_dksn"

        
        for source in self.df[group].unique():
            if not np.isnan(source):
                this_df = self.df[self.df[group] == source]
                this_MCdf = self.CNMC[self.CNMC[group] == source]
                fig = plt.figure(2, facecolor='white',figsize=(9,5))
                gs = plt.GridSpec(100,100,bottom=0.10,left=0.10,right=0.95,top=0.95)
                ax1 = fig.add_subplot(gs[0:100,0:90])
                ax2 = fig.add_subplot(gs[0:100,0:90], facecolor = "None")
                axa = fig.add_subplot(gs[0:100,90:100])
                cax = fig.add_subplot(gs[0:100,91:95])


                ax2.scatter(this_MCdf["chi"], this_MCdf["elevation"], s = 1, lw = 0, c = this_MCdf["source_key"], cmap = "RdBu_r")
                ou = self.final_out[self.final_out[group] == source]
                cb = ax1.scatter(ou["chi"], ou["elevation"], s = 50, c = ou["ksn"], lw = 0,marker = "o", label = "Knickpoint", cmap = "RdBu_r")

                ax2.set_xlabel(r'$\chi$')
                ax2.set_ylabel("elevation (m)")
                #ax2.set_xlim(0,1000)
                #ax2.scatter(this_df["deriv_delta_ksn"], this_df["KDE_delta_ksn"], s = 100, marker = "x", lw = 1, c = "k",label = "KDE d2(ksn)/d(chi)")
                ax1.xaxis.set_visible(False)
                ax1.yaxis.set_visible(False)

                ax1.set_xlim(ax2.get_xlim())
                ax1.set_ylim(ax2.get_ylim())
                ax1.set_title(group + " #" +str(source))
                ax1.legend(loc = 2)
                

                plt.colorbar(cb,cax = cax, ax = axa)
                axa.axis("off")

                plt.savefig(svdir+"Chi_profile_"+ str(group)+ "_" + str(source) + "_out_" + method+ ".png", dpi = 300)
                plt.clf()


    def map_of_knickpoint(self, color_extent_Mx = [], method = "deriv_ksn", cut_off = 0, utm_coord = False, color_extent_kp = []):
        """
            This function plot a map of the knickpoints in a latitude/longitude way

        """
        print("Now printing a map of the knickpoint and ksn")
        svdir = self.fpath+'raster_plots/'
        if not os.path.isdir(svdir):
            os.makedirs(svdir)

        # preprocessing data

        if(len(color_extent_Mx) >0):
            cmin = color_extent_Mx[0]
            cmax = color_extent_Mx[1]
        else:
            cmin = self.CNMC["m_chi"].min()
            cmax = self.CNMC["m_chi"].max()

        if(utm_coord):
            x_col = "X"
            y_col = "Y"
        else:
            x_col = "longitude"
            y_col = "latitude"


        if method == "deriv_ksn":
            method = "out_KDE_ksn"
        elif method == "deriv_delta_ksn":
            method = "out_KDE_dksn"

        fig = plt.figure(1, facecolor='black',figsize=(9,5))
        gs = plt.GridSpec(100,100,bottom=0.10,left=0.10,right=0.95,top=0.95)
        ax1 = fig.add_subplot(gs[0:100,0:100], facecolor = "black")


        scale = self.CNMC["drainage_area"] / self.CNMC["drainage_area"].max()
        ax1.scatter(self.CNMC[x_col], self.CNMC[y_col],cmap = "RdBu_r", s = scale, c = self.CNMC["m_chi"], vmin = cmin, vmax = cmax, label = r"$M_{\chi}$" )
        
        plotting_df = self.final_out
        
        if(cut_off >0):
            plotting_df = plotting_df[plotting_df["ksn"].abs()>=cut_off]
        
        if(len(color_extent_kp)>0):
            kpmin = color_extent_kp[0]
            kpmax = color_extent_kp[1]
        else:
            kpmin = plotting_df["ksn"].min()
            kpmax = plotting_df["ksn"].max()


        cb = ax1.scatter(plotting_df[x_col], plotting_df[y_col], marker = "o", s = 26, lw = 1, c = plotting_df["ksn"], label = "knickpoint", cmap = "autumn_r", vmin = kpmin, vmax = kpmax)
        ax1.legend()

        plt.colorbar(cb)

        figname = svdir+self.fprefix +"_mapofKP.png"
        plt.savefig(figname, dpi = 300)
        plt.clf()



        



































##############################################################################################################################################################
######################################## This is the class I will implement to provide nice plots callable from the command-line #############################
########################################   Deprecated, I'll rewrite it when the method will be finally confirmed and finalized   #############################
##############################################################################################################################################################

class KnickInfo(object):
    """
    This create a KnickInfo object, my plan is to incorporate and preprocess the informations for knickstuff analysis to save time when saving figures
    Author: BG - 14/11/2016
    """
    def __init__(self,fpath,fprefix, method = 'ksn',binning = 'general', outlier_detection = False, basin_list = [], weighting = True, source_key_selection = []):

        print("Preprocessing and loading your knickpoint/zone dataset")
        self.fpath = fpath
        self.fprefix = fprefix
        self.method = method
        self.binning = binning
        self.outlier = outlier_detection

        # Loading the data
        self.knickpoint_raw = Helper.ReadKnickpointCSV(self.fpath,self.fprefix)
        self.knickzone_raw = Helper.ReadKnickzoneCSV(self.fpath, self.fprefix)
        self.chanNet = Helper.ReadMChiSegCSV(self.fpath, self.fprefix,type = 'knickpoint')

        if(len(source_key_selection) >0):
            print("I am selecting the requested sources")
            self.knickpoint_raw = self.knickpoint_raw[self.knickpoint_raw["source_key"].isin(source_key_selection)]
            self.knickzone_raw = self.knickzone_raw[self.knickzone_raw["source_key"].isin(source_key_selection)]
            self.chanNet = self.chanNet[self.chanNet["source_key"].isin(source_key_selection)]


        if(len(basin_list) > 0):
            print("I am selecting your data for the following basins:")
            print(basin_list)

        # Processing the outlier
        if(outlier_detection):
            print("I am selecting outliers by %s, binned by %s"%(method,binning))
            if(weighting):
                self.knickzone_out = get_outliers_from_DF(self.knickzone_raw, method = 'Wg'+self.method, binning = self.binning)
                self.scaling_factor = 'Wg'+ method
                # self.knickzone_base_to_lip = self.get_base_to_lip(self.knickzone_out)
            else:
                self.knickzone_out = get_outliers_from_DF(self.knickzone_raw, method = self.method, binning = self.binning)
                self.scaling_factor = method
                self.knickzone_base_to_lip =self.get_base_to_lip(self.knickzone_out)
            self.knickpoint_out = get_outliers_from_DF(self.knickpoint_raw, method = self.method, binning = self.binning)
        else:
            self.scaling_factor = method
            # self.knickzone_base_to_lip = self.get_base_to_lip(self.knickzone_raw)
        # initializing knickzone info
                # generating the knickzones
        if(self.outlier):
            self.get_knickzone_segment_for_rasterplot(self.knickzone_out)
        else:
            self.get_knickzone_segment_for_rasterplot(self.knickzone_raw)

        tempDF = self.knickpoint_raw
        
        
        self.outlier_KDS = SUT.get_outlier_from_KernelDensityStuff(tempDF, column = "deriv_ksn", binning = "source_key", threshold = 6, method = "gaussian", sort_by = "chi")



    def raster_plot_knickpoint(self, size_format='ESURF', FigFormat='png'):
    
        # check if a directory exists for the chi plots. If not then make it.
        raster_directory = self.fpath+'raster_plots/'
        if not os.path.isdir(raster_directory):
            os.makedirs(raster_directory)

        # Set up fonts for plots
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['arial']
        rcParams['font.size'] = 8

        # set figure sizes based on format
        if size_format == "geomorphology":
            fig_width_inches = 6.25
        elif size_format == "big":
            fig_width_inches = 16
        else:
            fig_width_inches = 4.92126


        # get the rasters
        raster_ext = '.bil'
        BackgroundRasterName = self.fprefix+raster_ext
        HillshadeName = self.fprefix+'_hs'+raster_ext
        BasinsName = self.fprefix+'_AllBasins'+raster_ext

        
        # create the map figure
        MF = MapFigure(HillshadeName, self.fpath,coord_type="UTM_km", alpha = 0.7)
        
        # plot the basin outlines
        Basins = LSDP.GetBasinOutlines(self.fpath, BasinsName)
        MF.plot_polygon_outlines(Basins, linewidth=0.5)

        # add the channel network without color
        ChannelPoints = LSDP.LSDMap_PointData(self.chanNet, data_type = "pandas", PANDEX = True)
        MF.add_point_data(ChannelPoints,show_colourbar="False", scale_points=True,scaled_data_in_log= True, column_for_scaling='drainage_area',alpha=0.1,max_point_size = 0.5,min_point_size = 0.1,zorder=100)
        # Selecting the knickpoints

        if(self.outlier):
            KdfPoints = LSDP.LSDMap_PointData(self.knickpoint_out, data_type = "pandas", PANDEX = True)
        else:
            KdfPoints = LSDP.LSDMap_PointData(self.knickpoint_raw, data_type = "pandas", PANDEX = True)

        MF.add_point_data(KdfPoints,this_colourmap = "RdBu_r",column_for_plotting = "sign",show_colourbar="False", scale_points=True, scaled_data_in_log= False, column_for_scaling=self.method,alpha=1,max_point_size = 15,min_point_size = 1,zorder=200)

        if(self.outlier):
            ImageName = raster_directory+self.fprefix+"_Kp_"+self.method+'_'+self.binning+".png"
        else:
            ImageName = raster_directory+self.fprefix+"_Kp_"+self.method+".png"
        MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 500) # Save the figure

    

    def raster_plot_knickzone(self, size_format='ESURF', FigFormat='png', comparison_point = []):

        # check if a directory exists for the chi plots. If not then make it.
        raster_directory = self.fpath+'raster_plots/'
        if not os.path.isdir(raster_directory):
            os.makedirs(raster_directory)

        # Set up fonts for plots
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['arial']
        rcParams['font.size'] = 8

        # set figure sizes based on format
        if size_format == "geomorphology":
            fig_width_inches = 6.25
        elif size_format == "big":
            fig_width_inches = 16
        else:
            fig_width_inches = 4.92126

        # get the rasters
        raster_ext = '.bil'
        BackgroundRasterName = self.fprefix+raster_ext
        HillshadeName = self.fprefix+'_hs'+raster_ext
        BasinsName = self.fprefix+'_AllBasins'+raster_ext

        
        # create the map figure
        MF = MapFigure(HillshadeName, self.fpath,coord_type="UTM_km", alpha = 0.7)
        
        # plot the basin outlines
        Basins = LSDP.GetBasinOutlines(self.fpath, BasinsName)
        MF.plot_polygon_outlines(Basins, linewidth=0.5)

        # add the channel network without color
        ChannelPoints = LSDP.LSDMap_PointData(self.chanNet, data_type = "pandas", PANDEX = True)
        MF.add_point_data(ChannelPoints,show_colourbar="False", scale_points=True,scaled_data_in_log= True, column_for_scaling='drainage_area',alpha=0.1,max_point_size = 0.5,min_point_size = 0.1,zorder=100)

        for kz in self.knickzone_plot:
            if(kz[1].shape[0]>0):
                if(kz[1][0] == 1): # if positive knickzone
                    this_color_zone = "#EF0808"
                else:
                    this_color_zone = "#0017FF"

                thisPointData = LSDP.LSDMap_PointData(kz[0], data_type = "pandas", PANDEX = True)
                MF.plot_segment_of_knickzone(thisPointData, color = this_color_zone, lw = 1)

        # This part is specific in the case of comparison with Nelly's KZP algorithm, it uses a different method and produces nice results -> DOI: 10.1002/2017JF004250
        if(len(comparison_point)>0):
            base = LSDP.LSDMap_PointData(comparison_point[0], data_type = "pandas", PANDEX = True)
            lips = LSDP.LSDMap_PointData(comparison_point[1], data_type = "pandas", PANDEX = True)
            MF.add_point_data(base,column_for_plotting = "nature",this_colourmap = "Wistia", manual_size = 5, zorder = 200)
            MF.add_point_data(lips,column_for_plotting = "nature",this_colourmap = "Wistia_r", manual_size = 5, zorder = 200)

        if(self.outlier):
            ImageName = raster_directory+self.fprefix+"_Kz_"+self.method+'_'+self.binning+".png"
        elif(len(comparison_point)>0):
            ImageName = raster_directory+self.fprefix+"_Kz_"+self.method+'_'+self.binning+"_compare_to_nelly.png"
        else:
            ImageName = raster_directory+self.fprefix+"_Kz_"+self.method+".png"

        
        MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 500) # Save the figure


    def raster_plot_KDS(self, size_format='ESURF', FigFormat='png', comparison_point = []):

        # check if a directory exists for the chi plots. If not then make it.
        raster_directory = self.fpath+'raster_plots/'
        if not os.path.isdir(raster_directory):
            os.makedirs(raster_directory)

        # Set up fonts for plots
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['arial']
        rcParams['font.size'] = 8

        # set figure sizes based on format
        if size_format == "geomorphology":
            fig_width_inches = 6.25
        elif size_format == "big":
            fig_width_inches = 16
        else:
            fig_width_inches = 4.92126

        # get the rasters
        raster_ext = '.bil'
        BackgroundRasterName = self.fprefix+raster_ext
        HillshadeName = self.fprefix+'_hs'+raster_ext
        BasinsName = self.fprefix+'_AllBasins'+raster_ext

        
        # create the map figure
        MF = MapFigure(HillshadeName, self.fpath,coord_type="UTM_km", alpha = 0.7)
        
        # plot the basin outlines
        Basins = LSDP.GetBasinOutlines(self.fpath, BasinsName)
        MF.plot_polygon_outlines(Basins, linewidth=0.5)

        # add the channel network without color
        ChannelPoints = LSDP.LSDMap_PointData(self.chanNet, data_type = "pandas", PANDEX = True)
        MF.add_point_data(ChannelPoints,show_colourbar="False", scale_points=True,scaled_data_in_log= True, column_for_scaling='drainage_area',alpha=0.1,max_point_size = 0.5,min_point_size = 0.1,zorder=100)

        KDS_point = LSDP.LSDMap_PointData(self.outlier_KDS, data_type = "pandas", PANDEX = True)        
        MF.add_point_data(KDS_point, manual_size = 2, this_colourmap = "gray", zorder = 500)

        # This part is specific in the case of comparison with Nelly's KZP algorithm, it uses a different method and produces nice results -> DOI: 10.1002/2017JF004250
        if(len(comparison_point)>0):
            base = LSDP.LSDMap_PointData(comparison_point[0], data_type = "pandas", PANDEX = True)
            lips = LSDP.LSDMap_PointData(comparison_point[1], data_type = "pandas", PANDEX = True)
            MF.add_point_data(base,column_for_plotting = "nature",this_colourmap = "Wistia", manual_size = 5, zorder = 200)
            MF.add_point_data(lips,column_for_plotting = "nature",this_colourmap = "Wistia_r", manual_size = 5, zorder = 200)

        ImageName = raster_directory + "KALIBKDS.png"

        
        MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 500) # Save the figure
        
        


    def chi_profiles_knickzones(self,size_format='ESURF', FigFormat='png', comparison_point = []):
        """
        Plotting routines for knickzoe chi-spaced profiles

        """

        # check if a directory exists for the chi plots. If not then make it.
        raster_directory = self.fpath+'knickzone_river/'
        if not os.path.isdir(raster_directory):
            os.makedirs(raster_directory)

        # Set up fonts for plots
        label_size = 8
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['arial']
        rcParams['font.size'] = label_size

        #now plotting
        print("I am plotting one figure per river, it can take a while. If you are processing a large area, I would recommend to select main channels")
        if(self.outlier):
            kz = self.knickzone_out
        else:
            kz = self.knickzone_raw


        for hussard in kz["source_key"].unique():
            print hussard
            print kz["source_key"].unique()
            # make a figure with required dimensions
            if size_format == "geomorphology":
                fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))            
            elif size_format == "big":
                fig = plt.figure(1, facecolor='white',figsize=(16,9))            
            else:
                fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))

            # create the axis and its position
            ## axis 1: The Chi profile and the knickpoints
            gs = plt.GridSpec(100,100,bottom=0.15,left=0.15,right=0.95,top=0.95)
            ax = fig.add_subplot(gs[0:80,0:100])
            ## axis 2: The cumul axis
            gs = plt.GridSpec(100,100,bottom=0.15,left=0.15,right=0.95,top=0.95)
            ax2 = fig.add_subplot(gs[80:100,0:100])
            
            # Selecting the data for this river
            tKdf = self.knickpoint_raw[self.knickpoint_raw["source_key"] == hussard]
            tCdf = self.chanNet[self.chanNet["source_key"] == hussard]
            tKz = kz[kz["source_key"] == hussard]


            #Sorting by Chi values, not automatic since I am probably weridly using itrator to print the map in c++
            tCdf = tCdf.sort_values("chi")
            tKz = tKz.sort_values("Achi")
            tKdf = tKdf.sort_values("chi")
            #tKzdf = tKzdf.sort_values("chi")

            # Plotting the cumul ksn_variation

            ## then plotting
            ax2.plot(tKdf['chi'], tKdf[self.method],lw = 0.5, c = '#787878' )
            ax2.fill_between(tKdf["chi"],0,tKdf[self.method], color = "k", alpha = 0.9)
            ax.fill_between(tCdf["chi"],0,tCdf["elevation"], color = '#AAAAAA' , alpha = 0.8)
            for index, row in tKz.iterrows():
                # ax2.plot([row['Achi'],row['Bchi']],[row[cumul_col],row[cumul_col]], lw = 0.75, c = '#787878')
                if(row['sign'] == 1):
                    colotempolo = "#EF0808"
                else:
                    colotempolo = "#0055FF"
                tempdfforplotting = tCdf[tCdf["chi"]>=row['Achi']]
                tempdfforplotting = tempdfforplotting[tempdfforplotting["chi"]<=row['Bchi']]

                if(row['Achi'] ==row['Bchi']):
                    lw_temp = 0.5
                else:
                    lw_temp = 0
                ax.fill_between(tempdfforplotting["chi"],0,tempdfforplotting["elevation"], color = colotempolo , alpha = 0.5,lw = lw_temp)
            
                size = row[self.scaling_factor]
                size = size/tKz[self.scaling_factor].abs().max()*100 + 10
                ax.scatter(row["Achi"],row["Aelevation"], s = size, c = colotempolo, alpha = 0.5, lw =1)

            # Plotting the Chi profiles
            ax.plot(tCdf["chi"],tCdf["elevation"], lw = 1.2 , c ='#0089B9',alpha = 1,zorder = 7)
            ax.plot(tCdf["chi"],tCdf["segmented_elevation"], lw = 0.7 , c ='k',alpha = 1,zorder = 8)
            if(len(comparison_point) == 2):
                print("I am adding the comparison_point to the profiles")
                ax.scatter(comparison_point[0]["chi"][comparison_point[0]["source_key"] == hussard], comparison_point[0]["elevation"][comparison_point[0]["source_key"] == hussard], marker = "x", c = "k", s = 50,lw = 0.8, zorder = 500)
                ax.scatter(comparison_point[1]["chi"][comparison_point[1]["source_key"] == hussard], comparison_point[1]["elevation"][comparison_point[1]["source_key"] == hussard], marker = "+", c = "k", s = 50,lw = 0.8, zorder = 500)



            # Display options
            ## setting the same Chi xlimits to display on the same scale
            ax2.set_xlim(ax.get_xlim())
            
            # setting the elevation limits
            adjuster = 0.05 * (tCdf["elevation"].max()-tCdf["elevation"].min()) # This adjuster create larger ylimits for the elevation axis.  5 % from the min/max.
            ax.set_ylim([tCdf["elevation"].min()-adjuster,tCdf["elevation"].max()+adjuster])
            ## distance from the axis
            ax.yaxis.labelpad = 7
            ax2.yaxis.labelpad = 7
            ## Color of the label, ticks and axis
            ### Cumul axis
            # ax.yaxis.label.set_color('#787878')
            # ax.tick_params(axis='y', colors='#787878')
            # ax.spines['right'].set_color('#787878')
            # ax.spines['bottom'].set_visible(False)

            ## Name of the xlabels
            ax2.set_xlabel(r'$\chi$')
            ## Name of y labels

            ax2.set_ylabel(self.method, rotation = 90, fontsize = 7)
            ax.set_ylabel('Elevation (m)', rotation = 90, fontsize = 7)
            ## Disabling the xaxis of the cumul axis as this is the same than the firts one.
            ax.xaxis.set_ticks_position('none')
            ax.set_xticklabels([])
            ax2.xaxis.set_ticks_position('bottom')
            for axil in [ax,ax2]:
                for tick in axil.yaxis.get_major_ticks():
                    tick.label.set_fontsize(7)
            # Finally setting grids to test if this looks good
            ax.grid(color = 'k', linestyle = '-', linewidth = 0.5, alpha = 0.1)
            ax2.grid(color = 'k', linestyle = '-', linewidth = 0.5, alpha = 0.1)

            # Saving the figure
            ## Building the name, it has to be specific to avoid replacing files
            if(self.outlier):
                save_name = raster_directory + self.fprefix + "_Source" + str(hussard) + self.scaling_factor + "."+FigFormat
            else:
                save_name = raster_directory + self.fprefix + "_Source" + str(hussard) + "."+FigFormat
            plt.savefig(save_name, dpi = 400)
            print("test")

            # Clearing the figure to get ready for the new one
            plt.clf()

        # Printing done to tell people that this is done
        print("done")



    def get_knickzone_segment_for_rasterplot(self, knickzone_df):

   
        print("generating knickzone raster plotting informations")
        self.knickzone_plot = []
        for source in knickzone_df["source_key"].unique():
            temp_df = knickzone_df[knickzone_df['source_key'] == source]
            
            for index,row in temp_df.iterrows():
                this_kz_chi = self.chanNet[(self.chanNet["source_key"] == source) & (self.chanNet["chi"]<= row["Bchi"]) & (self.chanNet["chi"]>= row["Achi"])]
                this_kz_sign = np.ones(this_kz_chi.shape[0])*row["sign"]

                self.knickzone_plot.append([this_kz_chi,this_kz_sign])

        print("done")

    def get_base_to_lip_from_knickpoint(self, df):
        """
        Function that returns a list of kbase_to_lip knickzones, basic version that only deal with knickpoints
        """
        self.base_to_lip_from_knickpoint =  []
        id_btlkz_basic = 0
        for source in df["source_key"].unique():
            working_df = df[df["source_key"] == source]
            working_df = working_df.sort_values("chi")
            base = 0
            lips = 0
            for index,row in working_df.iterrows():
                if base == 0:
                    if(row["ksn"] > 0):
                        base = row["chi"]
                        base_df = row
                        
                elif(row["ksn"]<0 and base >0):
                    lips = row["chi"]
                    lips_df = row
                    self.base_to_lip_from_knickpoint.append(self.chanNet[(self.chanNet["source_key"] == source) & (self.chanNet["chi"]<=lips) &  (self.chanNet["chi"]>= base)])
                    base = 0
                    lips = 0

    def plot_base_to_lip_from_knickpoint_map(self, df = 0, comparison_point = [], size_format='ESURF', FigFormat='png'):

        # check if a directory exists for the chi plots. If not then make it.

        if(isinstance(df,int)):
            df = self.knickpoint_raw
        raster_directory = self.fpath+'raster_plots/'
        if not os.path.isdir(raster_directory):
            os.makedirs(raster_directory)

        self.get_base_to_lip_from_knickpoint(df)
        # Set up fonts for plots
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['arial']
        rcParams['font.size'] = 8

        # set figure sizes based on format
        if size_format == "geomorphology":
            fig_width_inches = 6.25
        elif size_format == "big":
            fig_width_inches = 16
        else:
            fig_width_inches = 4.92126

        # get the rasters
        raster_ext = '.bil'
        BackgroundRasterName = self.fprefix+raster_ext
        HillshadeName = self.fprefix+'_hs'+raster_ext
        BasinsName = self.fprefix+'_AllBasins'+raster_ext

        
        # create the map figure
        MF = MapFigure(HillshadeName, self.fpath,coord_type="UTM_km", alpha = 0.7)
        
        # plot the basin outlines
        Basins = LSDP.GetBasinOutlines(self.fpath, BasinsName)
        MF.plot_polygon_outlines(Basins, linewidth=0.5)

        # add the channel network without color
        ChannelPoints = LSDP.LSDMap_PointData(self.chanNet, data_type = "pandas", PANDEX = True)
        MF.add_point_data(ChannelPoints,show_colourbar="False", scale_points=True,scaled_data_in_log= True, column_for_scaling='drainage_area',alpha=0.5,max_point_size = 0.5,min_point_size = 0.1,zorder=100)
        for zones in self.base_to_lip_from_knickpoint:
            if(zones.shape[0]>0):
                base = LSDP.LSDMap_PointData(zones.iloc[[0]], data_type = "pandas", PANDEX = True)
                lips = LSDP.LSDMap_PointData(zones.iloc[[-1]], data_type = "pandas", PANDEX = True)
                MF.add_point_data(base,column_for_plotting = "basin_key",this_colourmap = "RdBu", manual_size = 7, zorder = 200)
                MF.add_point_data(lips,column_for_plotting = "basin_key",this_colourmap = "RdBu_r", manual_size = 7, zorder = 200)

        # This part is specific in the case of comparison with Nelly's KZP algorithm, it uses a different method and produces nice results -> DOI: 10.1002/2017JF004250
        if(len(comparison_point)>0):
            base = LSDP.LSDMap_PointData(comparison_point[0], data_type = "pandas", PANDEX = True)
            lips = LSDP.LSDMap_PointData(comparison_point[1], data_type = "pandas", PANDEX = True)
            MF.add_point_data(base,column_for_plotting = "nature",this_colourmap = "Wistia", manual_size = 5, zorder = 200)
            MF.add_point_data(lips,column_for_plotting = "nature",this_colourmap = "Wistia_r", manual_size = 5, zorder = 200)

        # if(self.outlier):
        #     ImageName = raster_directory+self.fprefix+"_Kz_"+self.method+'_'+self.binning+".png"
        # elif(len(comparison_point)>0):
        #     ImageName = raster_directory+self.fprefix+"_Kz_"+self.method+'_'+self.binning+"_compare_to_nelly.png"
        # else:
        #     ImageName = raster_directory+self.fprefix+"_Kz_"+self.method+".png"

        
        MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = raster_directory + "base_to_lip_from_knickpoint.png", FigFormat=FigFormat, Fig_dpi = 500) # Save the figure




    def plot_base_to_lip_from_knickpoint_profile(self, df = 0, comparison_point = [], size_format='ESURF', FigFormat='png'):

        if(isinstance(df,int)):
            df = self.knickpoint_raw
        raster_directory = self.fpath+'raster_plots/'
        if not os.path.isdir(raster_directory):
            os.makedirs(raster_directory)

        self.get_base_to_lip_from_knickpoint(df)

        for source in df["source_key"].unique():
            print(source)
            # make a figure with required dimensions
            if size_format == "geomorphology":
                fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))            
            elif size_format == "big":
                fig = plt.figure(1, facecolor='white',figsize=(16,9))            
            else:
                fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))

            # create the axis and its position
            ## axis 1: The Chi profile and the knickpoints
            gs = plt.GridSpec(100,100,bottom=0.15,left=0.15,right=0.95,top=0.95)
            ax = fig.add_subplot(gs[0:100,0:100])
            CN = self.chanNet[self.chanNet["source_key"] == source]
            ax.plot(CN["chi"], CN["elevation"],zorder = 100)
            for i in self.base_to_lip_from_knickpoint:
                if i.iloc[0]["source_key"] == source:
                    ax.scatter(i.iloc[0]["chi"],i.iloc[0]["segmented_elevation"], s = 10, c = "r", zorder = 150)
                    ax.scatter(i.iloc[-1]["chi"],i.iloc[-1]["segmented_elevation"], s = 10, c = "b", zorder = 150)

            if(len(comparison_point) == 2):
                print("I am adding the comparison_point to the profiles")
                ax.scatter(comparison_point[0]["chi"][comparison_point[0]["source_key"] == source], comparison_point[0]["elevation"][comparison_point[0]["source_key"] == source], marker = "x", c = "r", s = 20,lw = 0.8, zorder = 500)
                ax.scatter(comparison_point[1]["chi"][comparison_point[1]["source_key"] == source], comparison_point[1]["elevation"][comparison_point[1]["source_key"] == source], marker = "x", c = "y", s = 20,lw = 0.8, zorder = 500)


            plt.savefig(raster_directory +"profile_base_to_lip"+str(source)+".png",dpi = 300)
            plt.clf()

    def clean_enclaved_knickzones(self,kdf):
        """
        Remove all the included knickzone
        """
        # going through the sources
        out_df = pd.DataFrame(data=None, columns=kdf.columns)
        for source in kdf["knickzone_key"].unique():
            # Selecting the source in the dataset
            working_df = kdf[kdf["knickzone_key"] == source]
            working_df = working_df.sort_values("Achi")
            # now testing hte knickzone: I am identifying the size of the potential biggest knickzone
            min_chi_kz = working_df["Achi"].min()
            max_chi_kz = working_df["Bchi"].max()
            # if a knickzone correspond to that size, that's the one and I am selecting it, and only it
            if(working_df[(working_df["Achi"] == min_chi_kz) & (working_df["Bchi"] == max_chi_kz)].shape[0] == 1):
                out_df.append(working_df[(working_df["Achi"] == min_chi_kz) & (working_df["Bchi"] == max_chi_kz)])
            else:
                #This happens if there are several knickzones derived from a large first order knickzone

                while(working_df.shape[0]>0):
                    # selection of the largest knickzone
                    if(working_df["Achi"][working_df["length"] == working_df["length"].max()].shape[0]>1):
                        print("you have more than one maximum length knickzone with exactly the same size. Deal with it in the code, I am breacking")
                        quit()
                    min_chi_kz = working_df["Achi"][working_df["length"] == working_df["length"].max()]
                    max_chi_kz = working_df["Bchi"][working_df["length"] == working_df["length"].max()]
                    # removing the enclaved ones
                    working_df.drop(working_df[((working_df["Achi"] >= min_chi_kz) & (working_df["Bchi"] != max_chi_kz) & (working_df["Achi"] <= max_chi_kz) ) & 
                                    ((working_df["Bchi"] >= min_chi_kz) & (working_df["Achi"] != min_chi_kz) & (working_df["Bchi"] >= min_chi_kz) ) ])
                    # Saving the knickzone
                    out_df.append(working_df[working_df["length"] == working_df["length"].max()])
                    # removing the knickzone
                    working_df.drop(working_df[working_df["length"] == working_df["length"].max()])
                    # the loop will continue while the working df still contain knickzones

        return out_df




############################## Data (pre-)processing treatment ##########################################

def get_outliers_from_DF(df, method = "",binning = ""):
    """
    proxy function to detect the outliers according to a specified method
    
    param:
        df (pandas Dataframe): the dataframe originally read from the _KsnKs.csv, potentially preprocessed
        method (str): method of outlier detection. Basin, source, general TOCOMPLETE AS WE ADD METHOD

    return:
        Dataframe containing the outliers
    
    Author: BG - 05/10/2017
    """
    if(binning == "general"):
        df["general"] = pd.Series(np.ones(df.shape[0]),index = df.index)
    
      
    df = SUT.extract_outliers_by_header(df,data_column_name = method, header_for_group = binning, threshold = 10)

    return df










###########################################################################################################

#################### Plotting function ####################################################################



def map_knickpoint_standard(DataDirectory, fname_prefix, size_format='ESURF', FigFormat='png', basin_list = [], mancut = 0, outlier_detection_method = "None"):
    
    """
    This creates a basic knickpoint map

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_list (list): List of the basin ID you want.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        mancut (float): manual cutoff for the data selection.
        outlier_detection_method (str): determine the outlier detection method to select the right knickpoints

    Returns:
        Shaded relief plot with the basins outlines and the knickpoints sized by intensity

    Author: BG - 05/10/2017
    """


    # check if a directory exists for the chi plots. If not then make it.
    raster_directory = DataDirectory+'raster_plots/'
    if not os.path.isdir(raster_directory):
        os.makedirs(raster_directory)

    # Set up fonts for plots
    basls = basin_list
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # set figure sizes based on format
    if size_format == "geomorphology":
        fig_width_inches = 6.25
    elif size_format == "big":
        fig_width_inches = 16
    else:
        fig_width_inches = 4.92126


    # get the rasters
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext
    BasinsName = fname_prefix+'_AllBasins'+raster_ext

    
    # create the map figure
    MF = MapFigure(HillshadeName, DataDirectory,coord_type="UTM_km", alpha = 0.7)
    
    # plot the basin outlines
    Basins = LSDP.GetBasinOutlines(DataDirectory, BasinsName)
    MF.plot_polygon_outlines(Basins, linewidth=0.5)

    # add the channel network without color
    ChannelDF = Helper.ReadMChiSegCSV(DataDirectory,fname_prefix, type = "knickpoint")
    ChannelPoints = LSDP.LSDMap_PointData(ChannelDF, data_type = "pandas", PANDEX = True)
    MF.add_point_data(ChannelPoints,show_colourbar="False", scale_points=True,scaled_data_in_log= True, column_for_scaling='drainage_area',alpha=0.1,max_point_size = 0.5,min_point_size = 0.1,zorder=100)

    # add the knickpoints plots
    
    Kdf = Helper.ReadKnickpointCSV(DataDirectory,fname_prefix)
    # Selecting the basins in case you choose specific ones
    if(len(basls)>0):
        Kdf = Kdf[kdf["basin_key"].isin(basls)]

    # Sorting data in case of manual cut_off
    if(mancut>0):
        print("I am manually cutting off the data. If you need a automated outliers detection, switch mancut option off")
        Kdf = Kdf[Kdf["rad_diff"]> mancut]
    elif(outlier_detection_method != "None"):
        print("I will now select the outliers following the method " + outlier_detection_method)
        Kdf = get_outliers_from_DF(Kdf, method = outlier_detection_method)



    KdfPoints = LSDP.LSDMap_PointData(Kdf, data_type = "pandas", PANDEX = True)
    MF.add_point_data(KdfPoints,this_colourmap = "RdBu_r",column_for_plotting = "sign",show_colourbar="False", scale_points=True, scaled_data_in_log= False, column_for_scaling='rad_diff',alpha=1,max_point_size = 15,min_point_size = 1,zorder=200)

    #Saving and stuffs
    if(outlier_detection_method == "None"):
        outlier_detection_method = "raw"  
    ImageName = raster_directory+fname_prefix+"_knickpoint_"+ outlier_detection_method +"_map."+FigFormat
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 300) # Save the figure






def basic_hist(DataDirectory, fname_prefix ,basin_list = [] , size_format="ESURF", FigFormat=".png"):
    """
    Plot a simple histogram of the knickpoint repartition
    """
    print(" \n ########################## \n I am now going to print a basic histogram of your knickpoint in the area requested \n  ")

    from matplotlib.ticker import MaxNLocator

    # check if a directory exists for the sumarry plots. If not then make it.
    raster_directory = DataDirectory+'summary_plots/'
    if not os.path.isdir(raster_directory):
        os.makedirs(raster_directory)

    # loading the file:
    df = Helper.ReadKnickpointCSV(DataDirectory, fname_prefix)

    # Selecting the basin
    if(len(basin_list)>0):
        print("Selecting the basins %s" %(basin_list))
        df = df[df["basin_key"].isin(basin_list)]
    else:
        print("I am plotting for all the basins")

    # creating the figure

    plt.clf()
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))


    # Creating the axis
    
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.15,right=0.95,top=0.95)
    ax = fig.add_subplot(gs[0:100,0:100])

    # Plotting
    print("plotting ...")
    
    ls_baboty = []
    for i in df["basin_key"].unique():
        ls_baboty.append(df["rad_diff"][df["basin_key"] == i])

    n,bins, patch = ax.hist(ls_baboty,bins = 50, stacked  = True)
    n = np.array(n)

    # setting the yticks to be int

    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    

    #setting the labels
    ax.set_xlabel("Delta ksn")
    ax.set_ylabel("count")


    #Saving and stuffs   
    ImageName = raster_directory+fname_prefix+"_hist."+FigFormat
    plt.savefig(ImageName, dpi = 300) # Save the figure
    print(" done with saving your figure " + ImageName)
    plt.clf()





def chi_profile_knickpoint(DataDirectory, fname_prefix, size_format='ESURF', FigFormat='png', basin_list = [], mancut = 0, outlier_detection_method = "None", grouping = "basin_key", segments = True):
    
    """
    This creates a chi profiles with the knickpoint on top of the profile

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_list (list): List of the basin ID you want.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        mancut (float): manual cutoff for the data selection.
        outlier_detection_method (str): determine the outlier detection method to select the right knickpoints
        segments (bool): if segments is True, it plots the Mchi segmented elevation

    Returns:
        Shaded relief plot with the basins outlines and the knickpoints sized by intensity

    Author: BG - 05/10/2017
    """

    
    # check if a directory exists for the chi plots. If not then make it.
    raster_directory = DataDirectory+'chi_profile_knickpoint/'
    if(grouping == "source_key"):
        raster_directory = DataDirectory+'chi_profile_knickpoint_by_source/'
    if not os.path.isdir(raster_directory):
        os.makedirs(raster_directory)


    # Set up fonts for plots
    basls = basin_list
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    

    # add the channel network without color
    ChannelDF = Helper.ReadMChiSegCSV(DataDirectory,fname_prefix, type = "knickpoint")
    # add the knickpoints plots
    
    Kdf = Helper.ReadKnickpointCSV(DataDirectory,fname_prefix)

    segmented_elev = Helper.ReadMChiSegCSV(DataDirectory,fname_prefix, type = "knickpoint")
    # Selecting the basins in case you choose specific ones
    if(len(basls)>0):
        Kdf = Kdf[Kdf["basin_key"].isin(basls)]
        segmented_elev = segmented_elev[segmented_elev["basin_key"].isin(basls)]

    # Sorting data in case of manual cut_off
    if(mancut>0):
        print("I am manually cutting off the data. If you need a automated outliers detection, switch mancut option off")
        Kdf = Kdf[Kdf["rad_diff"]> mancut]
    elif(outlier_detection_method != "None"):
        print("I will now select the outliers following the method " + outlier_detection_method)
        Kdf = get_outliers_from_DF(Kdf, method = outlier_detection_method, binning = "general")


    #now plotting
    for hussard in Kdf[grouping].unique():
        print(hussard)
        # make a figure
        if size_format == "geomorphology":
            fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
            
        elif size_format == "big":
            fig = plt.figure(1, facecolor='white',figsize=(16,9))
            
        else:
            fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))

        # create the axis
        gs = plt.GridSpec(100,100,bottom=0.15,left=0.15,right=0.95,top=0.95)
        ax = fig.add_subplot(gs[0:100,0:100])

        # Selecting the data for this basin
        tcdf = ChannelDF[ChannelDF[grouping] == hussard] 
        tKdf = Kdf[Kdf[grouping] == hussard]
        tsegmented_elev = segmented_elev[segmented_elev[grouping] == hussard]
        # Setting the alpha, if segmented elevation is plotted, we want to see it under the chi plots
        if (segments):
            alo = 0.7
        else:
            alo = 1

        

        # Plotting the segmented elevation -  it can be computationally expensive depending on your number of segments
        if(segments):
            tcolseg = "#01DF01"
            for seg in tsegmented_elev["segment_number"].unique():
                
                # Setting the color of the segment and inverting it each turn
                if(tcolseg == "#01DF01"):
                    tcolseg = "#2E64FE"
                else:
                    tcolseg = "#01DF01"
                # selecting the unique segment I want
                teploseg = tsegmented_elev[tsegmented_elev["segment_number"] == seg]
                #plotting the segments
                ax.plot(teploseg["chi"],teploseg["segmented_elevation"], color = tcolseg, lw = 0.68)
        ## end with the segments

        # Plotting the chi profile
        ax.scatter(tcdf["chi"],tcdf["elevation"], c = tcdf["source_key"], cmap = "Accent", s = 1, lw = 0, alpha = alo)
        # Plotting the knickpoints
        sizel = abs(tKdf["rad_diff"]) * 200
        ax.scatter(tKdf["chi"],tKdf["elevation"], c = tKdf["sign"], s = sizel, alpha = 0.6, lw = 0.5,  cmap = "gnuplot", edgecolor = "k")

        # Details
        ax.set_xlabel("Chi")
        ax.set_ylabel("z")
        # saving details
        save_name = raster_directory + fname_prefix + "_" + grouping + str(hussard) + "_"+outlier_detection_method+"."+FigFormat

        plt.savefig(save_name, dpi = 600)
        plt.clf()


    print("done")

def chi_profile_knickzone_old(DataDirectory, fname_prefix, size_format='ESURF', FigFormat='png', basin_list = [], knickpoint_value = 'delta_ksn', river_length_threshold = 0, outlier_detection_method = '', outlier_detection_binning = ''):

    """
    This creates a chi profiles with the knickpoint on top of the profile, and the knickzones information in the back.

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_list (list): List of the basin ID you want.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        knickpoint_value (str): select which knickpoint/zone to display: 'delta_ksn' for the slope of the profile; 'ratio_ksn' for the slope of the ratio variations profile; 'natural' for the angle
        river_length_threshold (int or float): ignore the rivers below a certain length. TO CODE.
       
    Returns:
        Nothing, but save one figure for each rivers.

    Author: BG - 08/11/2017
    """

    
    # check if a directory exists for the chi plots. If not then make it.
    raster_directory = DataDirectory+'knickzone_river/'
    if not os.path.isdir(raster_directory):
        os.makedirs(raster_directory)

    # Setting the knickpoint/zone value to use 
    if(knickpoint_value == 'delta_ksn'):
        knickpoint_col = "diff"
        cumul_col = "ksn"
        suffix_method = "_dksn"
        ylabel_KZ = r'$ \Delta k_{sn}$'
        ylabel_der = r'$ \frac{d\sum \Delta k_{sn}}{d\chi} $'
    elif (knickpoint_value == 'ratio_ksn'):
        knickpoint_col = "ratio"
        cumul_col = "rksn"
        suffix_method = "_rksn"
        ylabel_KZ = r'$ \Delta ratio k_{sn}$'
        ylabel_der = r'$ \frac{d\sum \Delta ratio k_{sn}}{d\chi} $'
    elif (knickpoint_value == 'natural'):
        knickpoint_col = "rad_diff"
        cumul_col = "rad"
        suffix_method = "_angle"
        ylabel_KZ = r'$ \Delta \theta$'
        ylabel_der = r'$ \frac{d\sum \Delta \theta}{d\chi} $'

    else:
        print("Unvalid value for the knickpoint method, ")


    # Set up fonts for plots
    basls = basin_list
    label_size = 8
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # read the knickpoint informations
    Kdf = Helper.ReadKnickpointCSV(DataDirectory,fname_prefix)
    Cdf = Helper.ReadMChiSegCSV(DataDirectory, fname_prefix, type = "knickpoint")
    Kzdf = Helper.ReadKnickzoneCSV(DataDirectory,fname_prefix)

    # Selecting the basins in case you choose specific ones
    if(len(basls)>0):
        Kdf = Kdf[Kdf["basin_key"].isin(basls)]
        Cdf = Cdf[Cdf["basin_key"].isin(basls)]
        Kzdf = Kzdf[Kzdf["basin_key"].isin(basls)]

    if(outlier_detection_binning != ''):
        Kzdf = get_outliers_from_DF(Kzdf, method = outlier_detection_method, binning = outlier_detection_binning)

    #now plotting
    print("I am plotting one figure per river, it can take a while. If you are processing a large area, I would recommend to select main channels")
    #for hussard in Kdf["source_key"].unique():
    for hussard in [0,19]: #TEMPORARY TEST FOR COLUMBIA CA
        # make a figure with required dimensions
        if size_format == "geomorphology":
            fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))            
        elif size_format == "big":
            fig = plt.figure(1, facecolor='white',figsize=(16,9))            
        else:
            fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))

        # create the axis and its position
        ## axis 1: The Chi profile and the knickpoints
        gs = plt.GridSpec(100,100,bottom=0.15,left=0.15,right=0.95,top=0.95)
        ax = fig.add_subplot(gs[0:80,0:100])
        ## axis 2: The cumul axis
        gs = plt.GridSpec(100,100,bottom=0.15,left=0.15,right=0.95,top=0.95)
        ax2 = fig.add_subplot(gs[80:100,0:100])
        
        # Selecting the data for this river
        tKdf = Kdf[Kdf["source_key"] == hussard]
        tCdf = Cdf[Cdf["source_key"] == hussard]
        tKzdf = Kzdf[Kzdf["source_key"] == hussard]

        #Sorting by Chi values, not automatic since I am probably weridly using itrator to print the map in c++
        tKdf = tKdf.sort_values("chi")
        tCdf = tCdf.sort_values("chi")
        #tKzdf = tKzdf.sort_values("chi")

        # Plotting the cumul ksn_variation
        ## shifting first and initial value at 0 for the variations
        tKdf.iloc[0, tKdf.columns.get_loc('chi')] = tCdf["chi"].min()
        ## then plotting
        ax2.plot(tKdf['chi'], tKdf[knickpoint_col],lw = 0.5, c = '#787878' )
        ax2.fill_between(tKdf["chi"],0,tKdf[knickpoint_col], color = "k", alpha = 0.9)
        ax.fill_between(tCdf["chi"],0,tCdf["elevation"], color = '#AAAAAA' , alpha = 0.8)
        for index, row in tKzdf.iterrows():
            # ax2.plot([row['Achi'],row['Bchi']],[row[cumul_col],row[cumul_col]], lw = 0.75, c = '#787878')
            if(row['sign'] == 1):
                colotempolo = "#EF0808"
            else:
                colotempolo = "#0055FF"
            tempdfforplotting = tCdf[tCdf["chi"]>=row['Achi']]
            tempdfforplotting = tempdfforplotting[tempdfforplotting["chi"]<=row['Bchi']]

            if(row['Achi'] ==row['Bchi']):
                lw_temp = 0.5
            else:
                lw_temp = 0
            ax.fill_between(tempdfforplotting["chi"],0,tempdfforplotting["elevation"], color = colotempolo , alpha = 0.5,lw = lw_temp)

        # Plotting the Chi profiles
        ax.plot(tCdf["chi"],tCdf["elevation"], lw = 1.2 , c ='#0089B9',alpha = 1,zorder = 7)

        # sizing the points, casting between 10 and 100
        size = tKdf[knickpoint_col].abs()
        size = size/size.max()*100 + 10

        ax.scatter(tKdf["chi"],tKdf["elevation"], c = tKdf["sign"],cmap = 'RdBu_r', s = size, alpha = 0.75, lw = 0.5, edgecolor = "k", zorder = 10)

        # Display options
        ## setting the same Chi xlimits to display on the same scale
        ax2.set_xlim(ax.get_xlim())
        
        # setting the elevation limits
        adjuster = 0.05 * (tCdf["elevation"].max()-tCdf["elevation"].min()) # This adjuster create larger ylimits for the elevation axis.  5 % from the min/max.
        ax.set_ylim([tCdf["elevation"].min()-adjuster,tCdf["elevation"].max()+adjuster])
        ## distance from the axis
        ax.yaxis.labelpad = 7
        ax2.yaxis.labelpad = 7
        ## Color of the label, ticks and axis
        ### Cumul axis
        # ax.yaxis.label.set_color('#787878')
        # ax.tick_params(axis='y', colors='#787878')
        # ax.spines['right'].set_color('#787878')
        # ax.spines['bottom'].set_visible(False)

        ## Name of the xlabels
        ax2.set_xlabel(r'$\chi$')
        ## Name of y labels
        ax2.set_ylabel(ylabel_KZ, rotation = 90, fontsize = 7)
        ax.set_ylabel('Elevation (m)', rotation = 90, fontsize = 7)
        ## Disabling the xaxis of the cumul axis as this is the same than the firts one.
        ax.xaxis.set_ticks_position('none')
        ax.set_xticklabels([])
        ax2.xaxis.set_ticks_position('bottom')
        for axil in [ax,ax2]:
            for tick in axil.yaxis.get_major_ticks():
                tick.label.set_fontsize(7)
        # Finally setting grids to test if this looks good
        ax.grid(color = 'k', linestyle = '-', linewidth = 0.5, alpha = 0.1)
        ax2.grid(color = 'k', linestyle = '-', linewidth = 0.5, alpha = 0.1)

        # Saving the figure
        ## Building the name, it has to be specific to avoid replacing files
        save_name = raster_directory + fname_prefix + "_Source" + str(hussard) + suffix_method + '_' + outlier_detection_method +  "."+FigFormat
        plt.savefig(save_name, dpi = 400)

        # Clearing the figure to get ready for the new one
        plt.clf()

    # Printing done to tell people that this is done
    print("done")

def chi_profile_knickzone_old(DataDirectory, fname_prefix, size_format='ESURF', FigFormat='png', basin_list = [], knickpoint_value = 'delta_ksn', river_length_threshold = 0):

    """
    This creates a chi profiles with the knickpoint on top of the profile, and the knickzones information in the back.

    Args:
        DataDirectory (str): the data directory with the m/n csv files
        fname_prefix (str): The prefix for the m/n csv files
        basin_list (list): List of the basin ID you want.
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).
        FigFormat (str): The format of the figure. Usually 'png' or 'pdf'. If "show" then it calls the matplotlib show() command.
        knickpoint_value (str): select which knickpoint/zone to display: 'delta_ksn' for the slope of the profile; 'ratio_ksn' for the slope of the ratio variations profile; 'natural' for the angle
        river_length_threshold (int or float): ignore the rivers below a certain length. TO CODE.
       
    Returns:
        Nothing, but save one figure for each rivers.

    Author: BG - 08/11/2017
    """

    
    # check if a directory exists for the chi plots. If not then make it.
    raster_directory = DataDirectory+'knickzone_river/'
    if not os.path.isdir(raster_directory):
        os.makedirs(raster_directory)

    # Setting the knickpoint/zone value to use 
    if(knickpoint_value == 'delta_ksn'):
        knickpoint_col = "diff"
        cumul_col = "cumul_ksn"
        deriv_cumul = "deriv_cumul_ksn"
        suffix_method = "_dksn"
        ylabel_KZ = r'$\sum \Delta k_{sn}$'
        ylabel_der = r'$ \frac{d\sum \Delta k_{sn}}{d\chi} $'
    elif (knickpoint_value == 'ratio_ksn'):
        knickpoint_col = "ratio"
        cumul_col = "cumul_rksn"
        deriv_cumul = "deriv_cumul_rksn"
        suffix_method = "_rksn"
        ylabel_KZ = r'$\sum \Delta ratio k_{sn}$'
        ylabel_der = r'$ \frac{d\sum \Delta ratio k_{sn}}{d\chi} $'
    elif (knickpoint_value == 'natural'):
        knickpoint_col = "rad_diff"
        cumul_col = "cumul_rad"
        deriv_cumul = "deriv_cumul_rad"
        suffix_method = "_angle"
        ylabel_KZ = r'$\sum \Delta \theta$'
        ylabel_der = r'$ \frac{d\sum \Delta \theta}{d\chi} $'

    else:
        print("Unvalid value for the knickpoint method, ")


    # Set up fonts for plots
    basls = basin_list
    label_size = 8
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # read the knickpoint informations
    Kdf = Helper.ReadKnickpointCSV(DataDirectory,fname_prefix)
    Cdf = Helper.ReadMChiSegCSV(DataDirectory, fname_prefix, type = "knickpoint")

    # Selecting the basins in case you choose specific ones
    if(len(basls)>0):
        Kdf = Kdf[Kdf["basin_key"].isin(basls)]
        Cdf = Cdf[Cdf["basin_key"].isin(basls)]

    #now plotting
    print("I am plotting one figure per river, it can take a while. If you are processing a large area, I would recommend to select main channels")
    #for hussard in Kdf["source_key"].unique():
    for hussard in [0,19]: #TEMPORARY TEST FOR COLUMBIA CA
        # make a figure with required dimensions
        if size_format == "geomorphology":
            fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))            
        elif size_format == "big":
            fig = plt.figure(1, facecolor='white',figsize=(16,9))            
        else:
            fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))

        # create the axis and its position
        ## axis 1: The Chi profile and the knickpoints
        gs = plt.GridSpec(99,99,bottom=0.15,left=0.15,right=0.95,top=0.95)
        ax = fig.add_subplot(gs[0:33,0:99])
        ## axis 2: The cumul axis
        gs = plt.GridSpec(99,99,bottom=0.15,left=0.15,right=0.95,top=0.95)
        ax2 = fig.add_subplot(gs[33:66,0:99])
        ## axis 3: The derivative axis
        gs = plt.GridSpec(99,99,bottom=0.15,left=0.15,right=0.95,top=0.95)
        ax3 = fig.add_subplot(gs[66:99,0:99])
        
        # Selecting the data for this river
        tKdf = Kdf[Kdf["source_key"] == hussard]
        tCdf = Cdf[Cdf["source_key"] == hussard]

        #Sorting by Chi values, not automatic since I am probably weridly using itrator to print the map in c++
        tKdf = tKdf.sort_values("chi")
        tCdf = tCdf.sort_values("chi")

        # Plotting the cumul ksn_variation
        ## shifting first and initial value at 0 for the variations
        tKdf.iloc[0, tKdf.columns.get_loc('chi')] = tCdf["chi"].min()
        ## then plotting
        ax2.plot(tKdf["chi"],tKdf[cumul_col], lw = 0.75, c = '#787878')
        ax2.fill_between(tKdf["chi"],0,tKdf[cumul_col], color = "k", alpha = 0.3)

        # Plotting the Chi profiles
        ax.plot(tCdf["chi"],tCdf["elevation"], lw = 1.2 , c ='#0089B9',alpha = 1,zorder = 7)
        ax.scatter(tKdf["chi"],tKdf["elevation"], c = tKdf["sign"],cmap = 'RdBu_r', s = tKdf[knickpoint_col].abs(), alpha = 0.75, lw = 0.5, edgecolor = "k", zorder = 10)
        
        # Plotting the derivative
        ax3.plot(tKdf["chi"],tKdf[deriv_cumul], lw = 0.7, c = '#E70B0B',alpha = 0.7,zorder = 5)

        # Display options
        ## setting the same Chi xlimits to display on the same scale
        ax2.set_xlim(ax.get_xlim())
        ax3.set_xlim(ax.get_xlim())
        
        ## distance from the axis
        ax.yaxis.labelpad = 9
        ax2.yaxis.labelpad = 9
        ax3.yaxis.labelpad = 9
        ## Color of the label, ticks and axis
        ### Cumul axis
        # ax.yaxis.label.set_color('#787878')
        # ax.tick_params(axis='y', colors='#787878')
        # ax.spines['right'].set_color('#787878')
        # ax.spines['bottom'].set_visible(False)
        ### Deriv Axis
        # ax3.yaxis.label.set_color('#E70B0B')
        # ax3.tick_params(axis='y', colors='#E70B0B')
        # ax3.spines['right'].set_color('#E70B0B')
        # ax3.spines['top'].set_visible(True)
        ## Name of the xlabels
        ax3.set_xlabel(r'$\chi$')
        ## Name of y labels
        ax2.set_ylabel(ylabel_KZ, rotation = 0, fontsize = 7)
        ax3.set_ylabel(ylabel_der, rotation = 0, fontsize = 7)
        ax.set_ylabel('Elevation (m)', rotation = 0, fontsize = 7)
        ## Disabling the xaxis of the cumul axis as this is the same than the firts one.
        ax.xaxis.set_ticks_position('none')
        ax.set_xticklabels([])
        ax2.xaxis.set_ticks_position('none')
        ax2.set_xticklabels([])
        ax3.xaxis.set_ticks_position('bottom')
        for axil in [ax,ax2,ax3]:
            for tick in axil.yaxis.get_major_ticks():
                tick.label.set_fontsize(6)
        # ax2.patch.set_visible(False)
        # ax3.patch.set_visible(False)
        ## Transparence of the background for the chi axis
        # ax2.patch.set_alpha(0.1)
        # Finally setting grids to test if this looks good
        ax.grid(color = 'k', linestyle = '-', linewidth = 0.5, alpha = 0.1)
        ax2.grid(color = 'k', linestyle = '-', linewidth = 0.5, alpha = 0.1)
        ax3.grid(color = 'k', linestyle = '-', linewidth = 0.5, alpha = 0.1)

        # Saving the figure
        ## Building the name, it has to be specific to avoid replacing files
        save_name = raster_directory + fname_prefix + "_Source" + str(hussard) + suffix_method + "."+FigFormat
        plt.savefig(save_name, dpi = 400)

        # Clearing the figure to get ready for the new one
        plt.clf()

    # Printing done to tell people that this is done
    print("done")












#################################################################################################################










#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~ THE FOLLOWING FUNCTIONS ARE TESTS AND EARLIER WORK, NOT USED WITH THE PlotKnickpointAnalysis ~~~~~~~~~#
#~~~~~~~~~~~~ KEEP THEM AS THEY ARE THE BASE OF THE AUTOMATED FUNCTIONS AND FOR TESTING PURPOSES ~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################






# Test to keep saved, it was before separating the plots in three graphs.


    #     # Selecting the data for this river
    #     tKdf = Kdf[Kdf["source_key"] == hussard]
    #     tCdf = Cdf[Cdf["source_key"] == hussard]

    #     #Sorting by Chi values, not automatic since I am probably weridly using itrator to print the map in c++
    #     tKdf = tKdf.sort_values("chi")
    #     tCdf = tCdf.sort_values("chi")

    #     # Plotting the cumul ksn_variation
    #     ## shifting first and initial value at 0 for the variations
    #     tKdf.iloc[0, tKdf.columns.get_loc('chi')] = tCdf["chi"].min()
    #     ## then plotting
    #     ax.plot(tKdf["chi"],tKdf[cumul_col], lw = 0.75, c = '#787878')
    #     ax.fill_between(tKdf["chi"],0,tKdf[cumul_col], color = "k", alpha = 0.3)

    #     # Plotting the Chi profiles
    #     ax2.plot(tCdf["chi"],tCdf["elevation"], lw = 1.2 , c ='#0089B9',alpha = 1,zorder = 7)
    #     ax2.scatter(tKdf["chi"],tKdf["elevation"], c = tKdf["sign"],cmap = 'RdBu_r', s = tKdf[knickpoint_col].abs(), alpha = 0.75, lw = 0.5, edgecolor = "k", zorder = 10)
        
    #     # Plotting the derivative
    #     ax3.plot(tKdf["chi"],tKdf[deriv_cumul], lw = 0.7, c = '#E70B0B',alpha = 0.7,zorder = 5)

    #     # Display options
    #     ## setting the same Chi xlimits to display on the same scale
    #     ax.set_xlim(ax2.get_xlim())
    #     ax3.set_xlim(ax2.get_xlim())
    #     ## set the tick/label for the sum of delta ksn on the right of the plot
    #     ax.yaxis.set_ticks_position('right')
    #     ax.yaxis.set_label_position('right')
    #     ax3.yaxis.set_ticks_position('right')
    #     ax3.yaxis.set_label_position('right')
    #     ## distance from the axis
    #     ax.yaxis.labelpad = 12
    #     ax2.yaxis.labelpad = 7.5
    #     ax3.yaxis.labelpad = 12
    #     ## Color of the label, ticks and axis
    #     ### Cumul axis
    #     ax.yaxis.label.set_color('#787878')
    #     ax.tick_params(axis='y', colors='#787878')
    #     ax.spines['right'].set_color('#787878')
    #     ax.spines['bottom'].set_visible(False)
    #     ### Deriv Axis
    #     ax3.yaxis.label.set_color('#E70B0B')
    #     ax3.tick_params(axis='y', colors='#E70B0B')
    #     ax3.spines['right'].set_color('#E70B0B')
    #     ax3.spines['top'].set_visible(False)
    #     ## Name of the xlabels
    #     ax2.set_xlabel(r'$\chi$')
    #     ## Name of y labels
    #     ax.set_ylabel(ylabel_KZ, rotation = 0)
    #     ax3.set_ylabel(ylabel_der, rotation = 0)
    #     ax2.set_ylabel('Elevation (m)')
    #     ## Disabling the xaxis of the cumul axis as this is the same than the firts one.
    #     ax.xaxis.set_ticks_position('none')
    #     ax.set_xticklabels([])
    #     ax.patch.set_visible(False)
    #     ax3.patch.set_visible(False)
    #     ## Transparence of the background for the chi axis
    #     ax2.patch.set_alpha(0.1)
    #     # Finally setting grids to test if this looks good
    #     ax.grid(color = 'k', linestyle = '-', linewidth = 0.5, alpha = 0.1)
    #     ax3.grid(color = 'k', linestyle = '-', linewidth = 0.5, alpha = 0.1)

    #     # Saving the figure
    #     ## Building the name, it has to be specific to avoid replacing files
    #     save_name = raster_directory + fname_prefix + "_Source" + str(hussard) + suffix_method + "."+FigFormat
    #     plt.savefig(save_name, dpi = 400)

    #     # Clearing the figure to get ready for the new one
    #     plt.clf()

    # # Printing done to tell people that this is done
    # print("done")




def plot_knickpoint_elevations(PointData, DataDirectory, basin_key=0, kp_threshold=0,
                               FigFileName='Image.pdf', FigFormat='pdf', size_format='ESURF', kp_type = "rad_diff"):
    """
    Function to create a plot of knickpoint elevation vs flow distance for each
    basin. Knickpoints are colour-coded by source node, and the marker size represents
    the magnitude of the knickpoint.

    Args:
        PointData: the LSDMap_PointData object with the knickpoint information
        DataDirectory (str): the data directory for the knickpoint file
        csv_name (str): name of the csv file with the knickpoint information
        basin_key (int): key to select the basin of interest
        kp_threshold (int): threshold knickpoint magnitude, any knickpoint below this will be removed (This option may be removed soon)
        kp_type (string): switch between diff and ratio data
        FigFileName (str): The name of the figure file
        FigFormat (str): format of output figure, can be 'pdf' (default), 'png', 'return', or 'show'
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).

    Returns:
        Plot of knickpoint elevations against flow distance

    Author: FJC
    """
    #PointData = LSDMap_PD.LSDMap_PointData(kp_csv_fname)
    # thin out small knickpoints
    KPData = PointData
    #KPData.ThinData(kp_type,kp_threshold)

    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
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

    # get the data

    # Soert the dataset to the basin key
    KPData.selectValue('basin_key', value = basin_key, operator = "==")

    elevation = KPData.QueryData('elevation')

    flow_distance = KPData.QueryData('flow distance')

    magnitude = KPData.QueryData(kp_type)
    print("For the plotting, if you want to manage the scale, " +kp_type + " max is "+ str(magnitude.max()) +" and min is " + str(magnitude.min()))

    source = KPData.QueryData('source_key')




    #colour by source - this is the same as the script to colour channels over a raster,
    # (BasicChannelPlotGridPlotCategories) so that the colour scheme should match
    # make a color map of fixed colors
    NUM_COLORS = len(np.unique(source))
    this_cmap = plt.cm.Set1
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    plt.cm.ScalarMappable(norm=cNorm, cmap=this_cmap)
    channel_data = [x % NUM_COLORS for x in source]

    # normalise magnitude for plotting
    min_size = np.min(magnitude)
    max_size = np.max(magnitude)
    #normSize = [100*((x - min_size)/(max_size - min_size)) for x in magnitude]
    normSize = 100 * (magnitude - min_size)/(max_size - min_size)

    # now get the channel profiles that correspond to each knickpoint source and basin
    # PointData.ThinDataFromKey('basin_key',basin_key)
    # PointData.ThinDataSelection('source_key',maskSource)
    #
    # channel_elev = PointData.QueryData('elevation')
    # channel_elev = [float(x) for x in channel_elev]
    # channel_dist = PointData.QueryData('flow_distance')
    # channel_dist = [float(x) for x in channel_dist]

    # now plot the knickpoint elevations and flow distances
    #ax.scatter(channel_dist, channel_elev, s=0.1, c='k')
    ax.scatter(flow_distance, elevation, c = channel_data, cmap=this_cmap, s = normSize, lw=0.5,edgecolors='k',zorder=100)
    ax.set_xlabel('Flow distance (m)')
    ax.set_ylabel('Elevation (m)')

    # return or show the figure
    print("The figure format is: " + FigFormat)
    if FigFormat == 'show':
        plt.show()
    elif FigFormat == 'return':
        return fig
    else:
        save_fmt = FigFormat
        plt.savefig(DataDirectory+FigFileName,format=save_fmt,dpi=500)
        fig.clf()


def plot_diff_ratio(PointData, DataDirectory, saveName = "Basic_diff_ratio", save_fmt = ".png", size_format = "ESURF", log_data = False):
    """
    Basic plot to have a general view of the knickpoints: diff against ratio colored by elevation

    Args:
        PointData: A PointData object
        DataDirectory: Where the data is saved
        saveName: save name

    returns:
        Nothing, sorry.
    Author: BG
    """
    plt.clf()
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
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

    diff = PointData.QueryData("diff")
    ratio = PointData.QueryData("ratio")

    if(log_data):
        print("I am logging the data")
        diff = np.log10(diff)
        ratio = np.log10(ratio)

    elevation =PointData.QueryData("elevation")
    ax.scatter(diff,ratio, s=0.5, lw = 0, c = elevation)
    ax.set_xlabel("Diff")
    ax.set_ylabel("Ratio")

    plt.savefig(DataDirectory+saveName+save_fmt,dpi=500)

def plot_basic_DA(PointData, DataDirectory, saveName = "Basic_DA", save_fmt = ".png", size_format = "ESURF", log_data = False):
    """
    Basic plot to have a general view of the knickpoints: drainage area against ratio and diff colored by elevation

    Args:
        PointData: A PointData object
        DataDirectory: Where the data is saved
        saveName: save name

    returns:
        Nothing, sorry.
    Author: BG
    """
    plt.clf()
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(2, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(2, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(2, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    diff = PointData.QueryData("diff")
    ratio = PointData.QueryData("ratio")
    DA = PointData.QueryData("drainage area")

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax1 = fig.add_subplot(gs[10:50,10:95])
    ax2 = fig.add_subplot(gs[50:100,10:95])



    if(log_data):
        print("I am logging the data")
        diff = np.log10(diff)
        ratio = np.log10(ratio)

    elevation = PointData.QueryData("elevation")
    DA = np.log10(DA)
    ax1.scatter(DA,ratio, s=0.7, lw = 0, c = elevation)
    ax1.set_ylabel("Ratio")
    ax1.tick_params(axis = 'x', length = 0, width = 0, labelsize = 0)
    ax1.spines['bottom'].set_visible(False)
    ax2.scatter(DA,diff,s=0.7, lw = 0, c = elevation)
    ax2.set_ylabel("Diff")
    ax2.set_xlabel("Drainage area")
    #ax2.set_xticks([1,2,3,4,5,6,7])
    ax2.tick_params(axis = 'x', labelsize = 6)
    ax1.set_xticks([4,5,6,7,8,9,10])
    ax2.set_xticks([4,5,6,7,8,9,10])
    ax2.set_xticklabels([ur"$10^{4}$",ur"$10^{5}$",ur"$10^{6}$",ur"$10^{7}$",ur"$10^{8}$",ur"$10^{9}$",ur"$10^{10}$"])

    plt.savefig(DataDirectory+saveName+save_fmt,dpi=500)

def plot_basic_FD(PointData, DataDirectory, saveName = "Basic_FD", save_fmt = ".png", size_format = "ESURF", log_data = False):
    """
    Basic plot to have a general view of the knickpoints: flow distance against ratio and diff colored by elevation

    Args:
        PointData: A PointData object
        DataDirectory: Where the data is saved
        saveName: save name

    returns:
        Nothing, sorry.
    Author: BG
    """
    plt.clf()
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(2, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(2, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(2, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    diff = PointData.QueryData("diff")
    ratio = PointData.QueryData("ratio")
    FD = PointData.QueryData("flow distance")

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax1 = fig.add_subplot(gs[10:50,10:95])
    ax2 = fig.add_subplot(gs[50:100,10:95])



    if(log_data):
        print("I am logging the data")
        diff = np.log10(diff)
        ratio = np.log10(ratio)

    elevation = PointData.QueryData("elevation")

    ax1.scatter(FD,ratio, s=0.7, lw = 0, c = elevation)
    ax1.set_ylabel("Ratio")
    ax1.tick_params(axis = 'x', length = 0, width = 0, labelsize = 0)
    ax1.spines['bottom'].set_visible(False)
    ax2.scatter(FD,diff,s=0.7, lw = 0, c = elevation)
    ax2.set_ylabel("Diff")
    ax2.set_xlabel("Flow distance")

    #ax2.tick_params(axis = 'x', labelsize = 6)
    #ax1.set_xticks([4,5,6,7,8,9,10])
    #ax2.set_xticks([4,5,6,7,8,9,10])
    #ax2.set_xticklabels([ur"$10^{4}$",ur"$10^{5}$",ur"$10^{6}$",ur"$10^{7}$",ur"$10^{8}$",ur"$10^{9}$",ur"$10^{10}$"])

    plt.savefig(DataDirectory+saveName+save_fmt,dpi=500)

def plot_basic_Z(PointData, DataDirectory, saveName = "Basic_Z", save_fmt = ".png", size_format = "ESURF", log_data = False):
    """
    Basic plot to have a general view of the knickpoints: flow distance against ratio and diff colored by elevation

    Args:
        PointData: A PointData object
        DataDirectory: Where the data is saved
        saveName: save name

    returns:
        Nothing, sorry.
    Author: BG
    """
    plt.clf()
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(2, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(2, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(2, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    diff = PointData.QueryData("diff")
    ratio = PointData.QueryData("ratio")
    Z = PointData.QueryData("elevation")

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax1 = fig.add_subplot(gs[10:50,10:95])
    ax2 = fig.add_subplot(gs[50:100,10:95])



    if(log_data):
        print("I am logging the data")
        diff = np.log10(diff)
        ratio = np.log10(ratio)

    sign = PointData.QueryData("sign")

    ax1.scatter(Z,ratio, s=0.7, lw = 0, c = sign)
    ax1.set_ylabel("Ratio")
    ax1.tick_params(axis = 'x', length = 0, width = 0, labelsize = 0)
    ax1.spines['bottom'].set_visible(False)
    ax2.scatter(Z,diff,s=0.7, lw = 0, c = sign)
    ax2.set_ylabel("Diff")
    ax2.set_xlabel("elevation")

    #ax2.tick_params(axis = 'x', labelsize = 6)
    #ax1.set_xticks([4,5,6,7,8,9,10])
    #ax2.set_xticks([4,5,6,7,8,9,10])
    #ax2.set_xticklabels([ur"$10^{4}$",ur"$10^{5}$",ur"$10^{6}$",ur"$10^{7}$",ur"$10^{8}$",ur"$10^{9}$",ur"$10^{10}$"])

    plt.savefig(DataDirectory+saveName+save_fmt,dpi=500)

def plot_outliers_x_vs_diff_ratio(PointData, DataDirectory,x_col = "elevation", saveName = "Outliers", save_fmt = ".png", size_format = "ESURF", log_data = False, ylim_ratio = [], ylim_diff = []):
    """
    Basic plot to have a general view of the knickpoints: flow distance against ratio and diff colored by elevation

    Args:
        PointData: A PointData object
        DataDirectory: Where the data is saved
        saveName: save name

    returns:
        Nothing, sorry.
    Author: BG
    """
    # Merging the dictionnary
    if(isinstance(PointData, dict)):
        print("Your data is a dictionnary of dataframes, let me create a PointTool object that contains all of these.")
        PointData = pd.concat(PointData)
        PointData = LSDMap_PD.LSDMap_PointData(PointData,data_type ="pandas", PANDEX = True)


    plt.clf()
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(2, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(2, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(2, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax1 = fig.add_subplot(gs[5:45,10:95])
    ax2 = fig.add_subplot(gs[55:100,10:95])


    diff = PointData.QueryData("diff")
    ratio = PointData.QueryData("ratio")
    Z = PointData.QueryData(x_col)
    PointData.selectValue('diff_outlier',value = True, operator = "==")
    PointData.selectValue('ratio_outlier',value = True, operator = "==")

    diffout = PointData.QueryData("diff")
    ratioout = PointData.QueryData("ratio")
    Z_out = PointData.QueryData(x_col)
    if(log_data):
        print("I am logging the data")

        diff = np.log10(diff)
        ratio = np.log10(ratio)
        if(isinstance(diffout, list) == False):
            diffout = [diffout]
            ratioout = [ ratioout]

        for i in range(len(diffout)):

            diffout[i]= np.log10(diffout[i])

            ratioout[i]= np.log10(ratioout[i])

    sign = PointData.QueryData("sign")

    ax1.scatter(Z,ratio, s=1.5, lw = 0, c = "gray")
    ax1.scatter(Z_out,ratioout, s=1.5, lw = 0, c = sign)
    ax1.set_ylabel("Ratio")
    ax1.tick_params(axis = 'x', length = 0, width = 0, labelsize = 0)
    ax1.spines['bottom'].set_visible(False)
    ax2.scatter(Z,diff,s=1.5, lw = 0, c = "gray")
    ax2.scatter(Z_out,diffout,s=1.5, lw = 0, c = sign)
    ax2.set_ylabel("Diff")
    ax2.set_xlabel(x_col)
    if(ylim_diff != []):
        ax2.set_ylim(ylim_diff[0],ylim_diff[1])
    if ylim_ratio != []:
        ax.set_ylim(ylim_ratio[0],ylim_ratio[1])

    #ax2.tick_params(axis = 'x', labelsize = 6)
    #ax1.set_xticks([4,5,6,7,8,9,10])
    #ax2.set_xticks([4,5,6,7,8,9,10])
    #ax2.set_xticklabels([ur"$10^{4}$",ur"$10^{5}$",ur"$10^{6}$",ur"$10^{7}$",ur"$10^{8}$",ur"$10^{9}$",ur"$10^{10}$"])

    plt.savefig(DataDirectory+saveName+save_fmt,dpi=500)

def plot_outliers_vs_others(PointData, DataDirectory, saveName = "Basic_diff_ratio", save_fmt = ".png", size_format = "ESURF", log_data = False):
    """
    Basic plot to have a general view of the knickpoints: diff against ratio colored by elevation

    Args:
        PointData: A PointData object or a dictionnary of dataframe containing outliers and none outliers
        DataDirectory: Where the data is saved
        saveName: save name

    returns:
        Nothing, sorry.
    Author: BG
    """
    # Merging the dictionnary
    if(isinstance(PointData, dict)):
        print("Your data is a dictionnary of dataframes, let me create a PointTool object that contains all of these.")
        PointData = pd.concat(PointData)
        PointData = LSDMap_PD.LSDMap_PointData(PointData,data_type ="pandas", PANDEX = True)


    plt.clf()
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
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


    diff = PointData.QueryData("diff")
    ratio = PointData.QueryData("ratio")

    PointData.selectValue('diff_outlier',value = True, operator = "==")
    PointData.selectValue('ratio_outlier',value = True, operator = "==")

    diffout = PointData.QueryData("diff")
    ratioout = PointData.QueryData("ratio")
    elevation = PointData.QueryData("elevation")
    if(log_data):
        print("I am logging the data")

        diff = np.log10(diff)
        ratio = np.log10(ratio)
        if(isinstance(diffout, list) == False):
            diffout = [diffout]
            ratioout = [ ratioout]

        for i in range(len(diffout)):

            diffout[i]= np.log10(diffout[i])

            ratioout[i]= np.log10(ratioout[i])
    print("Now plotting the outliers vs the non-outliers")
    ax.scatter(diff,ratio, s=0.5, lw = 0, c = 'gray')
    ax.scatter(diffout,ratioout, s = 0.5,c = elevation,lw = 0)
    ax.set_xlabel("Diff")
    ax.set_ylabel("Ratio")

    plt.savefig(DataDirectory+saveName+save_fmt,dpi=500)
    print("Your figure is " +DataDirectory+saveName+save_fmt)





############################ Mapping scripts ######################################




def map_custom():
    """
    Testing function to plot custom maps before creating real function for mapping routines

    Args:
        Yes.
    returns:
        No.
    Author:
        BG
    """
    ###### Parameters ######
    Directory = "/home/s1675537/PhD/DataStoreBoris/GIS/Data/Carpathian/knickpoint/" # reading directory
    wDirectory = Directory # writing directory
    Base_file = "Buzau" # It will be the cabkground raster. Each other raster you want to drap on it will be cropped to its extents including nodata
    csv_file = Directory + "test_sign_m_chi.csv" # Name of your point file, add a similar line with different name if you have more than one point file
    DrapeRasterName = "Buzau_hs.bil" # if you want to drap a raster on your background one. Just add a similar line in case you want another raster to drap and so on
    wname = "sign_test" # name of your output file
    dpi = 500 # Quality of your output image, don't exceed 900
    fig_size_inches = 7 # Figure size in Inches

    ##### Now we can load and plot the data

    BackgroundRasterName = Base_file + ".bil" # Ignore this line
    thisPointData = LSDP.LSDMap_PointData(csv_file, PANDEX = True) # Load the point file #1, add a similar line with different name if you have more than one point file.

    plt.clf() # Ignore this line

    MF = MapFigure(BackgroundRasterName, Directory,coord_type="UTM_km") # load the background raster

    MF.add_drape_image(DrapeRasterName,Directory, # Calling the function will add a drapped raster on the top of the background one
                        colourmap = "gray", # colormap used for this raster, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                        alpha=0.5, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                        show_colourbar = False, # Well, this one is explicit I think
                        colorbarlabel = "Colourbar") # Name of your Colourbar, it might bug though



    MF.add_point_data( thisPointData, # this function plot the requested point file using the lat/long column in the csv file
                       column_for_plotting = "m_chi_sign",  # Column used to color the data
                       this_colourmap = "cubehelix", # Colormap used, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                       colorbarlabel = "Colourbar", # Label
                       scale_points = False, # All the point will have the same size if False
                       column_for_scaling = "None", # If scale point True, you can scale the size of your points using one of the columns
                       scaled_data_in_log = False, # If scale point True, you can log the scaling
                       max_point_size = 5, # max size if scale point True again
                       min_point_size = 0.5, # You should be able to guess that one now
                       coulor_log = False, # do you want a log scale for your colorbar ?
                       coulor_manual_scale = [], #Do you want to manually limit the scale of your colorbar? if not let is false
                       manual_size = 0.5, # If none of above is choosen but you want to put another value than 0.5 to scale your point
                       alpha = 1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                       minimum_log_scale_cut_off = -10) # you probably won't need this

    ImageName = wDirectory+str(int(clock.time()))+wname+".png" # Ignore this
    ax_style = "Normal" # Ignore this
    MF.save_fig(fig_width_inches = fig_size_inches, FigFileName = ImageName, axis_style = ax_style, Fig_dpi = dpi) # Save the figure

def map_knickpoint_sign(PointData, DataDirectory, Raster_base_name, HS_name = "none",Time_in_name = False, river_network = "none", saveName = "none", size = 2, outliers = 'none'):
    """
    Will create a map of the knickpoint simply colored by sign.

    Args:
        PointData (PointTools object)
        DataDirectory (str): directory where the data will be saved and loaded.
        Raster_base_name (str): Base name of your files without the .bil
        HS_name (str): name of your Hillshade file, by default baseName + _hs.bil like LSDTT create it
        Time_in_name (bool): Option to add timing info in the nae of the figure. Can be useful if you test loads of parameters and you want to be sure that your files names are different (but awful).
    returns:
        No, but creates a map named map_knickpoint_sign.png
    Author:
        BG
    """
    ###### Parameters ######
    if(isinstance(PointData, dict)):
        print("Your data is a dictionnary of dataframes, let me create a PointTool object that contains all of these.")
        PointData = pd.concat(PointData)
        PointData = LSDMap_PD.LSDMap_PointData(PointData,data_type ="pandas", PANDEX = False)
    if(outliers != 'none' ):
        PointData.selectValue(outliers, operator = "==", value = True)
        print PointData

    Directory = DataDirectory # reading directory
    wDirectory = Directory # writing directory
    Base_file = Raster_base_name # It will be the cabkground raster. Each other raster you want to drap on it will be cropped to its extents including nodata
    if(saveName == "none"):
        wname = "map_knickpoint_sign" # name of your output file
    else:
        wname = saveName
    dpi = 500 # Quality of your output image, don't exceed 900
    fig_size_inches = 7 # Figure size in Inches
    if(HS_name == "none"):
        HS_name = Raster_base_name+("_hs.bil")
    DrapeRasterName = HS_name # if you want to drap a raster on your background one. Just add a similar line in case you want another raster to drap and so on

    ##### Now we can load and plot the data

    BackgroundRasterName = Base_file + ".bil" # Ignore this line
    plt.clf() # Ignore this line

    MF = MapFigure(BackgroundRasterName, Directory,coord_type="UTM_km", NFF_opti = True) # load the background raster

    MF.add_drape_image(DrapeRasterName,Directory, # Calling the function will add a drapped raster on the top of the background one
                        colourmap = "gray", # colormap used for this raster, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                        alpha=0.5, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                        show_colourbar = False, # Well, this one is explicit I think
                        colorbarlabel = "Colourbar", # Name of your Colourbar, it might bug though
                        NFF_opti = True)


    if(isinstance(river_network,LSDP.LSDMap_PointData)):
        MF.add_point_data( river_network, # this function plot the requested point file using the lat/long column in the csv file
                           column_for_plotting = "none",  # Column used to color the data
                           this_colourmap = "cubehelix", # Colormap used, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                           colorbarlabel = "Colourbar", # Label
                           scale_points = False, # All the point will have the same size if False
                           column_for_scaling = "None", # If scale point True, you can scale the size of your points using one of the columns
                           scaled_data_in_log = False, # If scale point True, you can log the scaling
                           max_point_size = 5, # max size if scale point True again
                           min_point_size = 0.5, # You should be able to guess that one now
                           coulor_log = False, # do you want a log scale for your colorbar ?
                           coulor_manual_scale = [], #Do you want to manually limit the scale of your colorbar? if not let is false
                           manual_size = 0.5, # If none of above is choosen but you want to put another value than 0.5 to scale your point
                           alpha = 1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                           minimum_log_scale_cut_off = -10) # you probably won't need this

    MF.add_point_data( PointData, # this function plot the requested point file using the lat/long column in the csv file
                       column_for_plotting = "sign",  # Column used to color the data
                       this_colourmap = "cubehelix", # Colormap used, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                       colorbarlabel = "Colourbar", # Label
                       scale_points = False, # All the point will have the same size if False
                       column_for_scaling = "None", # If scale point True, you can scale the size of your points using one of the columns
                       scaled_data_in_log = False, # If scale point True, you can log the scaling
                       max_point_size = 5, # max size if scale point True again
                       min_point_size = 0.5, # You should be able to guess that one now
                       coulor_log = False, # do you want a log scale for your colorbar ?
                       coulor_manual_scale = [], #Do you want to manually limit the scale of your colorbar? if not let is false
                       manual_size = size, # If none of above is choosen but you want to put another value than 0.5 to scale your point
                       alpha = 1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                       minimum_log_scale_cut_off = -10 # you probably won't need this
                      )
    if(Time_in_name):
        ImageName = wDirectory+str(int(clock.time()))+wname+".png" # Ignore this
    else:
        ImageName = wDirectory+wname+".png" # Ignore this
    ax_style = "Normal" # Ignore this
    MF.save_fig(fig_width_inches = fig_size_inches, FigFileName = ImageName, axis_style = ax_style, Fig_dpi = dpi) # Save the figure

def map_knickpoint_diff_sized_colored_ratio(PointData, DataDirectory, Raster_base_name, HS_name = "none",Time_in_name = False, river_network = "none", saveName = "none", log = False):
    """
    Will create a map of the knickpoint simply colored by sign.

    Args:
        PointData (PointTools object)
        DataDirectory (str): directory where the data will be saved and loaded.
        Raster_base_name (str): Base name of your files without the .bil
        HS_name (str): name of your Hillshade file, by default baseName + _hs.bil like LSDTT create it
        Time_in_name (bool): Option to add timing info in the nae of the figure. Can be useful if you test loads of parameters and you want to be sure that your files names are different (but awful).
    returns:
        No, but creates a map named map_knickpoint_sign.png
    Author:
        BG
    """
    ###### Parameters ######
    Directory = DataDirectory # reading directory
    wDirectory = Directory # writing directory
    Base_file = Raster_base_name # It will be the cabkground raster. Each other raster you want to drap on it will be cropped to its extents including nodata
    if(saveName == "none"):
        wname = "map_knickpoint_diff_sized_colored_ratio" # name of your output file
    else:
        wname = saveName
    dpi = 500 # Quality of your output image, don't exceed 900
    fig_size_inches = 7 # Figure size in Inches
    if(HS_name == "none"):
        HS_name = Raster_base_name+("_hs.bil")
    DrapeRasterName = HS_name # if you want to drap a raster on your background one. Just add a similar line in case you want another raster to drap and so on

    ##### Now we can load and plot the data

    BackgroundRasterName = Base_file + ".bil" # Ignore this line
    plt.clf() # Ignore this line

    MF = MapFigure(BackgroundRasterName, Directory,coord_type="UTM_km") # load the background raster

    MF.add_drape_image(DrapeRasterName,Directory, # Calling the function will add a drapped raster on the top of the background one
                        colourmap = "gray", # colormap used for this raster, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                        alpha=0.5, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                        show_colourbar = True, # Well, this one is explicit I think
                        colorbarlabel = "Colourbar") # Name of your Colourbar, it might bug though



    if(isinstance(river_network,LSDP.LSDMap_PointData)):
        MF.add_point_data( river_network, # this function plot the requested point file using the lat/long column in the csv file
                           column_for_plotting = "none",  # Column used to color the data
                           this_colourmap = "cubehelix", # Colormap used, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                           colorbarlabel = "Colourbar", # Label
                           scale_points = False, # All the point will have the same size if False
                           column_for_scaling = "None", # If scale point True, you can scale the size of your points using one of the columns
                           scaled_data_in_log = False, # If scale point True, you can log the scaling
                           max_point_size = 5, # max size if scale point True again
                           min_point_size = 0.5, # You should be able to guess that one now
                           coulor_log = False, # do you want a log scale for your colorbar ?
                           coulor_manual_scale = [], #Do you want to manually limit the scale of your colorbar? if not let is false
                           manual_size = 0.1, # If none of above is choosen but you want to put another value than 0.5 to scale your point
                           alpha = 1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                           minimum_log_scale_cut_off = -10) # you probably won't need this

    MF.add_point_data( PointData, # this function plot the requested point file using the lat/long column in the csv file
                       column_for_plotting = "ratio",  # Column used to color the data
                       this_colourmap = "cubehelix", # Colormap used, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                       colorbarlabel = "Ratio", # Label
                       scale_points = True, # All the point will have the same size if False
                       column_for_scaling = "diff", # If scale point True, you can scale the size of your points using one of the columns
                       scaled_data_in_log = log, # If scale point True, you can log the scaling
                       max_point_size = 20, # max size if scale point True again
                       min_point_size = 0.5, # You should be able to guess that one now
                       coulor_log = log, # do you want a log scale for your colorbar ?
                       coulor_manual_scale = [], #Do you want to manually limit the scale of your colorbar? if not let is false
                       manual_size = 0.5, # If none of above is choosen but you want to put another value than 0.5 to scale your point
                       alpha = 1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                       minimum_log_scale_cut_off = -10) # you probably won't need this

    if(Time_in_name):
        ImageName = wDirectory+str(int(clock.time()))+wname+".png" # Ignore this
    else:
        ImageName = wDirectory+wname+".png" # Ignore this
    ax_style = "Normal" # Ignore this
    MF.save_fig(fig_width_inches = fig_size_inches, FigFileName = ImageName, axis_style = ax_style, Fig_dpi = dpi) # Save the figure






def DEPRECATED_map_knickpoint_standard(PointData, DataDirectory, Raster_base_name, HS_name = "none",Time_in_name = False, river_network = "none", saveName = "none", log = False):
    """
    Will create a map of the knickpoint sized by their delta and colored by size.

    Args:
        PointData (PointTools object)
        DataDirectory (str): directory where the data will be saved and loaded.
        Raster_base_name (str): Base name of your files without the .bil
        HS_name (str): name of your Hillshade file, by default baseName + _hs.bil like LSDTT create it
        Time_in_name (bool): Option to add timing info in the nae of the figure. Can be useful if you test loads of parameters and you want to be sure that your files names are different (but awful).
    returns:
        No, but creates a map named map_knickpoint_sign.png
    Author:
        BG
    """
    ###### Parameters ######
    Directory = DataDirectory # reading directory
    wDirectory = Directory # writing directory
    Base_file = Raster_base_name # It will be the cabkground raster. Each other raster you want to drap on it will be cropped to its extents including nodata
    if(saveName == "none"):
        wname = "map_knickpoint_std" # name of your output file
    else:
        wname = saveName
    dpi = 500 # Quality of your output image, don't exceed 900
    fig_size_inches = 7 # Figure size in Inches
    if(HS_name == "none"):
        HS_name = Raster_base_name+("_hs.bil")
    DrapeRasterName = HS_name # if you want to drap a raster on your background one. Just add a similar line in case you want another raster to drap and so on

    ##### Now we can load and plot the data

    BackgroundRasterName = Base_file + ".bil" # Ignore this line
    plt.clf() # Ignore this line

    MF = MapFigure(BackgroundRasterName, Directory,coord_type="UTM_km") # load the background raster

    MF.add_drape_image(DrapeRasterName,Directory, # Calling the function will add a drapped raster on the top of the background one
                        colourmap = "gray", # colormap used for this raster, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                        alpha=0.5, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                        show_colourbar = True, # Well, this one is explicit I think
                        colorbarlabel = "Colourbar") # Name of your Colourbar, it might bug though



    if(isinstance(river_network,LSDP.LSDMap_PointData)):
        MF.add_point_data( river_network, # this function plot the requested point file using the lat/long column in the csv file
                           column_for_plotting = "none",  # Column used to color the data
                           this_colourmap = "cubehelix", # Colormap used, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                           colorbarlabel = "Colourbar", # Label
                           scale_points = False, # All the point will have the same size if False
                           column_for_scaling = "None", # If scale point True, you can scale the size of your points using one of the columns
                           scaled_data_in_log = False, # If scale point True, you can log the scaling
                           max_point_size = 5, # max size if scale point True again
                           min_point_size = 0.5, # You should be able to guess that one now
                           coulor_log = False, # do you want a log scale for your colorbar ?
                           coulor_manual_scale = [], #Do you want to manually limit the scale of your colorbar? if not let is false
                           manual_size = 0.1, # If none of above is choosen but you want to put another value than 0.5 to scale your point
                           alpha = 1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                           minimum_log_scale_cut_off = -10) # you probably won't need this

    MF.add_point_data( PointData, # this function plot the requested point file using the lat/long column in the csv file
                       column_for_plotting = "sign",  # Column used to color the data
                       this_colourmap = "cubehelix", # Colormap used, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                       colorbarlabel = "diff", # Label
                       scale_points = True, # All the point will have the same size if False
                       column_for_scaling = "diff", # If scale point True, you can scale the size of your points using one of the columns
                       scaled_data_in_log = log, # If scale point True, you can log the scaling
                       max_point_size = 20, # max size if scale point True again
                       min_point_size = 5, # You should be able to guess that one now
                       coulor_log = log, # do you want a log scale for your colorbar ?
                       coulor_manual_scale = [], #Do you want to manually limit the scale of your colorbar? if not let is false
                       manual_size = 0.5, # If none of above is choosen but you want to put another value than 0.5 to scale your point
                       alpha = 1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                       minimum_log_scale_cut_off = -10) # you probably won't need this

    if(Time_in_name):
        ImageName = wDirectory+str(int(clock.time()))+wname+".png" # Ignore this
    else:
        ImageName = wDirectory+wname+".png" # Ignore this
    ax_style = "Normal" # Ignore this
    MF.save_fig(fig_width_inches = fig_size_inches, FigFileName = ImageName, axis_style = ax_style, Fig_dpi = dpi) # Save the figure


def plot_pdf_diff_ratio(df, DataDirectory, saveName = "pdf_diff_ratio", save_fmt = ".png", size_format = "ESURF",  xlim =[]):
    """
    Basic plot to have a general view of the knickpoints: flow distance against ratio and diff colored by elevation

    Args:
        PointData: A PointData object
        DataDirectory: Where the data is saved
        saveName: save name

    returns:
        Nothing, sorry.
    Author: BG
    """
    plt.clf()
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(2, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(2, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(2, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35


    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax1 = fig.add_subplot(gs[10:50,10:95])
    ax2 = fig.add_subplot(gs[50:100,10:95])


    ax1.scatter(df["ratio"],norm.pdf(df["ratio"]),lw =0, s = 1, c = "red")
    ax1.set_ylabel("Ratio")
    ax1.tick_params(axis = 'x', length = 0, width = 0, labelsize = 0)
    ax1.spines['bottom'].set_visible(False)
    ax2.scatter(df["diff"],norm.pdf(df["diff"]),lw =0, s = 1, c = "red")
    ax2.set_ylabel("Diff")
    ax2.set_xlabel("PDF")


    #'###### Setting the limits
    if(xlim != []):
        ax2.set_xlim(xlim[0],xlim[1])
        ax1.set_xlim(xlim[0],xlim[1])



    #ax2.tick_params(axis = 'x', labelsize = 6)
    #ax1.set_xticks([4,5,6,7,8,9,10])
    #ax2.set_xticks([4,5,6,7,8,9,10])
    #ax2.set_xticklabels([ur"$10^{4}$",ur"$10^{5}$",ur"$10^{6}$",ur"$10^{7}$",ur"$10^{8}$",ur"$10^{9}$",ur"$10^{10}$"])

    plt.savefig(DataDirectory+saveName+save_fmt,dpi=500)

def violin_by_bin(ldf, DataDirectory, saveName = "Violin", column = "elevation", size_format = "ESURF"):

    """
    Will plot violin from a list of bins. NOT READY YET.

    Author: BG

    matplotlib description:
        Violin plots are similar to histograms and box plots in that they show
    an abstract representation of the probability distribution of the
    sample. Rather than showing counts of data points that fall into bins
    or order statistics, violin plots use kernel density estimation (KDE) to
    compute an empirical distribution of the sample. That computation
    is controlled by several parameters. This example demonstrates how to
    modify the number of points at which the KDE is evaluated (``points``)
    and how to modify the band-width of the KDE (``bw_method``).

    For more information on violin plots and KDE, the scikit-learn docs
    have a great section: http://scikit-learn.org/stable/modules/density.html
    """

    plt.clf()
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(2, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(2, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(2, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35


    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
    ax1 = fig.add_subplot(gs[10:50,10:95])
    ax2 = fig.add_subplot(gs[50:100,10:95])


    ax2.set_ylabel("Ratio")
    ax1.set_ylabel("Diff")
    ax2.set_xlabel(column)
    plt.savefig(DataDirectory+saveName+save_fmt,dpi=500)


def pdf_from_bin(ldf, DataDirectory, saveName = "BasicPDF_", column = "elevation", size_format = "ESURF" ):

    """
    Produce some simple pdf plots from a list of pandas dataframe.

    Arg:

    Returns: nothing, but produce a plot.

    Author: BG
    """

    for inch in ldf:
        plt.clf()
        label_size = 10
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['arial']
        rcParams['font.size'] = label_size

        # make a figure
        if size_format == "geomorphology":
            fig = plt.figure(2, facecolor='white',figsize=(6.25,3.5))
            l_pad = -40
        elif size_format == "big":
            fig = plt.figure(2, facecolor='white',figsize=(16,9))
            l_pad = -50
        else:
            fig = plt.figure(2, facecolor='white',figsize=(4.92126,3.5))
            l_pad = -35



        gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=1.0,top=1.0)
        ax1 = fig.add_subplot(gs[10:50,10:95])
        ax2 = fig.add_subplot(gs[50:100,10:95])

        ax2.scatter(ldf[inch]["diff"],norm.pdf(ldf[inch]["diff"]), s = 1.5, lw = 0)
        ax1.scatter(ldf[inch]["ratio"],norm.pdf(ldf[inch]["ratio"]), s = 1.5, lw = 0)

        ax2.set_ylabel("PDF (Diff)")
        ax1.set_ylabel("PDF (Ratio)")
        ax2.set_xlabel("Diff/ratio binned by " + column + "_" + inch)
        plt.savefig(DataDirectory+saveName+inch+"_"+column+".png",dpi=500)


def pdf_from_bin_one_col(ldf, DataDirectory, saveName = "BasicPDF_", column = "elevation", size_format = "ESURF", pdf_col = "diff", combine_diff_sign = False, argsort = False ):

    """
    Produce some simple pdf plots from a dict of pandas dataframe.

    Arg:

    Returns: nothing, but produce a plot.

    Author: BG
    """




    for inch in ldf:
        plt.clf()
        label_size = 10
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['arial']
        rcParams['font.size'] = label_size



        if(combine_diff_sign):
            ldf[inch]["diff"][ldf[inch]["sign"] == -1] = -ldf[inch]["diff"][ldf[inch]["sign"] == -1]

        data = np.array(ldf[inch][pdf_col].values)


        # make a figure
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
        ax1 = fig.add_subplot(gs[10:100,5:95])

        print(data.shape)
        if(data.shape[0]>0):
            ax1.hist(data, 100, normed=1, facecolor='green', alpha=0.75)






        ax1.set_ylabel("PDF")
        ax1.set_xlabel("elevation by " + pdf_col)
        ax1.set_xlim(-100,100)
        plt.savefig(DataDirectory+saveName+inch+"_"+column+".png",dpi=500)

def plot_2d_density_map(dataframe, DataDirectory, columns = ["drainage area", "diff"], bin = 50,   saveName = "BasicPDF_", size_format = "ESURF",):

    """
    Plots a 2d histogram or density plot or heatmap depending how you name it of two variables.

    Args:
        dataframe: a Pandas dataframe
        columns (list of str): The x,y columns to plot
        bin (int): number of bins

    returns:
        Nothing yet, plot a figure.
    """

    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.5))
        l_pad = -35

    X_data = dataframe[columns[0]]
    Y_data = dataframe[columns[1]]




if __name__ == "__main__":
    print("Do not use this file as a script, refer to LSDMT documentation for instructions")




#
