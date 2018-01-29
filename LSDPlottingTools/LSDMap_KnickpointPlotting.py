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
from matplotlib.patches import Rectangle
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

    ######################################################################################################################
    ####################### A first set of general functions to prepare the data/figures #################################
    ######################################################################################################################


                                            #                  .
                                            #              /\ /l
                                            #             ((.Y(!
                                            #              \ |/
                                            #              /  6~6,
                                            #              \ _    +-.
                                            #               \`-=--^-'
                                            #                \ \
                                            #               _/  \
                                            #              (  .  Y
                                            #             /"\ `--^--v--.
                                            #            / _ `--"T~\/~\/
                                            #           / " ~\.  !
                                            #     _    Y      Y./'
                                            #    Y^|   |      |~~7
                                            #    | l   |     / ./'
                                            #    | `L  | Y .^/~T
                                            #    |  l  ! | |/| |          
                                            #    | .`\/' | Y | !
                                            #    l  "~   j l j_L______
                                            #     \,____{ __"~ __ ,\_,\_
                                            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    def get_fig_right_size(self, size = "esurf", n_axis =1 ,facecolor = 'white'):
        """
        return a matplotlib fig object presized
        param:
            size = size code (esurf,geomorphology,big or tuple/array of size)
            n_axis: number of axis
        """

        # Cheching the type of input

        if (isinstance(size,str)):
            if size.lower() == "geomorphology":
                fig = plt.figure(n_axis, facecolor = facecolor, figsize=(6.25,3.5))            
            elif size.lower() == "big":
                fig = plt.figure(n_axis, facecolor = facecolor, figsize=(16,9))            
            elif size.lower() == "esurf":
                fig = plt.figure(n_axis, facecolor = facecolor, figsize=(4.92126,3.5))
            else:
                print("I did not understood your format input (%s), I am defaulting to esurf." %(size))
                fig = plt.figure(n_axis, facecolor = facecolor, figsize=(4.92126,3.5))
        if ((isinstance(size,tuple)) or (isinstance(size,list))):
            if len(size) == 2:
                fig = plt.figure(n_axis, facecolor = facecolor, figsize=(size[0], size[1]))
            else:
                print("I did not understood your format input (%s), I am defaulting to esurf." %(size))
                fig = plt.figure(n_axis, facecolor = facecolor, dpi = 600) 

        return fig



    ######################################################################################################################
    ########################################### The plotting figures #####################################################
    ######################################################################################################################

                    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$$$$$$**$$$$$$$$$**$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$$$$"   ^$$$$$$F    *$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$$$     z$$$$$$L    ^$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$$$    e$$$$$$$$$e  J$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$$$eee$$$$$$$$$$$$$e$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$$b$$$$$$$$$$$$$$$$$$*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$)$$$$P"e^$$$F$r*$$$$F"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$d$$$$  "z$$$$"  $$$$%  $3$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$*"""*$$$  .$$$$$$ z$$$*   ^$e*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$"     *$$ee$$$$$$$$$$*"     $$$C$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$.      "***$$"*"$$""        $$$$e*$$$$$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$b          "$b.$$"          $$$$$b"$$$$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$c.         """            $$$$$$$^$$$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$e..                     $$$$$$$$^$$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$$$$$$eeee..            J$$$$$$$$b"$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$r          z$$$$$$$$$$r$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"         z$$$$$**$$$$$^$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$$$$$$$$$*"          z$$$P"   ^*$$$ $$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$$$$$$*"           .d$$$$       $$$ $$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$$$$$"           .e$$$$$F       3$$ $$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$$$$$.         .d$$$$$$$         $PJ$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$$$$$$eeeeeeed$*""""**""         $\$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$                  $d$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$.                 $$$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$e.              d$$$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$eeeeeee$$$$$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


    def print_ksn_profile(self,size = "big", format = "png", x_axis = "chi", knickpoint = True, title = "none", label_size = 8, facecolor = 'white',
        size_of_ksn = 4, legend = True, size_of_TVD_ksn = 3):

        """
        print a plot for each source keys selected with ksn value function to Chi or Flow Distance.
        param:
            size: the size of the figure, default big.
            format: format of the output: "png", "svg" or "show".
            x_axis: The coordinates to print the data: "chi" or "flow distance"
            knickpoint: True or False to display it or not.
            title: "none" for no title, "auto" for displaying the source key.
            size_of_ksn: size of the m_chi (ksn) points before processing
            legend: if True, will plot the legend
        Author: B.G. - 25/01/2018
        """

        # check if a directory exists for the chi plots. If not then make it.
        out_directory = self.fpath+'river_plots/'
        if not os.path.isdir(out_directory):
            print("I am creating the river_plot/ directory to save your figures")
            os.makedirs(out_directory)

        
        # Set up fonts for plots
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['Liberation Sans'] # Liberation Sans is a free alternative to Arial. Albeit being quite universal, Arial is propietary. #PRAISE_FREE_AND_OPENSOURCE
        rcParams['font.size'] = label_size


        for sources in self.df_SK["source_key"].unique():

            # Select the data
            this_df_SK = self.df_SK[self.df_SK["source_key"] == sources]
            this_df_kp = self.df_kp[self.df_kp["source_key"] == sources]
            this_df_kp_raw = self.df_kp_raw[self.df_kp_raw["source_key"] == sources]
            this_df_river = self.df_river[self.df_river["source_key"] == sources]

            
            # Create a figure with required dimensions
            n_axis = 1
            fig = self.get_fig_right_size(size = size, n_axis =1 , facecolor = facecolor)

            # create the axis using the gridspec method
            gs = plt.GridSpec(100,100,bottom=0.15,left=0.10,right=0.90,top=0.95)
            ax2 = fig.add_subplot(gs[0:100,0:100], facecolor = "None")
            gs = plt.GridSpec(100,100,bottom=0.15,left=0.10,right=0.90,top=0.95)
            ax1 = fig.add_subplot(gs[0:100,0:100], facecolor = "None")


            # plotting the ksn
            ## not processed (in the back and quite transparent)
            ax1.scatter(this_df_river[x_axis],this_df_river["m_chi"], s = size_of_ksn, c = "r", lw =0, alpha = 0.3, label = "ksn (before TVD)")
            ax1.scatter(this_df_river[x_axis],this_df_river["TVD_ksn"], s = size_of_TVD_ksn, c ="k", lw =0, alpha = 1, label = "ksn (TVD)")
            ## Getting the extents of this first plot to apply it to the knickpoint one
            this_xlim = ax1.get_xlim()

            # Label
            if(x_axis == "chi"):
                xlab = r"$\chi$"
            elif(x_axis == "flow_distance"):
                xlab = "Glow distance (m)"
            ax1.set_xlabel(xlab)
            ax1.set_ylabel(r"$k_{sn}$")



            if(knickpoint):
                this_df_kp_pos = this_df_kp[this_df_kp["delta_ksn"]>0]
                this_df_kp_neg = this_df_kp[this_df_kp["delta_ksn"]<0]
                ax2.scatter(this_df_kp_pos[x_axis], this_df_kp_pos["delta_ksn"], marker = "s", s = 3, c = "#E79A00")
                ax2.scatter(this_df_kp_neg[x_axis], this_df_kp_neg["delta_ksn"], marker = "s", s = 3, c = "#2939FF")

            # Adapting hte extents 
            ax2.set_xlim(this_xlim)

            # Title
            if(title.lower() == "auto"):
                this_title = "source %s" %(sources)
            elif(title.lower() != "none"):
                this_title = title

            if(title.lower() != "none"):
                extra = ax1.add_patch(Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0, label = this_title))



            # Legend
            ax1.legend(loc = 1) # 1 = upper right
            ax2.xaxis.set_visible(False)
            if(knickpoint):
                ax2.yaxis.set_ticks_position("right")
                ax2.yaxis.set_label_position("right")
                ax2.set_ylabel(r"$Knickpoint \/ \Delta k_{sn}$") # the \/ add a space in between, the mathematical expression compiling can screw this a bit
            else:
                ax2.yaxis.set_visible(False)

            # Saving the figure
            plt.savefig(out_directory + self.fprefix+"_ksn_source%s_%s.%s"%(sources,x_axis,format), dpi = 500)
            plt.clf()
            # switching to the next figure



 

        # End of this function

    def print_river_profile(self,size = "big", format = "png", x_axis = "chi", knickpoint = True, title = "none", label_size = 8, facecolor = 'white',
        size_of_ksn = 4, legend = True, size_of_TVD_ksn = 3):
        print("TODO")

    def print_map_of_kp(self,size = "big", format = "png", x_axis = "chi", knickpoint = True, title = "none", label_size = 8, facecolor = 'white',
        size_of_ksn = 4, legend = True, size_of_TVD_ksn = 3):
        print("TODO")



#
#
#
#
#
#
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# End of the file, I am adding these blank lines before because my text editor delete them when I save but I hate working on the extreme bottom of my screen


#                      ____,------------------,______
#                  ___/    \            /            \_____
#               __/         \__________/              \___ \___
# ,^------.____/\           /          \              /   `----\_
# | (O))      /  \_________/            \____________/         \ \
# \_____,--' /   /         \            /            \          \ \
#   \___,---|___/_______,---`----------'----,_________\__________\_\
#             /  :__________________________/  :___________________/
#            /   :          /   :          /   :          /   :
#           /    :         /    :         /    :         /    :
#       (~~~     )     (~~~     )     (~~~     )     (~~~     )
#        ~~~~~~~~       ~~~~~~~~       ~~~~~~~~       ~~~~~~~~
