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


    def __init__(self, fpath,fprefix, basin_key = [], source_key = [], min_length = 0, cut_off_val = 0):

        print("Let me first preprocess and check your files")

        # Loading the attributes
        self.fpath = fpath # the path of your file : /home/your/path/
        self.fprefix = fprefix # the common prefix of all your files

        # Loading the files

        print("Loading the knickpoint-related files")
        
        try:
            self.df_river = Helper.ReadMChiSegCSV(self.fpath, self.fprefix, type = "knickpoint") # Contains the river info
            self.df_rivraw = Helper.ReadMChiSegCSV(self.fpath, self.fprefix, type = "knickpoint") # Contains the river info (will not be thinned by your selection choices)
            self.df_kp_raw = Helper.ReadKnickpointCSV(self.fpath, self.fprefix, ftype = "raw") # Contains the raw knickpint info (before TVD or else) -> Debugging purposes
            self.df_kp = Helper.ReadKnickpointCSV(self.fpath, self.fprefix) # Contains the knickpoint location and informations
            self.df_SK = Helper.readSKKPstats(self.fpath, self.fprefix) # Contains few metrics per river keys
        except IOError:
            print("I didnae find your knickpoint related files make sure that:")
            print("- You ran the knickpoint analysis")
            print("- Your path is good and finishing by '/' e.g. /home/name/kazakhstan/ ")

            quit()


        print("Managing the data:")

        if(cut_off_val > 0):
            print("I am removing the knickpoints with a magnitude below %s" %(cut_off_val))
            self.df_kp = self.df_kp[self.df_kp["delta_ksn"].abs() >= cut_off_val]


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

        


        # dealing with no datas

        self.df_river["m_chi"][self.df_river["m_chi"] == -9999] = 0

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

    def get_figwidth_right_size(self, size = "esurf"):
        """
        return a matplotlib fig object presized
        param:
            size = size code (esurf,geomorphology,big or tuple/array of size)
            n_axis: number of axis
        """

        # Cheching the type of input

        if (isinstance(size,str)):
            if size.lower() == "geomorphology":
                wsize = 6.25            
            elif size.lower() == "big":
                wsize = 16           
            elif size.lower() == "esurf":
                wsize = 4.92126
            else:
                print("I did not understood your format input (%s), I am defaulting to esurf." %(size))
                wsize = 4.92126
        if ((isinstance(size,tuple)) or (isinstance(size,list))):
            if len(size) == 2:
                wsize = size[0]
            else:
                print("I did not understood your format input (%s), I am defaulting to esurf." %(size))
                wsize = 4.92126

        return wsize



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
            this_df_kp = this_df_kp[this_df_kp["out"] == 1]
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
                xlab = "Distance from the outlet (m)"
            ax1.set_xlabel(xlab)
            ax1.set_ylabel(r"$k_{sn}$")



            if(knickpoint):
                this_df_kp_pos = this_df_kp[this_df_kp["delta_ksn"]>0]
                this_df_kp_neg = this_df_kp[this_df_kp["delta_ksn"]<0]
                ax2.scatter(this_df_kp_pos[x_axis], this_df_kp_pos["delta_ksn"], marker = "s", s = 5, c = "#E79A00")
                ax2.scatter(this_df_kp_neg[x_axis], this_df_kp_neg["delta_ksn"], marker = "s", s = 5, c = "#2939FF")

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
            ax1.legend(loc = 0) # 1 = upper right 0 - best choice
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
        size_of_river = 0.5, legend = True, size_of_TVD_ksn = 3):

        """
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
            # this_df_kp = this_df_kp[this_df_kp["out"] == 1]
            this_df_kp_pos = this_df_kp[this_df_kp["sign"] == 1]
            this_df_kp_neg = this_df_kp[this_df_kp["sign"] == -1]

            this_df_kp_raw = self.df_kp_raw[self.df_kp_raw["source_key"] == sources]
            this_df_river = self.df_river[self.df_river["source_key"] == sources]

            
            # Create a figure with required dimensions
            n_axis = 1
            fig = self.get_fig_right_size(size = size, n_axis =1 , facecolor = facecolor)

            gs = plt.GridSpec(100,100,bottom=0.15, left=0.10, right=0.95, top=0.95)
            ax1 = fig.add_subplot(gs[0:100,0:100], facecolor = "None")

            #plot the long/Chi profile
            cb1 = ax1.scatter(this_df_river[x_axis], this_df_river["elevation"], s = size_of_river, c = "blue")
            cb1 = ax1.scatter(this_df_river[x_axis], this_df_river["segmented_elevation"], s = size_of_river/2, c = "k", alpha = 0.5)

            ax1.scatter(this_df_kp_pos[x_axis], this_df_kp_pos["elevation"], s = 30, lw = 0, marker = "^", c = this_df_kp_pos["delta_ksn"], cmap = "RdBu_r", vmin = this_df_kp["delta_ksn"].min(), vmax = this_df_kp["delta_ksn"].max() , alpha = 0.95)
            ax1.scatter(this_df_kp_neg[x_axis], this_df_kp_neg["elevation"], s = 30, lw = 0, marker = "v", c = this_df_kp_neg["delta_ksn"], cmap = "RdBu_r", vmin = this_df_kp["delta_ksn"].min(), vmax = this_df_kp["delta_ksn"].max() , alpha = 0.95)
            ax1.scatter(this_df_kp_pos[x_axis], this_df_kp_pos["elevation"], s = 30, lw = 0.5, marker = "^", facecolor = "none", edgecolor = "k", cmap = "RdBu_r", vmin = this_df_kp["delta_ksn"].min(), vmax = this_df_kp["delta_ksn"].max() , alpha = 0.95)
            ax1.scatter(this_df_kp_neg[x_axis], this_df_kp_neg["elevation"], s = 30, lw = 0.5, marker = "v", facecolor = "none", edgecolor = "k", cmap = "RdBu_r", vmin = this_df_kp["delta_ksn"].min(), vmax = this_df_kp["delta_ksn"].max() , alpha = 0.95)

            cb2 = ax1.scatter(this_df_kp[x_axis], this_df_kp["elevation"], s = 0, marker = "^", c = this_df_kp["delta_ksn"], cmap = "RdBu_r")


            if(x_axis == "chi"):
                ax1.set_xlabel(r"$\chi$")
            else:
                ax1.set_xlabel("Distance from the outlet (m)")

            # Title
            if(title.lower() == "auto"):
                this_title = "source %s" %(sources)
            elif(title.lower() != "none"):
                this_title = title

            if(title.lower() != "none"):
                extra = ax1.add_patch(Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0, label = this_title))

            ax1.legend([extra],[this_title], loc = 0) # 1 = upper right 0 - best choice

            plt.colorbar(cb2)


            # Saving the figure
            plt.savefig(out_directory + self.fprefix + "_source%s_%s.%s"%(sources,x_axis,format), dpi = 500)
            plt.clf()
            # switching to the next figure







    def print_map_of_kp(self,size = "big", format = "png", black_bg = False, scale_points = False, label_size = 8):

            # check if a directory exists for the chi plots. If not then make it.
        raster_directory = self.fpath+'raster_plots/'
        if not os.path.isdir(raster_directory):
            os.makedirs(raster_directory)

        # Set up fonts for plots
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['Liberation Sans'] # Liberation Sans is a free alternative to Arial. Albeit being quite universal, Arial is propietary. #PRAISE_FREE_AND_OPENSOURCE
        rcParams['font.size'] = label_size

        # set figure sizes based on format
        fig_width_inches = self.get_figwidth_right_size(size = size)

        # get the rasters
        raster_ext = '.bil'
        BackgroundRasterName = self.fprefix+raster_ext
        HillshadeName = self.fprefix+'_hs'+raster_ext
        BasinsName = self.fprefix+'_AllBasins'+raster_ext

        
        # create the map figure
        MF = MapFigure(HillshadeName, self.fpath, coord_type="UTM_km", alpha = 0.7)
        if(black_bg):
            MF.add_drape_image(HillshadeName,self.fpath,colourmap = "gray",alpha=1,colour_min_max = [10000,10001],modify_raster_values=False,old_values=[], new_values=[],NFF_opti = True)

        
        # plot the basin outlines
        Basins = LSDP.GetBasinOutlines(self.fpath, BasinsName)
        MF.plot_polygon_outlines(Basins, linewidth = 0.5)

        # add the channel network without color
        rivnet = LSDP.LSDMap_PointData(self.df_river, data_type = "pandas", PANDEX = True)
        rawriv = LSDP.LSDMap_PointData(self.df_rivraw, data_type = "pandas", PANDEX = True)
        MF.add_point_data(rawriv, show_colourbar="False", scale_points=True,scaled_data_in_log= True, column_for_scaling='drainage_area',alpha=0.1,max_point_size = 0.4,min_point_size = 0.2,zorder=100)
        MF.add_point_data(rivnet, column_for_plotting = "m_chi", show_colourbar="False", scale_points=True,scaled_data_in_log= True, column_for_scaling='drainage_area',alpha=0.1,max_point_size = 1.5,min_point_size = 0.4,zorder=100)

        # add the knickpoints plots

        kp_pos = LSDP.LSDMap_PointData(self.df_kp[self.df_kp["sign"] == 1], data_type = "pandas", PANDEX = True)
        kp_neg = LSDP.LSDMap_PointData(self.df_kp[self.df_kp["sign"] == -1], data_type = "pandas", PANDEX = True)

        MF.add_point_data(kp_pos,this_colourmap = "RdBu_r",colour_manual_scale = [0,self.df_kp["delta_ksn"].max()], marker ="^", column_for_plotting = "delta_ksn", color_abs = True ,show_colourbar=True, colorbarlabel = r'$\Delta k_{sn}$', colourbar_location = "bottom", scale_points = scale_points, scaled_data_in_log= False, column_for_scaling = 'delta_ksn',scale_in_absolute = True ,alpha=1,max_point_size = 15,min_point_size = 1,zorder=200,manual_size = 20)
        MF.add_point_data(kp_neg,this_colourmap = "RdBu_r",colour_manual_scale = [0,self.df_kp["delta_ksn"].max()], marker ="v", column_for_plotting = "delta_ksn", color_abs = True ,show_colourbar="False", scale_points = scale_points, scaled_data_in_log= False, column_for_scaling = 'delta_ksn',scale_in_absolute = True ,alpha=1,max_point_size = 15,min_point_size = 1,zorder=200,manual_size = 20)

        if(black_bg):
            suffix = "dark"
        else:
            suffix = "hs"
        ImageName = raster_directory+self.fprefix+"_ksnkp_map_%s."%(suffix) + format
        MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat = format, Fig_dpi = 500) # Save the figure
        plt.clf()





    def print_histogram(self,size = "big", format = "png", label_size = 8, n_bin = 'auto', facecolor = "white", grid = True):
        """
        This figure print an histogram of the knickpoint repartition for the selected basins/sources

        """


        # check if a directory exists for the chi plots. If not then make it.
        out_directory = self.fpath+'stat_plots/'
        if not os.path.isdir(out_directory):
            print("I am creating the river_plot/ directory to save your figures")
            os.makedirs(out_directory)

        
        # Set up fonts for plots
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['Liberation Sans'] # Liberation Sans is a free alternative to Arial. Albeit being quite universal, Arial is propietary. #PRAISE_FREE_AND_OPENSOURCE
        rcParams['font.size'] = label_size


        # Create a figure with required dimensions
        n_axis = 1
        fig = self.get_fig_right_size(size = size, n_axis =1 , facecolor = facecolor)

        gs = plt.GridSpec(100,100,bottom=0.15, left=0.10, right=0.95, top=0.95)
        ax1 = fig.add_subplot(gs[0:100,0:100], facecolor = "white")

        if(grid):
             ax1.grid(ls = 'dotted', lw = 0.1, c = "k", zorder = 5)

        ax1.hist(self.df_kp["delta_ksn"], bins = n_bin, fc = "#848484", lw = 0.5, edgecolor = "k", zorder = 10)


        ax1.set_xlabel(r"$\Delta k_{sn}$")
        ax1.set_ylabel("n knickpoints")

        plt.savefig(out_directory + self.fprefix + "_kp_hist.%s"%(format), dpi = 500)
        plt.clf()


    def print_box_and_whisker(self,size = "big", format = "png", label_size = 8, binning = 'source_key', facecolor = "white", grid = True):
        """
        This figure print an histogram of the knickpoint repartition for the selected basins/sources

        """


        # check if a directory exists for the chi plots. If not then make it.
        out_directory = self.fpath+'stat_plots/'
        if not os.path.isdir(out_directory):
            print("I am creating the river_plot/ directory to save your figures")
            os.makedirs(out_directory)

        
        # Set up fonts for plots
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['Liberation Sans'] # Liberation Sans is a free alternative to Arial. Albeit being quite universal, Arial is propietary. #PRAISE_FREE_AND_OPENSOURCE
        rcParams['font.size'] = label_size


        # Create a figure with required dimensions
        n_axis = 1
        fig = self.get_fig_right_size(size = size, n_axis =1 , facecolor = facecolor)

        gs = plt.GridSpec(100,100,bottom=0.15, left=0.10, right=0.95, top=0.95)
        ax1 = fig.add_subplot(gs[0:100,0:100], facecolor = "white")

        if(grid):
            ax1.grid(ls = 'dotted', lw = 0.1, c = "k", zorder = 5)


        # aggregating the data
        data_to_plot = []
        data_name = []
        n_data = []

        for bing in self.df_kp[binning].unique():
            data_to_plot.append(self.df_kp["delta_ksn"][self.df_kp[binning] == bing])
            data_name.append(str(bing))
            n_data.append(self.df_kp["delta_ksn"][self.df_kp[binning] == bing].shape[0])



        bp = ax1.boxplot(data_to_plot, labels = data_name, patch_artist = True)

        ## change outline color, fill color and linewidth of the boxes
        for box in bp['boxes']:
            # change outline color
            box.set( color='k', linewidth=1.5)
            # change fill color
            box.set( facecolor = '#848484' )

        ## change color and linewidth of the whiskers
        for whisker in bp['whiskers']:
            whisker.set(color='k', linewidth=1)

        ## change color and linewidth of the caps
        for cap in bp['caps']:
            cap.set(color='#2C2C2C', linewidth=1)

        ## change color and linewidth of the medians
        for median in bp['medians']:
            median.set(color='#E0E0E0', linewidth=1)

        ## change the style of fliers and their fill
        for flier in bp['fliers']:
            flier.set(marker='+', color='#A0A0A0', alpha=0.5)

        ax1.set_ylabel(r"$\Delta k_{sn}$")

        if(binning == "source_key"):
            xlabo = "Source keys"
        elif(binning == "basin_key"):
            xlabo = "Basin keys"
        else:
            xlabo = binning

        ax1.set_xlabel(xlabo)


        plt.savefig(out_directory + self.fprefix + "_kp_baw_%s.%s"%(binning,format), dpi = 500)
        plt.clf()

            





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
