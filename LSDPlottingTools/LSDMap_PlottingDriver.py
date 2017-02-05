# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 10:20:18 2017

@author: smudd
"""

import LSDPlottingTools.LSDMap_ChiPlotting as LSDMap_CP
import LSDPlottingTools.LSDMap_BasicPlotting as LSDMap_BP
import LSDPlottingTools.LSDMap_OSystemTools as LSDOst
import os
import sys
import ast

class LSDMap_PlottingDriver(object):
    
    # The constructor: it needs a filename to read    
    def __init__(self,FileName):
        """This is the LSDMap_PlottingDriver object. 
        
        It manages plotting functions via a parameter file. 
        The parameter file syntax is similar to LSDTopoTools driver file syntax       
        
        Args:
            Filename (str): The name of the parameter file (with path and extension).
            
            Parameter files are made up of keywords followed by semicolons. 
            The keyword file_prefix *MUST* be included.
            
        Author: SMM
        """
        
        print("The filename is: "+FileName)         
        
        # This gets the path of the parameter file
        file_path = LSDOst.GetPath(FileName)
        self.FilePath = file_path        
        print("The path is: "+file_path)        


        #See if the parameter file exists
        if os.access(FileName,os.F_OK):
            this_file = open(FileName, 'r')
            lines = this_file.readlines()

            param_dict = {}
            for line in lines:            
                this_line =LSDOst.RemoveEscapeCharacters(line)
                                
                print("This line is: "+ this_line)
                if len(this_line) != 0:
                    if not this_line[0] == '#':
    
                        # get the parameters
                        split_line = this_line.split(':')
                    
                        if len(split_line) == 2:
                            # get rid of the spaces
                            remove_spaces = LSDOst.RemoveWhitespace(split_line[1])
                
                            # add this to the dict
                            param_dict[split_line[0]]=remove_spaces 
                        else:
                            print("This line, "+this_line+", isn't formatted properly, you need a colon after the parameter. Ignoring.")
                    else:
                        print("This line is a comment, the line is: ")
                        print(this_line)
                    
                    
            self.parameter_dict = param_dict
        else:
            print("Parameter file does not exit! I cannot proceed. Check your filenname.")
            sys.exit()

        # Now get out the file prefix
        self.FilePrefix = "None"
        if "file_prefix" in param_dict:
            this_prefix = param_dict["file_prefix"]          
            self.FilePrefix = LSDOst.RemoveWhitespace(this_prefix)
            
            print("The file prefix is: "+self.FilePrefix)
        else:
            print("Parameter file does does not contain the keyword 'file_prefix': Check your parameter file. I must exit now.")
            sys.exit()            
        
        self.base_faster_fname = str(self.FilePath+os.sep+self.FilePrefix+".bil")
        self.hs_fname = str(self.FilePath+os.sep+self.FilePrefix+"_hs.bil")
        self.basin_fname = str(self.FilePath+os.sep+self.FilePrefix+"_AllBasins.bil")
        self.chi_raster_fname = str(self.FilePath+os.sep+self.FilePrefix+"_chi_coord.bil")
        self.chi_csv_fname = str(self.FilePath+os.sep+self.FilePrefix+"_MChiSegmented.csv")
        self.basic_chi_csv_fname = str(self.FilePath+os.sep+self.FilePrefix+"_chi_coord_basins.csv")
        self.basin_csv_fname = str(self.FilePath+os.sep+self.FilePrefix+"_AllBasinsInfo.csv")
        
        
        # Now load and parse parameters
        self.create_plotting_switches()
        self.create_default_parameters()
        self.parse_parameter_dict()

    def create_plotting_switches(self):
        """This creates a dict containing plotting switches. It is later read to decide which plots to make
        
        Author: SMM
        """
        plotting_switches = {}
        
        # Basic plots
        plotting_switches["BasicDensityPlot"] = False
        plotting_switches["BasicDrapedPlotGridPlot"] = False
        plotting_switches["DrapedOverHillshade"] = False
        plotting_switches["DrapedOverFancyHillshade"] = False
        plotting_switches["BasinsOverFancyHillshade"] = False
        
        # Chi plots
        plotting_switches["BasicChiPlotGridPlot"] = False
        plotting_switches["BasicChiCoordinatePlot"] = False
        plotting_switches["BasicChannelPlotGridPlotCategories"] = False
        plotting_switches["ChiProfiles"] = False
        plotting_switches["StackedChiProfiles"] = False
        plotting_switches["StackedProfilesGradient"] = False 
        
        self.plotting_switches = plotting_switches


    def create_default_parameters(self):
        """This creates a dict containing plotting parameters. 
        It contains all the default parameters and is used to parse the parameter dict that is read from file. 
        
        Author: SMM
        """
        str_default_parameters = {}
        num_default_parameters = {}
        bool_default_parameters = {}
        
        
        # Add parameters for plotting here. 
        str_default_parameters["base_cmap"] = 'gray'
        str_default_parameters["cbar_label"] = 'Elevation in metres'
        str_default_parameters["FigFormat"] = "png"
        str_default_parameters["DrapeName"] = "None"
        str_default_parameters["drape_cmap"] = 'gray'        
        str_default_parameters["size_format"] = "esurf"        
        str_default_parameters["FigFileName"] = "None"  
        str_default_parameters["chan_net_csv"] = "None"
 

        num_default_parameters["clim_val"] = (0,0)
        num_default_parameters["drape_alpha"] = 0.6
        num_default_parameters["source_chi_threshold"] = 10
        num_default_parameters["grouped_basin_list"] = []
        num_default_parameters["basin_order_list"] = []
        num_default_parameters["spread"] = 10
        num_default_parameters["basin_rename_list"] = []
        num_default_parameters["elevation_threshold"] = 0
        num_default_parameters["source_thinning_threshold"]= 0
        
        bool_default_parameters["label_sources"] = False
        bool_default_parameters["is_log"] = False
        bool_default_parameters["plot_M_chi"]= False
        bool_default_parameters["plot_segments"] = False

        
        
        self.num_default_parameters = num_default_parameters
        self.str_default_parameters = str_default_parameters
        self.bool_default_parameters = bool_default_parameters        
        
    def parse_parameter_dict(self):
        """This parses the parameter dict. 
        
        Its purpose is to read both settings and switches for what sort of plots to be made.
        
        Author: SMM
        """
        
        # First we find the switches        
        for key in self.plotting_switches:
            if key in self.parameter_dict:
                print("I found the switch: "+str(key)+ " in the parameter file, the value is:"+self.parameter_dict[key])
                
                # check to see if the string is true
                if self.parameter_dict[key] in ["True","true","t","T"]: 
                    self.plotting_switches[key] = True
                
                
        # now the boolean parameters
        for key in self.bool_default_parameters:
            if key in self.parameter_dict:
                print("I found the parameter: "+str(key)+ " in the parameter file.")
                
                if self.parameter_dict[key] in ["True","true","t","T"]:
                    self.bool_default_parameters[key] = True
                else:
                    self.bool_default_parameters[key] = False

        # now the parameters that are numbers (this can include lists)
        for key in self.num_default_parameters:
            if key in self.parameter_dict:
                print("I found the parameter: "+str(key)+ " in the parameter file.")
                self.num_default_parameters[key]  = ast.literal_eval(self.parameter_dict[key])               
                print("The value is: "+str(self.num_default_parameters[key]))

        # now the parameters that are strings
        for key in self.str_default_parameters:
            if key in self.parameter_dict:
                print("I found the parameter: "+str(key)+ " in the parameter file.")
                self.str_default_parameters[key]  = self.parameter_dict[key] 
                print("The string is: "+self.str_default_parameters[key])
                
        # now concatenate the dicts
        self.plotting_parameters = {}
        self.plotting_parameters.update(self.bool_default_parameters)
        self.plotting_parameters.update(self.str_default_parameters)
        self.plotting_parameters.update(self.num_default_parameters)
        
        # reset the defaults, just to be safe
        self.create_default_parameters()
                
    def plot_data(self):
        """This is the bit that actually plots the data.
        
        """
        import LSDPlottingTools.LSDMap_PointTools as LSDMap_PD

            
        if self.plotting_switches["BasicDensityPlot"]: 
            
            # Check to see if there is a filename. If not set a default file name
            if self.plotting_parameters["FigFileName"] == "None":
                self.plotting_parameters["FigFileName"] = self.FilePath+os.sep+self.FilePrefix+"BDP."+self.plotting_parameters["FigFormat"]     
            
            print("Hey there partner, I am making a grid plot.")
            LSDMap_BP.BasicDensityPlot(self.base_faster_fname, 
                                       self.plotting_parameters["base_cmap"],
                                       self.plotting_parameters["cbar_label"],
                                       self.plotting_parameters["clim_val"],
                                       self.plotting_parameters["FigFileName"],
                                       self.plotting_parameters["FigFormat"],
                                       self.plotting_parameters["size_format"],
                                       self.plotting_parameters["is_log"])
                
            
        if self.plotting_switches["BasicDrapedPlotGridPlot"]:
            # Check to see if there is a filename. If not set a default file name
            if self.plotting_parameters["FigFileName"] == "None":
                self.plotting_parameters["FigFileName"] = self.FilePath+os.sep+self.FilePrefix+"BDPDPG."+self.plotting_parameters["FigFormat"]     
 
            LSDMap_BP.BasicDrapedPlotGridPlot(self.base_faster_fname, 
                                       self.plotting_parameters["DrapeName"],
                                       self.plotting_parameters["base_cmap"],
                                       self.plotting_parameters["drape_cmap"],                                
                                       self.plotting_parameters["cbar_label"],
                                       self.plotting_parameters["clim_val"],
                                       self.plotting_parameters["drape_alpha"],
                                       self.plotting_parameters["FigFileName"],
                                       self.plotting_parameters["FigFormat"])
                                              
        if self.plotting_switches["BasinsOverFancyHillshade"]:
            # Check to see if there is a filename. If not set a default file name
            if self.plotting_parameters["FigFileName"] == "None":
                self.plotting_parameters["FigFileName"] = self.FilePath+os.sep+self.FilePrefix+"Basins."+self.plotting_parameters["FigFormat"]     
            
            print("Type of base raster is: "+str(type(self.base_faster_fname)))
            print("The chan net csv is: " +str(self.plotting_parameters["chan_net_csv"]))
            print("The drape cmap is: "+self.plotting_parameters["drape_cmap"])
            
            thisPointData = LSDMap_PD.LSDMap_PointData(self.basin_csv_fname)
            
            if self.plotting_parameters["chan_net_csv"] == "None":
                chanPointData = "None"
            else:
                chanPointData = LSDMap_PD.LSDMap_PointData(self.plotting_parameters["chan_net_csv"])
            
            
            LSDMap_BP.BasinsOverFancyHillshade(self.base_faster_fname, 
                                       self.hs_fname,
                                       self.basin_fname,
                                       self.basin_csv_fname,
                                       thisPointData,                                      
                                       self.plotting_parameters["base_cmap"],
                                       self.plotting_parameters["drape_cmap"],                                
                                       self.plotting_parameters["clim_val"],
                                       self.plotting_parameters["drape_alpha"],
                                       self.plotting_parameters["FigFileName"],
                                       self.plotting_parameters["FigFormat"],
                                       self.plotting_parameters["elevation_threshold"],
                                       self.plotting_parameters["grouped_basin_list"], 
                                       self.plotting_parameters["basin_rename_list"],
                                       self.plotting_parameters["spread"],
                                       chanPointData,
                                       self.plotting_parameters["label_sources"],
                                       self.plotting_parameters["source_chi_threshold"],
                                       self.plotting_parameters["size_format"])

        if self.plotting_switches["BasicChiCoordinatePlot"]:
            print("I am plotting a basic chi plot!")
            
            # Check to see if there is a filename. If not set a default file name
            if self.plotting_parameters["FigFileName"] == "None":
                self.plotting_parameters["FigFileName"] = self.FilePath+os.sep+self.FilePrefix+"Chi."+self.plotting_parameters["FigFormat"]     

                
                
            thisBasinData = LSDMap_PD.LSDMap_PointData(self.basin_csv_fname)
            chi_drape_cname = 'CMRmap_r'
            #chi_drape_cname = 'brg_r'
            cbar_lablel = "$\chi$ (m)"

            LSDMap_CP.BasicChiCoordinatePlot(self.hs_fname, 
                                       self.chi_raster_fname,
                                       self.basic_chi_csv_fname,                                      
                                       self.plotting_parameters["base_cmap"],
                                       chi_drape_cname,
                                       cbar_lablel,
                                       self.plotting_parameters["clim_val"],
                                       self.plotting_parameters["basin_order_list"],
                                       thisBasinData,
                                       self.basin_fname,
                                       self.plotting_parameters["drape_alpha"],
                                       self.plotting_parameters["FigFileName"],
                                       self.plotting_parameters["FigFormat"],
                                       self.plotting_parameters["size_format"])            
            
        if self.plotting_switches["ChiProfiles"]:
            print("I am plotting a basic chi profile plot!")
            
            # Check to see if there is a filename. If not set a default file name
            if self.plotting_parameters["FigFileName"] == "None":
                self.plotting_parameters["FigFileName"] = self.FilePath+os.sep+self.FilePrefix+"ChiProfile."+self.plotting_parameters["FigFormat"]     

                
                
            thisBasinData = LSDMap_PD.LSDMap_PointData(self.basin_csv_fname)
            chi_drape_cname = 'CMRmap_r'
            #chi_drape_cname = 'brg_r'
            cbar_lablel = "$\chi$ (m)"
            
            print("The csv filename is: "+ self.chi_csv_fname)
            
            LSDMap_CP.ChiProfiles(self.chi_csv_fname,
                                       self.plotting_parameters["FigFileName"],
                                       self.plotting_parameters["FigFormat"],
                                       self.plotting_parameters["basin_order_list"],      
                                       self.plotting_parameters["basin_rename_list"],      
                                       self.plotting_parameters["label_sources"], 
                                       self.plotting_parameters["elevation_threshold"], 
                                       self.plotting_parameters["source_thinning_threshold"], 
                                       self.plotting_parameters["plot_M_chi"],                                       
                                       self.plotting_parameters["plot_segments"],
                                       self.plotting_parameters["size_format"])            
            