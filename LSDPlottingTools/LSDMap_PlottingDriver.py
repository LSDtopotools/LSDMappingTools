# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 10:20:18 2017

@author: smudd
"""

from . import LSDMap_BasicPlotting as LSDMap_BP
from . import LSDMap_OSystemTools as LSDOst
import os
import sys


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

                # get the parameters
                split_line = this_line.split(':')
                
                # get rid of the spaces
                remove_spaces = LSDOst.RemoveWhitespace(split_line[1])
            
                # add this to the dict
                param_dict[split_line[0]]=remove_spaces 

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
        
        self.base_faster_fname = self.FilePath+os.sep+self.FilePrefix+".bil"
        self.hs_fname = self.FilePath+os.sep+self.FilePrefix+"_hs.bil"
        self.chi_csv_fname = self.FilePath+os.sep+self.FilePrefix+"_MChi_segmented.csv"
        self.basin_csv_fname = self.FilePath+os.sep+self.FilePrefix+"_AllBasinsInfo.csv"
        
        
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
        plotting_switches["BasicDensityPlotGridPlot"] = False
        plotting_switches["BasicDrapedPlotGridPlot"] = False
        plotting_switches["DrapedOverHillshade"] = False
        plotting_switches["DrapedOverFancyHillshade"] = False
        plotting_switches["BasinsOverFancyHillshade"] = False
        
        # Chi plots
        plotting_switches["BasicChiPlotGridPlot"] = False
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
        default_parameters = {}
        
        # Add parameters for plotting here. 
        default_parameters["base_cmap"] = 'gray'
        default_parameters["cbar_label"] = 'Elevation in metres'
        default_parameters["clim_val"] = (0,0)
        default_parameters["FigFormat"] = "png"
        default_parameters["drape_alpha"] = 0.6
        default_parameters["size_format"] = "esurf"
        default_parameters["source_chi_threshold"] = 10
        default_parameters["grouped_basin_list"] = []
        default_parameters["spread"] = 10
        default_parameters["label_sources"] = False
        default_parameters["FigFileName"] = "None"
        
        
        self.default_parameters = default_parameters
        
        
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
                
                
        # now the parameters
        for key in self.default_parameters:
            if key in self.parameter_dict:
                print("I found the parameter: "+str(key)+ " in the parameter file.")
                self.default_parameters[key] = self.parameter_dict[key] 
                
    def plot_data(self):
        """This is the bit that actually plots the data.
        
        """
        
        if self.plotting_switches["BasicDensityPlot"]:
            LSDMap_BP.BasicDensityPlot(self.base_faster_fname, 
                                       self.default_parameters["base_cmap"],
                                       self.default_parameters["cbar_label"],
                                       self.default_parameters["clim_val"])
            
        if self.plotting_switches["BasicDensityPlotGridPlot"]: 
            
            # Check to see if there is a filename. If not set a default file name
            if self.default_parameters["FigFileName"] == "None":
                self.default_parameters["FigFileName"] = self.FilePath+os.sep+self.FilePrefix+"BDPG."+self.default_parameters["FigFormat"]     
            
            print("Hey there partner, I am making a grid plot.")
            LSDMap_BP.BasicDensityPlotGridPlot(self.base_faster_fname, 
                                       self.default_parameters["base_cmap"],
                                       self.default_parameters["cbar_label"],
                                       self.default_parameters["clim_val"],
                                       self.default_parameters["FigFileName"],
                                       self.default_parameters["FigFormat"])
                
                