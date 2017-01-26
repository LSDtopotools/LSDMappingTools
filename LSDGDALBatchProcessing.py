# -*- coding: utf-8 -*-
"""
Created on Wed Sep 09 16:52:27 2015

@author: smudd
"""

import numpy as np
from glob import glob
import LSDOSystemTools as LSDost
import os
import shutil
import subprocess

# This function looks for all the files of a certain format in a directory and
# then translates them into a new format into a subdirectory named
# after the target format. 
def GDALBatchConvert(DataDirectory,raster_format,target_format):
 
    NewDataDirectory = LSDost.ReformatSeperators(DataDirectory)   
    DataDirectory = LSDost.AppendSepToDirectoryPath(NewDataDirectory)

    # Check the target format   
    if target_format == "ENVI":
        target_extension = ".bil"
    elif target_format == "EHdr":
        target_extension = ".bil"
    elif target_format == "GTiff":
        target_extension = ".tiff"
    else:
        print("You have not selcted a valid raster format!")
        print("Options are ENVI, EHdr and GTiff")
        target_extension = "NULL"   

    # now make a directory   
    if target_extension != "NULL":
        
        target_directory = DataDirectory+target_format

        if not os.access(target_directory,os.F_OK):
            print("Making path: ")
            os.mkdir(target_directory)
            print("I made a directory: " + target_directory)
        else:
            print("Path: " +target_directory+" already exists.")     

    # Now check the source format   
    if raster_format == "ENVI":
        raster_extension = ".bil"
    elif raster_format == "EHdr":
        raster_extension = ".bil"
    elif raster_format == "GTiff":
        raster_extension = ".tif"
    else:
        print("You have not selcted a valid raster format!")
        print("Options are ENVI, EHdr and GTiff")
        raster_extension = "NULL"

    # find all the dataset of the source format  
    print("The data directory is: " + DataDirectory)
    print("The raster extension is: " + raster_extension)
    if raster_extension != "NULL":
        for FileName in glob(DataDirectory+"*"+raster_extension):
            print("found file: " + FileName)
            subprocess.call(['gdalinfo',FileName])
        
def GDALBatchMerge(DataDirectory,merge_subfolder_name,merge_filename,raster_format,target_format):
    
    NewDataDirectory = LSDost.ReformatSeperators(DataDirectory)   
    DataDirectory = LSDost.AppendSepToDirectoryPath(NewDataDirectory)

    # get the name of the data directory into which the file should be merged
    merge_DataDirectory = DataDirectory+merge_subfolder_name
    mDataDriectory = LSDost.AppendSepToDirectoryPath(merge_DataDirectory)

    # make the directory
    if not os.access(mDataDriectory,os.F_OK):
        print("Making path: ")
        os.mkdir(mDataDriectory)
        print("I made a directory: " + mDataDriectory)
    else:
        print("Path: " +mDataDriectory+" already exists.")             

    # Check the source format   
    if raster_format == "ENVI":
        raster_extension = ".bil"
    elif raster_format == "EHdr":
        raster_extension = ".bil"
    elif raster_format == "GTiff":
        raster_extension = ".tif"
    else:
        print("You have not selcted a valid raster format!")
        print("Options are ENVI, EHdr and GTiff")
        raster_extension = "NULL"

    # Check the target format. Default is geotiff  
    if target_format == "ENVI":
        target_extension = ".bil"
    elif target_format == "EHdr":
        target_extension = ".bil"
    elif target_format == "GTiff":
        target_extension = ".tif"
    else:
        print("You have not selcted a valid raster format!")
        print("Defaulting to GTiff")
        target_format == "GTiff"        
        target_extension = ".tif"
    
    # set the name of the target file
    target_FileName = mDataDriectory+merge_filename+target_extension

    # find all the dataset of the source format  
    print("The data directory is: " + DataDirectory)
    print("The raster extension is: " + raster_extension)
    if raster_extension != "NULL":
        
        # Set up the list for holding command prompt commands        
        command_prompt = []
        command_prompt.append("gdal_merge.py")
        command_prompt.append("-of")
        command_prompt.append(target_format)
        command_prompt.append("-o")
        command_prompt.append(target_FileName)
                
        
        for FileName in glob(DataDirectory+"*"+raster_extension):
            print("found file: " + FileName)
            command_prompt.append(FileName)            
            
        print("The subprocess call is: ")
        print(command_prompt)         
        subprocess.call(command_prompt)
