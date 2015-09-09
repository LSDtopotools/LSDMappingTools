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

# This function looks for all the fiules of a certain format in a directory and
# then translates them into a new format into a subdirectory named
# after the target format. 
def batch_convert(DataDirectory,raster_format,target_format):
 
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
        print "You have not selcted a valid raster format!"
        print "Options are ENVI, EHdr and GTiff"
        target_extension = "NULL"   

    # now make a directory   
    if target_extension != "NULL":
        
        target_directory = DataDirectory+target_format

        if not os.access(target_directory,os.F_OK):
            print "Making path: "
            os.mkdir(target_directory)
            print "I made a directory: " + target_directory
        else:
            print "Path: " +target_directory+" already exists."     

    # Now check the source format   
    if raster_format == "ENVI":
        raster_extension = ".bil"
    elif raster_format == "EHdr":
        raster_extension = ".bil"
    elif raster_format == "GTiff":
        raster_extension = ".tiff"
    else:
        print "You have not selcted a valid raster format!"
        print "Options are ENVI, EHdr and GTiff"
        raster_extension = "NULL"

    # find all the dataset of the source format  
    print "The data directory is: " + DataDirectory
    print "The raster extension is: " + raster_extension
    if raster_extension != "NULL":
        for FileName in glob(DataDirectory+"*"+raster_extension):
            print "found file: " + FileName
            subprocess.call(['gdalinfo',FileName])
        
        
        
if __name__ == "__main__":
    DataDirectory = "/home/smudd/SMMDataStore/analysis_for_papers/Chile_paper/SRTM_30metre"
    raster_format = "EHdr"
    batch_convert(DataDirectory,raster_format)           