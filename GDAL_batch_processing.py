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

def batch_convert(DataDirectory,raster_format):
 
    NewDataDirectory = LSDost.ReformatSeperators(DataDirectory)   
    DataDirectory = LSDost.AppendSepToDirectoryPath(NewDataDirectory)
   
     
   
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