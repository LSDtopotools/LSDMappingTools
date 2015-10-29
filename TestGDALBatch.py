# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 09:28:53 2015

@author: smudd
"""

#==============================================================================
# These are some scripts for testing the functionality of LSDMappingTools
#==============================================================================
# -*- coding: utf-8 -*-
"""
Created on 9 July 2015

@author: smudd
"""

import numpy as np
import LSDGDALBatchProcessing as GDALbatch

def TestGDALBatch():

    # Set up the test folders and parameters
    DataDirectory = "/home/smudd/SMMDataStore/analysis_for_papers/Chile_paper/SRTM_30metre/test"
    raster_format = "EHdr"
    target_format = "GTiff"
    merge_subfolder_name  = "test_merge"
    merge_filename = "t_merge"

    GDALbatch.GDALBatchMerge(DataDirectory,merge_subfolder_name,merge_filename,raster_format,target_format)
 
    
if __name__ == "__main__":
    TestGDALBatch()     