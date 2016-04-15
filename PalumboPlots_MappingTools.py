#==============================================================================
# These are some scripts for testing the functionality of LSDMappingTools
#==============================================================================
# -*- coding: utf-8 -*-
"""
Created on Tue May 05 14:08:16 2015

@author: smudd
"""

import numpy as np
import LSDPlottingTools as LSDP

def PalumboPlots():
    DataDirectory = "T://analysis_for_papers//Cosmo_paper//Palumbo_shield_plots//"
    Filename = "SpawnedBasin_07C13-(Q8).bil"
    Drapename = "SpawnedBasin_07C13-(Q8)_BASINS.bil"

    ThisFile = DataDirectory+Filename
    DrapeFile = DataDirectory+Drapename

    LSDP.BasicDrapedPlotGridPlot(ThisFile, DrapeFile, thiscmap='terrain',drape_cmap='gray',
                            colorbarlabel='Elevation in meters',clim_val = (2000,5000),
                            drape_alpha = 0.4,FigFileName = 'Basin.pdf',FigFormat = 'show')

def ShieldPlots_coarse():
    DataDirectory = "T://analysis_for_papers//Cosmo_paper//Palumbo_shield_plots//"
    Filename = "SpawnedBasin_07C13-coarse_SH.bil"

    ThisFile = DataDirectory+Filename

    LSDP.BasicDensityPlotGridPlot(ThisFile, thiscmap='summer',colorbarlabel='Topographic shielding coarse',
                             clim_val = (0.85,1),FigFileName = 'Shield_coarse.pdf', FigFormat = 'show')

def ShieldPlots_fine():
    DataDirectory = "T://analysis_for_papers//Cosmo_paper//Palumbo_shield_plots//"
    Filename = "SpawnedBasin_07C13-fine_SH.bil"

    ThisFile = DataDirectory+Filename

    LSDP.BasicDensityPlotGridPlot(ThisFile, thiscmap='summer',colorbarlabel='Topographic shielding fine',
                             clim_val = (0.85,1),FigFileName = 'Shield_fine.pdf', FigFormat = 'show')


if __name__ == "__main__":
    #fit_weibull_from_file(sys.argv[1]) 
    #TestNewMappingTools2() 
    #ResetErosionRaster()
    PalumboPlots()
    ShieldPlots_coarse()
    ShieldPlots_fine()
    