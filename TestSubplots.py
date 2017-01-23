#==============================================================================
# Make a publication ready field sites figure of all the hillshades in a folder
# Will also read in any map images in EPS format and create subplots to show the
# location of the field sites.

# FJC 22/12/16
#==============================================================================
import matplotlib
matplotlib.use('Agg')
import LSDPlottingTools as LSDP

def make_field_sites_figure(DataDirectory):
    N_HSFiles = 4
    NRows = 3
    NCols = 2
    n_target_ticks = 6
    
    LSDP.field_sites(DataDirectory, N_HSFiles, NRows, NCols, n_target_ticks)
    
def make_flood_maps(DataDirectory):
    LSDP.multiple_flood_maps(DataDirectory)
    
def flood_maps_with_shapefile(DataDirectory):
    LSDP.flood_maps_with_shapefile(DataDirectory)

def flood_maps_shapefile_zoom(DataDirectory):
    LSDP.flood_maps_shapefile_with_zoom(DataDirectory)
    
if __name__ == "__main__":
    DataDirectory = "/home/s0923330/Datastore/Python/floodplain_initiation/results_comparison/field_data/"
    #make_field_sites_figure(DataDirectory)
    make_flood_maps(DataDirectory)
    #flood_maps_with_shapefile(DataDirectory)
    #flood_maps_shapefile_zoom(DataDirectory)
    

    