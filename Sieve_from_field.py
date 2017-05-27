#Importing the stuff we need
import pandas as bamboo_bears
import numpy as np
from matplotlib import pyplot as plt
from LSDMapFigure.PlottingRaster import MapFigure
from LSDMapFigure.PlottingRaster import BaseRaster
from LSDPlottingTools import colours as lsdcolours
from LSDPlottingTools import init_plotting_DV
import LSDPlottingTools as LSDP


# the cumulative grainsize curve and calculate the D50/84
def grainsize_curve_from_sieving(name_csv, directory):
    # Importing the data
    df = bamboo_bears.read_csv(directory+name_csv, sep =",") #reading the file
    data = np.array(df.values) # Converting to numpy array
    data_out = data[:4]

    # data_processing
    power_third = float(1./3.)

    #calculation of the missing Diameters
    for i in range(data.shape[0]):
        if(data[i,0] == -1):
            data[i,0] = 1000*np.power((6*data[i,1]/(2650*np.pi)),(power_third)) # D = [6M/2650]1/3
            if(data[i,0] < 80):
                data_out[3,1] += data[i,1]
            else:
                data_out = np.insert(data_out,data_out.shape[0],data[i],axis = 0)

    # Calculation of the cumulative curve
    ## sorting
    data_out = data_out[data_out[:,0].argsort()]
    ## cumulation
    data_cumul = data_out
    for i in range(1,data_out.shape[0]):
        data_cumul[i,1] += data_cumul[i-1,1]
    ## percentages
    maximum = data_cumul[:,1].max()
    data_cumul[:,1] = data_cumul[:,1]*100/maximum
    data_cumul = np.insert(data_cumul,0,[0,0],axis = 0)

    # determination of D50
    ## determining the closest values to D50
    for i in range(data_cumul.shape[0]):
        if(data_cumul[i,1]>=50):
            index_50_plus = i
            index_50_minus = i-1
            break
    ## getting the line equation
    m_fifty = (data_cumul[index_50_plus,1] - data_cumul[index_50_minus,1])/(data_cumul[index_50_plus,0] - data_cumul[index_50_minus,0])
    b_fifty = data_cumul[index_50_plus,1] - m_fifty*data_cumul[index_50_plus,0]
    d_fifty = (50 - b_fifty)/m_fifty

    # determination of D84
    ## determining the closest values to D84
    for i in range(data_cumul.shape[0]):
        if(data_cumul[i,1]>=84):
            index_84_plus = i
            index_84_minus = i-1
            break
    ## getting the line equation
    m_eighty_four = (data_cumul[index_84_plus,1] - data_cumul[index_84_minus,1])/(data_cumul[index_84_plus,0] - data_cumul[index_84_minus,0])
    b_eighty_four = data_cumul[index_84_plus,1] - m_eighty_four*data_cumul[index_84_plus,0]

    d_eighty_four = (84 - b_eighty_four)/m_eighty_four
    # Plotting the result
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot(data_cumul[:,0],data_cumul[:,1],lw =1.2,  color = "black", marker = "o", markersize = 4, aa = True)
    ax.set_xlabel("Diameters in mm")
    ax.set_xlim(0,data_cumul[:,0].max())
    ax.set_ylabel("Percentage Mass finer")
    ax.set_ylim(0,100)
    ax.annotate('D50: ' + str(d_fifty), xy=(data_cumul[:,0].max() - 20,90))
    ax.annotate('D84: ' + str(d_eighty_four), xy=(data_cumul[:,0].max() - 20,85))
    plt.savefig(directory+name_csv[:-3]+"png",dpi = 300)
    plt.clf()
    return d_fifty, d_eighty_four
#
# Compile and sort the data and create a synthesis figure per river
def compile_several_d(tab_def,names, Directory):

###### sorting by river
    tempTable_number =[]
    tempTable_D50 = []
    tempTable_D84 =[]
    for i in range(len(tab_def)):
        if i == 0:
            tempTable_number.append(tab_def[i][0])
            tempTable_D50.append(tab_def[i][1])
            tempTable_D84.append(tab_def[i][2])
        else:
            if(tab_def[i][3] == tab_def[i-1][3] and i != len(tab_def)-1):
                tempTable_number.append(tab_def[i][0])
                tempTable_D50.append(tab_def[i][1])
                tempTable_D84.append(tab_def[i][2])
            else:
                ####plotting when the River change
                if(i == len(tab_def)-1): # to add the last data
                    tempTable_number.append(tab_def[i][0])
                    tempTable_D50.append(tab_def[i][1])
                    tempTable_D84.append(tab_def[i][2])
                #conveting to numpy because it is easier to use
                ploti0 = np.array(tempTable_number)
                ploti1 = np.array(tempTable_D50)
                ploti2 = np.array(tempTable_D84)
                fig, ax = plt.subplots(figsize=(12, 8))
                ax.plot(ploti0,ploti1,lw =1.2,  color = "black", marker = "o", markersize = 4, aa = True, label = "D50")
                ax.plot(ploti0,ploti2,lw =1.2,  color = "red", marker = "s", markersize = 4, aa = True, label = "D84")
                ax.set_xlabel("Sieves")
                ax.set_title(tab_def[i-1][3])
                ax.set_ylabel("Percentage Mass finer")
                ax.set_ylim(0,100)
                plt.savefig(Directory+"Sieving_compiled_"+tab_def[i-1][3]+".png",dpi = 300)
                plt.clf()
                # going to the next river
                tempTable_number =[]
                tempTable_D50 = []
                tempTable_D84 =[]
                tempTable_number.append(tab_def[i][0])
                tempTable_D50.append(tab_def[i][1])
                tempTable_D84.append(tab_def[i][2])



    plt.clf()
# Function to load parametersfrom the general file
def load_general_sieving_file(Directory, fname):
    df = bamboo_bears.read_csv(Directory + fname, sep =",")
    return df

# Create a new csv file with the D50/84 added in case you want it for later
def print_new_csv_with_D(Directory,dataf,dfi,deig):
    dataf["D50"] = bamboo_bears.Series(dfi, index = dataf.index)
    dataf["D84"] = bamboo_bears.Series(deig, index = dataf.index)
    indexator = 0
    tempo = []
    for i in range(dataf.river.values.shape[0]):
        if (i == 0):
            tempo.append(indexator)
        else:
            if(dataf.river.values[i] == dataf.river.values[i-1]):
                tempo.append(indexator)
            else:
                indexator +=1
                tempo.append(indexator)

    dataf["river_ID"] = bamboo_bears.Series(tempo, index = dataf.index)
    dataf.to_csv(Directory+"sieving_processed.csv", index = False)
    return "sieving_processed.csv"

if __name__ == '__main__':
    print("I will process your data from field sieving")
################ Parameters to changes ###########
    Directory = "C:/Data/Field_2017/Grainsize/Sieveing/" # Path to your sieving files
    gname = "sieving.csv" # Name of the general sieving file
    maps = True # Do you want to create a map with the grainsize plotted?
    D_map = "C:/Data/ORDI_BORIS_ROUMANIE/Romania/" # path to the main raster
    N_map = "Romania_study2.bil" #name of the raster
    HS_name = "Romania_study_hs2.bil" # ' name of the hs raster'
    river_file = "Knickpoint_MChiSegmented.csv" # name of the river network csv file, any file that contain at least the latitude longitude of the river points

##################### Do not change the following if you don't know why
    file_info = load_general_sieving_file(Directory, gname)
    tab_d = []
    pdf = []
    pde = []
    file_info.csv_file_name.values
    print("Let me load and preprocess your "+str(len(file_info.csv_file_name.values))+ " files")
    for i in range(len(file_info.csv_file_name.values)):
        print("File "+  str(i+1))
        tdf,tde = grainsize_curve_from_sieving(file_info.csv_file_name.values[i],Directory)
        tab_d.append([i+1,tdf,tde,file_info.river.values[i]])
        pdf.append(tdf)
        pde.append(tde)
        print("Your D50 is: %s and your D84 is %s" %(tdf,tde))
        print("Done plus figure created in "+Directory+file_info.csv_file_name.values[i][:-3]+"png")

    print("I am now compiling the data and creatng the D50/84 evolution figure")
    compile_several_d(tab_d,file_info.csv_file_name.values,Directory)
    print("Your compiled filed is named "+ Directory+ "Sieving_compiled.png")
################ printing the maps
    if(maps):
        print("I will now create a map with the sieving data sized by D50")
        file_info = print_new_csv_with_D(Directory,file_info,pdf,pde)

        thisPointData = LSDP.LSDMap_PointData(Directory + file_info) # Load the point file #1, add a similar line with different name if you have more than one point file.
        riverPointData = LSDP.LSDMap_PointData(D_map + river_file)
        plt.clf() # Ignore this line

        MF = MapFigure(N_map, D_map,coord_type="UTM_km") # load the background raster
        MF.add_drape_image(HS_name,D_map, # Calling the function will add a drapped raster on the top of the background one
                            colourmap = "gray", # colormap used for this raster, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                            alpha=0.5, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                            show_colourbar = False, # Well, this one is explicit I think
                            colorbarlabel = "Colourbar") # Name of your Colourbar, it might bug though
        MF.add_point_data( riverPointData, # this function plot the requested point file using the lat/long column in the csv file
                           column_for_plotting = "None",  # Column used to color the data
                           this_colourmap = "RdBu_r", # Colormap used, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                           colorbarlabel = "River", # Label
                           scale_points = False, # All the point will have the same size if False
                           column_for_scaling = "None", # If scale point True, you can scale the size of your points using one of the columns
                           scaled_data_in_log = False, # If scale point True, you can log the scaling
                           max_point_size = 12, # max size if scale point True again
                           min_point_size = 1, # You should be able to guess that one now
                           coulor_log = False, # do you want a log scale for your colorbar ?
                           coulor_manual_scale = [], #Do you want to manually limit the scale of your colorbar? if not let is false
                           manual_size = 0.25, # If none of above is choosen but you want to put another value than 0.5 to scale your point
                           alpha = 0.2, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                           minimum_log_scale_cut_off = -10) # you probably won't need this
        MF.add_point_data( thisPointData, # this function plot the requested point file using the lat/long column in the csv file
                           column_for_plotting = "river_ID",  # Column used to color the data
                           this_colourmap = "RdBu_r", # Colormap used, see http://matplotlib.org/users/colormaps.html for examples, put _r at the end of a colormap to get the reversed version
                           colorbarlabel = "River", # Label
                           scale_points = True, # All the point will have the same size if False
                           column_for_scaling = "D50", # If scale point True, you can scale the size of your points using one of the columns
                           scaled_data_in_log = False, # If scale point True, you can log the scaling
                           max_point_size = 70, # max size if scale point True again
                           min_point_size = 15, # You should be able to guess that one now
                           coulor_log = False, # do you want a log scale for your colorbar ?
                           coulor_manual_scale = [], #Do you want to manually limit the scale of your colorbar? if not let is false
                           manual_size = 0.5, # If none of above is choosen but you want to put another value than 0.5 to scale your point
                           alpha = 1, # transparency of this specific layer, 0 for fully transparent (why not) and 1 for fully opaque
                           minimum_log_scale_cut_off = -10) # you probably won't need this
        ax_style = "Normal" # Ignore this
        MF.save_fig(fig_width_inches = 10, FigFileName = Directory + "Map_Sieve.png", axis_style = ax_style, Fig_dpi = 500) # Save the figure
