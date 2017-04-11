import matplotlib
matplotlib.use('Agg')
import os
import subprocess
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import pandas
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import matplotlib

# Parameters

## Axis
Chi_plot = False # If true, x axis will be Chi instead of flow distance

##### Ignore this, this is lazyness ####
########################################
if(Chi_plot):                          #
    Chi_plot_s = "Chi"                 #
else:                                  #
    Chi_plot_s = "Flow distance"       #
########################################

grid_plot = False # plot a grid

x_label_txt = "%s (m)" %(Chi_plot_s) # ignore this
y_label_txt = "Elevation (m)" # Y axis label
fig_title = "" # Figure title in case you want to display it


## Other options
ext = (".png") # ignore this
pdpi = 500 # Resolution of the figure
median_leg = False # Will adjust the max_mChi legend in regards to the median
median_mult = 6 # multiplicator applied to the median to fix the max m_Chi
reductor = 1 # reduction of the mChi legend
underdrapping = 75 # thinckness of the underlying geology color plot
space_between_river = 50 # space between the river plot and the geology
alphaGeol = 75 # transparency of the geology in percent I guess???
Chi_sort = False # if true, use the following threshold to exclude the outlier Chi values
Chi_threshold = 2500 # will erase the datapoints over this threshold in case you have really weird data
horizontal_geology = False # will make the geology as horizontal bar with vertical line rather than underlying the main river
y_horizontal_value = 50 # position of this bar in case of horizontal geology
size_point = 2 # Size of the points
automated_legend = True # will automatically plot the legend from the geologic_keys csv file
separated_legend = True # will print the legend in a separated file
basin_key = 0 # Your basin key

if(automated_legend):
    rDirectory3 = "/home/s1675537/PhD/DataStoreBoris/GIS/Data/Carpathian/Geology_tes/" # geologic key directory path
    geologic_keys = "Carpathians_WGS84_new_lithokey.csv" # file that host the geological keys
    rDirectory4 = "/home/s1675537/PhD/DataStoreBoris/GIS/Data/Carpathian/Geology_tes/" # geologic key signification directory path
    glimolofile = "glimology.csv" # file that host the geological keys signification

####### your file
rDirectory1 = "/home/s1675537/PhD/DataStoreBoris/GIS/Data/Carpathian/Basins/small_Olt/"
baseName1 = "smOlt"
fname1 = baseName1+"_geol.csv"

fdist_col1 = 7 # column that host the flow distance info
elev_col1 = 6 # column that host the elevation
mchi_col1 = 9 # column that host the mchi
geol_col1 = 8 # column that host the geology key (it has to be number)
sources_col1 = 11 # column that host the source key (it has to be number)
Chi_col1 = 4 # column that host the Chi value
basin_key_col1 = 3 # column that host the basin keys

## Writing directory and name
wDirectory = "/home/s1675537/PhD/DataStoreBoris/presentation/Poster_EGU_2017/Poster/Figure/Geologic_Figures/"
wname = "TESTEMMA_"+baseName1




###############################
######## Computing part #######
###############################

if(Chi_plot):
    fdist_col1 = Chi_col1
    print("The X axis will be the Chi value")
else:
    print("The X axis will be the Flow distance")

print("Checking the output path")
# Check if the writing path exists and creates it if necessary
if (os.path.isfile(wDirectory) == False):
    print("Ignore the following error message (if any)")
    system_call = "mkdir " + wDirectory
    p = subprocess.call(system_call, shell=True)


if(os.path.isfile(rDirectory1+fname1)):


    print("I am now loading the data")

    # Most efficient way to load csv file so far, but require pandas
    df = pandas.read_csv(rDirectory1 + fname1, sep=",")
    X = np.array(df.values)
    max1 = X[:,fdist_col1].max()

    ##### Sort the data from a basin
    X = X[X[:,basin_key_col1] == basin_key]

    if(automated_legend):
        df = pandas.read_csv(rDirectory3 + geologic_keys, sep=",")
        GeoKey = df.values
        df = pandas.read_csv(rDirectory4 + glimolofile, sep=",")
        glimolokey = df.values


    print("I am now processing the data")

    ## I am processing the Xaxis for the right side: I am still thinking of how to do it for Chi

    if(Chi_plot and Chi_sort):
        print("I am excluding Chi values over your threshold")
        X = X[X[:,fdist_col1] < Chi_threshold]

    ##### Identification of the main channel for each plot to underdrap (is it a word?? if not it should be) the geology ###
    print("I am selecting the main channel on each side and I am underdrapping the corresponding geology")
    unique,pos = np.unique(X[:,sources_col1],return_inverse=True)
    counts = np.bincount(pos)
    maxpos = counts.argmax()
    main_source1 = maxpos

    geo1 = X[X[:,sources_col1] == main_source1]



    cGeology = geo1[:,geol_col1]
    yGeology = geo1[:,elev_col1]

    if(horizontal_geology):
        oldYGeology = yGeology
        yGeology = yGeology/yGeology + y_horizontal_value

    else:
        yGeology = yGeology - underdrapping
    xGeology = geo1[:,fdist_col1]

    print("The sources ID of these is %s " %(main_source1))


    ##### Elevation against flow distance/Chi colored by mChi #####
    # appending the two sides
    flow_distance = X[:,fdist_col1]
    elevation = X[:,elev_col1]
    mchi = X[:,mchi_col1]

    # Positionning, If you wanna change it
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_position([0.08, 0.125, 0.85, 0.80])

    # Plotting the geology (unexpectidly painful)
    lastGeol = -666
    thisGeol = cGeology[0]
    tempPlotx = np.array([])
    tempPloty = np.array([])
    tempPlotc = np.array([])
    GeoKey_used = []

    cNorm  = matplotlib.colors.Normalize(vmin=GeoKey[:,1].min(), vmax=GeoKey[:,1].max())
    scalarMap = cm.ScalarMappable(norm = cNorm, cmap= plt.cm.Accent)
    max_y = elevation.max()+elevation.max()/20
    min_y = 0+elevation.max()/20
    min_x = 0-flow_distance.max()/20
    max_x = flow_distance.max()+flow_distance.max()/20
    interval = max_y + min_y
    for i in range(xGeology.shape[0]):

        thisGeol = cGeology[i]
        if(thisGeol == lastGeol and i<(xGeology.shape[0]-1)):
            tempPlotx = np.append(tempPlotx, xGeology[i])
            tempPloty = np.append(tempPloty, yGeology[i])
            tempPlotc = np.append(tempPlotc, scalarMap.to_rgba(cGeology[i]))

        else:
            if(tempPlotc.shape[0]>4): # (the rgb file has a minium size to be viable)

                ax.fill_between(tempPlotx,tempPloty - space_between_river,tempPloty + underdrapping - space_between_river, color = tempPlotc, alpha = alphaGeol, lw = 0)
            if(automated_legend):
                temp_bool = True
                for j in GeoKey_used:
                    if j == cGeology[i]:
                        temp_bool = False

                if(temp_bool):
                    GeoKey_used.append(int(cGeology[i]))

                tempPlotx = np.array([])
                tempPloty = np.array([])
                tempPlotc = np.array([])
                if(horizontal_geology):

                    ax.axvline(x=xGeology[i], ymin = (y_horizontal_value)/interval, ymax = (min_y+oldYGeology[i])/interval , linewidth=0.5, color='k', ls="-")
                    print("Geological contact at %s / %s" )%(xGeology[i],oldYGeology[i])




        lastGeol = thisGeol

    #### Dealing with the legend ####
    if(automated_legend):
        print("I am now creating the legend")
        geolegend = []
        geoNumLeg = []
        # Selected the used geological key, avoiding the unused
        for i in range(GeoKey.shape[0]):
            for j in GeoKey_used:
                if(j == GeoKey[i][1]):
                    for h in range(glimolokey.shape[0]):

                        if(GeoKey[i][0] == glimolokey[h][0]):

                            geolegend.append(glimolokey[h])
                            geoNumLeg.append(GeoKey[i][1])

        #plotting parameters
        rx = (min_x + (max_x-min_x)*0.05)
        ry = (max_y - (max_y-min_y)*0.1)
        rw = (min_x + (max_x-min_x)*0.10)
        rh = ((max_y-min_y)*0.05)
        # plotting
        if(separated_legend == False):
            for i in range(len(geolegend)):
                ax.add_patch(patches.Rectangle((rx, ry),rw, rh,facecolor = scalarMap.to_rgba(geoNumLeg[i]))) # rectangle
                ax.annotate(geolegend[i][1], xy=(rx + rw + rw*0.25, ry + rh*0.33)) # text
                ry = ry-(rh/2 + rh)


        print("The detected geological are: ")
        print(geolegend)
        print("Done, if the legend contains too much element, It might be unusable")

    #Plotting the other stuffs
    print("Plotting the rest of the data")
    if(median_leg):
        obj = ax.scatter(flow_distance, elevation, c = mchi,cmap=plt.cm.jet, s = size_point, lw=0, norm=matplotlib.colors.Normalize(vmin=0, vmax=median_mult*np.median(mchi)))
    else:
        obj = ax.scatter(flow_distance, elevation, c = mchi,cmap=plt.cm.jet, s = size_point, lw=0, norm=matplotlib.colors.Normalize(vmin=0, vmax=mchi.max()*reductor))
    # m_Chi colorbar
    cbaxes = fig.add_axes([0.1, 0.03, 0.8, 0.02], frameon = False )
    cb = plt.colorbar(obj, cax = cbaxes,orientation="horizontal")




    print("Saving the plot, It might take times (increase with your dpi)")
    # tidy up the figure
    ax.grid(grid_plot)
    ax.set_title(fig_title)
    ax.set_xlabel(x_label_txt)
    ax.set_ylabel(y_label_txt)
    ax.set_ylim(0-elevation.max()/20,elevation.max()+elevation.max()/20)
    ax.set_xlim(0-flow_distance.max()/20,flow_distance.max()+flow_distance.max()/20)
    if(Chi_plot):
        plt.savefig(wDirectory+wname+"_mChi_Chi_Elevation_2sides"+ext,dpi=pdpi)
    else:
        plt.savefig(wDirectory+wname+"_mChi_FlowDistance_Elevation_2sides"+ext,dpi=pdpi)
    plt.clf()

    ####### separated_legend ######
    if(separated_legend):
        matplotlib.rcParams["font.size"] = 6
        print("I am now printing the geological legend in a separate file (switch off the separated_legend if you don't want it)")
        fig, ax = plt.subplots(figsize=(4, 3)) # creating the new figure that will host the legend
        ax.set_position([0, 0, 1, 1]) # the ax should take the whole figure
        ax.tick_params(axis='both',width =0 )
        rx = 0.1
        ry = 0.9
        rw = 0.05
        rh = 0.03
        for i in range(GeoKey[:,1].shape[0]):
            for h in range(glimolokey.shape[0]):
                if(GeoKey[i][0] == glimolokey[h][0]):
                    ax.add_patch(patches.Rectangle((rx, ry),rw, rh, facecolor = scalarMap.to_rgba(GeoKey[i][1]))) # rectangle
                    ax.annotate(glimolokey[h][1], xy=(rx + rw + rw*0.25, ry + rh*0.33)) # text
                    ry = ry-(rh/2 + rh)
        print("Saving the legend")
        plt.savefig(wDirectory+"Geologic_Legend"+ext,dpi=pdpi)


        print("Done")
else:
    print("Check your file, I cannot find it")
