# small tutorial to use pandas to create box and whisker plots
import matplotlib
matplotlib.use("Agg")

from matplotlib import pyplot as plt
import pandas as pd

# first step: read the file you want the information
## you'll have to adapt the path of your file
df = pd.read_csv("/home/s1675537/PhD/LSDTopoData/Chrystiann_Boris/movern_0.35/corsica_Mchi_Geology.csv")

# your df (dataframe) now contains all the column and data of your csv file. 
# if you need to check the following command print the name of each columns:
print(df.columns.values)

# I assume you have your file with all the geology correspondances ?
# you need to isolate the data bins you want
# imagine I want to plot the m_chi repartition for my lithologies 5, 11 and 13/15 grouped together: I need a list of the arrays

df_5 = df[df["geology"] == 5] # the logic is: df_5 is equal to the entire df WHERE the column "geology" is 5
# df_5 now contains all the information of my original file, but just when the column geology was 5
df_11 = df[df["geology"] == 11] # same logic

df_13_15 = df[df["geology"].isin([13,15])] # the logic is slightly different when you want to check several values


# I want now to gather everything in a "list of arrays" to plot:
data_to_plot = [df_5["m_chi"], df_11["m_chi"], df_13_15["m_chi"]] # I only want the column m_chi
# you want to name it as well, in the same order
names = ["litho 5", "litho 11", "litho 13/15"]

# ok now we are ready to plot:

# Create a figure
fig = plt.figure(1, figsize=(9, 6))

# Create an axes
ax = fig.add_subplot(111)

# Create the boxplot
bp = ax.boxplot(data_to_plot, labels = names)

# Save the figure
fig.savefig('simple_boxplot.png', bbox_inches='tight')
