import pandas as bamboo_bears
import numpy as np
from matplotlib import pyplot as plt

# Importing the data
df = bamboo_bears.read_csv("sieve_2.csv", sep =",") #reading the file
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
plt.show()
