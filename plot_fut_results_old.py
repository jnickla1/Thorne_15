import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import numpy as np
import pdb

# Initialize lists to store Edyrs15 values for each ssp_id
ssp_126 = []
ssp_245 = []
ssp_370 = []

ssp_126e = []
ssp_245e = []
ssp_370e = []


sel_methods = ['CGWL10y_sfUKCP', 'EBMKF_ta2','Kal_flexLin' ,'removeGreensfx',]
#already sorted
sel_methods_colors = [ "#CC3311", "#999933","#88CCEE" ,"#882255"]



# Define the file pattern and loop through all matching files
file_pattern = "/Results/current_fut_statistics_fut_ESM1-2-LR_SSP{ssp_id}_constVolc{run}.csv"
for ssp_id in [126, 245, 370]:
    for run in range(50):
        file_name = os.getcwd() + file_pattern.format(ssp_id=ssp_id, run=run)
        if os.path.exists(file_name):  # Check if the file exists
            # Read the CSV file
            df = pd.read_csv(file_name)
            
            # Filter rows with the desired method_name
            filtered = df[df['method_name'].isin(sel_methods)]
            df_sorted = filtered.sort_values(by='method_name') 

            # Extract Edyrs15 values and add to the corresponding list
            if ssp_id == 126:
                ssp_126.append(df_sorted['Edyrs15'].values)
                ssp_126e.append(df_sorted['eEdyrs15'].values)
                if (abs(df_sorted.loc[df_sorted['method_name'] == 'EBMKF_ta2', 'Edyrs15'].values[0]) > 10):
                    print(file_name)
            elif ssp_id == 245:
                ssp_245.append(df_sorted['Edyrs15'].values)
                ssp_245e.append(df_sorted['eEdyrs15'].values)
            elif ssp_id == 370:
                ssp_370.append(df_sorted['Edyrs15'].values)
                ssp_370e.append(df_sorted['eEdyrs15'].values)

file_pattern2 = "/Results/current_fut_statistics_fut_NorESM_RCP45_{volc_id}{run}.csv"

v_45 = []
v_45e = []
cv_45 = []
cv_45e = []

for volc_id in ["Volc","VolcConst"]:
    for run in range(50):
        file_name = os.getcwd() + file_pattern2.format(volc_id=volc_id, run=run)
        if os.path.exists(file_name):  # Check if the file exists
            # Read the CSV file
            df = pd.read_csv(file_name)
            
            # Filter rows with the desired method_name
            filtered = df[df['method_name'].isin(sel_methods)]
            df_sorted = filtered.sort_values(by='method_name') 

            # Extract Edyrs15 values and add to the corresponding list
            if volc_id == "Volc":
                v_45.append(df_sorted['Edyrs15'].values)
                v_45e.append(df_sorted['eEdyrs15'].values)
                
            elif volc_id == "VolcConst":
                cv_45.append(df_sorted['Edyrs15'].values)
                cv_45e.append(df_sorted['eEdyrs15'].values)

# Create a figure with 3 subplots for the smoothed histograms
fig, axs = plt.subplots(5, 2, figsize=(10, 10), sharex=True)

# Define a helper function to plot soothed histograms
def plot_smoothed_histogram(ax, data0, title , leg=False):
    for i in range(np.shape(data0)[1]):
        data = data0[:,i]
        kde = gaussian_kde(data)
        x_range = np.linspace(min(data)-1.5, max(data)+1.5, 1000)
        ax.plot(x_range, kde(x_range), color=sel_methods_colors[i], label = sel_methods[i], lw=2)
        #ax.hist(data,density=True,color=sel_methods_colors[i],alpha=0.3)
        #ax.fill_between(x_range, kde(x_range), color=sel_methods_colors[i], alpha=0.2)
    ax.set_title(title, fontsize=14)
    ax.grid(True)
    ax.tick_params(axis='x', labelbottom=True)
    if leg:
        ax.legend(loc='best')

# Plot smoothed histograms for each SSP
plot_smoothed_histogram(axs[2,0], np.array(ssp_126), "MPI SSP126 in-member", )
plot_smoothed_histogram(axs[1,0], np.array(ssp_245), "MPI SSP245 in-member", )
plot_smoothed_histogram(axs[0,0], np.array(ssp_370), "MPI SSP370 in-member",leg=True)
plot_smoothed_histogram(axs[3,0], np.array(cv_45), "NorESM RCP45 ConstVolc in-member",)
plot_smoothed_histogram(axs[4,0], np.array(v_45), "NorESM RCP45 Volc in-member")

plot_smoothed_histogram(axs[2,1], np.array(ssp_126e), "MPI SSP126 ensemble", )
plot_smoothed_histogram(axs[1,1], np.array(ssp_245e), "MPI SSP245 ensemble", )
plot_smoothed_histogram(axs[0,1], np.array(ssp_370e), "MPI SSP370 ensemble",leg=True)
plot_smoothed_histogram(axs[3,1], np.array(cv_45e), "NorESM RCP45 ConstVolc ensemble",)
plot_smoothed_histogram(axs[4,1], np.array(v_45e), "NorESM RCP45 Volc ensemble")

axs[4,0].set_xlim(-10,10)
# Set x-axis label
axs[4,0].set_xlabel("Error in Crossing Instant (years)", fontsize=12)
axs[4,1].set_xlabel("Error in Crossing Instant (years)", fontsize=12)

# Adjust layout and show the plot
plt.tight_layout()
plt.show()
