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

sel_methods = ['CGWL10y_sfUKCP', 'EBMKF_ta3',]
sel_methods_colors = [ "#CC3311", "#999933" ]

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
                if (df_sorted['Edyrs15'].values == 15).any():
                    print(file_name)
            elif ssp_id == 245:
                ssp_245.append(df_sorted['Edyrs15'].values)
            elif ssp_id == 370:
                ssp_370.append(df_sorted['Edyrs15'].values)

# Create a figure with 3 subplots for the smoothed histograms
fig, axs = plt.subplots(3, 1, figsize=(6, 10), sharex=True)

# Define a helper function to plot smoothed histograms
def plot_smoothed_histogram(ax, data0, title , leg=False):
    for i in range(np.shape(data0)[1]):
        data = data0[:,i]
        kde = gaussian_kde(data)
        x_range = np.linspace(min(data), max(data), 1000)
        ax.plot(x_range, kde(x_range), color=sel_methods_colors[i], label = sel_methods[i], lw=2)
        #ax.fill_between(x_range, kde(x_range), color=sel_methods_colors[i], alpha=0.2)
    ax.set_title(title, fontsize=14)
    ax.grid(True)
    if leg:
        ax.legend(loc='best')

# Plot smoothed histograms for each SSP
plot_smoothed_histogram(axs[0], np.array(ssp_126), "SSP126", )
plot_smoothed_histogram(axs[1], np.array(ssp_245), "SSP245", )
plot_smoothed_histogram(axs[2], np.array(ssp_370), "SSP370", leg=True)

# Set x-axis label
axs[2].set_xlabel("Error in Crossing Instant", fontsize=14)

# Adjust layout and show the plot
plt.tight_layout()
plt.show()
