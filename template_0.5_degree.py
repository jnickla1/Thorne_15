import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

# Example data: Replace this with actual crossing year data for each method
methods = ['Method 1', 'Method 2', 'Method 3', 'Method 4']
crossing_years = {
    'Method 1': np.random.normal(1985, 3, 100),
    'Method 2': np.random.normal(1987, 4, 100),
    'Method 3': np.random.normal(1986, 2, 100),
    'Method 4': np.random.normal(1988, 3, 100)
}
central_estimates = {method: np.mean(years) for method, years in crossing_years.items()}
standard_errors = {method: np.std(years) / np.sqrt(len(years)) for method, years in crossing_years.items()}

# Set up the figure
fig, (ax, hist_ax)  = plt.subplots(2, 1, gridspec_kw={'height_ratios': [2, 1] ,'left':0.28 , 'right':0.9 , 'top':0.88, 'bottom':0.12}, sharex=True)
ax.axvline(1986, color='gray', linestyle='--', linewidth=2, label="Gold Standard (1986)")

# Plot each method’s violin plot
for i, method in enumerate(methods):
    # Violin plot
    sns.violinplot(
        x=crossing_years[method],
        y=i/10, 
        orient='h', 
        inner=None, 
        color="skyblue", 
        ax=ax,
        linewidth=0.8, 
        width=0.5,
        #scale="width",

    )
    # Central estimate with error bars
    central_year = central_estimates[method]
    se = standard_errors[method]
    ax.errorbar(
        central_year, i, xerr=se, fmt='o', color='blue', capsize=3, label="Central Estimate" if i == 0 else ""
    )
    ax.scatter(central_year, i, color='red', s=50, label="Interpolated Year" if i == 0 else "")

# Formatting
ax.set_yticks(range(len(methods)))
ax.set_yticklabels(methods)
ax.set_xlabel("Year")
ax.set_title("Crossing Years for 0.5°C Above Preindustrial by Method")
ax.legend(loc='upper right')

# Histogram of crossing years at the bottom

all_central_estimates = [central_estimates[method] for method in methods]
hist_ax.hist(all_central_estimates, bins=np.arange(1975, 2000, 1), color="skyblue", edgecolor="gray")
hist_ax.set_ylabel("Count")
hist_ax.set_xlabel("Year")
#hist_ax.set_title("Histogram of Central Estimate Crossing Years")

plt.show()
