import matplotlib.pyplot as plt
import numpy as np
# Create a figure and a set of subplots
fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(12, 7))

# Generate fictitious temperature data (for demonstration)
np.random.seed(0)
nyrs=24
years = np.arange(0, nyrs)  # Example time range (to be hidden)
temperature = np.concatenate([np.zeros(int(nyrs/2)),np.linspace(0,1,int(nyrs/2))])
temperature = temperature + np.random.normal(0, 0.1, len(years))

all_method_names = ["5yr lagging", "10yr lagging",
                    "OLS Fit", "Theil Sen", "Hinge Fit","Quartic","Bayes CP",
                    "11yr offset", "15yr fit end","LOWESS-linear", "Butterworth"
                    "Autoregression",
                    "Remove El Nino", "Remove Greens Functions"]

impl_methods = ["5yr lagging", "OLS Fit", "Hinge Fit", "11yr offset", "LOWESS-linear","Butterworth", "Autoregression", "Remove El Nino"]

# Loop over each subplot (method) to configure the layout and base plot
for i, ax in enumerate(axes.flat):
    # Scatter plot of the base temperature data
    ax.scatter(years, temperature, color='lightgrey', label='Temperature Data')

    # Hide the x-axis labels but keep the axis visible
    ax.set_xticks([])
    
    # Only label the y-axis on the left-most column
    if i % 4 == 0:
        ax.set_ylabel('GMST Anomaly (°C)')
    else:
        ax.set_yticks([])  # Hide y-axis labels for non-leftmost panes
    
    # Set a title for each method (to be customized later)
    ax.set_title(impl_methods[i])



ax = axes[0, 0]  # First pane
for i in range(5, len(years)):
    lag_mean = np.mean(temperature[i-5:i])
    color = 'mediumorchid' #if (i - 5) % 5 == 1 else 'thistle'  # Darker every 5th line
    if (i - 5) % 5 == 3:
        ax.plot([years[i-5], years[i]], [lag_mean, lag_mean], color=color)
        ax.plot(years[i], lag_mean, 's', color=color, markersize=6)
    else:
        ax.plot(years[i], lag_mean, 's', color=color, markersize=3)

from sklearn.linear_model import LinearRegression

# Method 2: Ordinary Least Squares (OLS) Trendline
ax = axes[0, 1]  # Second pane

# Fit the OLS trendline
X = years.reshape(-1, 1)
model = LinearRegression().fit(X, temperature)
trendline = model.predict(X)

# Plot the OLS trendline
ax.plot(years, trendline, color='deepskyblue', label='OLS Trendline')


# Method 3: Ordinary Least Squares (OLS) Trendline
ax = axes[0, 2]  # Second pane

# Fit the OLS trendline
X = years.reshape(-1, 1)
model = LinearRegression().fit(X[12:], temperature[12:])
trendline = model.predict(X[12:])

# Plot the OLS trendline
ax.plot(years[12:], trendline, color='blue', label='OLS Trendline')
ax.plot(years[0:12], np.zeros(12), color='darkgrey', label='OLS Trendline')


ax = axes[0, 3]  # Fourth pane
offset = 50/nyrs*2 * 0.11
for i in range(11, len(years)):
    lag_mean = np.mean(temperature[i-11:i])
    color = 'orange' #if (i - 5) % 5 == 1 else 'thistle'  # Darker every 5th line
    if (i - 11) % 11 == 1:
        ax.plot([years[i-11], years[i]], [lag_mean, lag_mean], color=color)
        ax.plot([years[i], years[i]], [lag_mean, lag_mean+offset], color=color)
        ax.plot(years[i], lag_mean+offset, 's', color=color, markersize=6)
    else:
        ax.plot(years[i], lag_mean+offset, 's', color=color, markersize=3)


# Add a legend outside the panes (customize legend entries later)
handles, labels = ax.get_legend_handles_labels()
#fig.legend(handles, labels, loc='upper right')

# Adjust layout to minimize spacing between panes
plt.tight_layout(rect=[0, 0, 0.9, 1])

# Display the figure
plt.show()
