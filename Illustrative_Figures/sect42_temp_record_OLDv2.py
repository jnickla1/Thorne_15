import matplotlib.pyplot as plt
import numpy as np
# Create a figure and a set of subplots
npanels2=5
fig, axes = plt.subplots(nrows=2, ncols=npanels2, figsize=(14, 7))

# Generate fictitious temperature data (for demonstration)
data = np.genfromtxt(open(os.path.expanduser("~/")+"Dropbox/KalmanFilterPython/toyKFmodelData8c.csv", "rb"),dtype=float, delimiter=',')
years=data[:,0]
years[0]=1850
nyrs = len(years)
temperature = data[:, 1]
all_method_names = ["5yr lagging", "10yr lagging",
                    "OLS Fit", "Theil Sen", "Hinge Fit","Quartic","Bayes CP",
                    "11yr offset", "15yr fit end","LOWESS-linear", "Butterworth"
                    "Autoregression",
                    "Remove El Nino", "Remove Greens Functions"]

impl_methods = ["20yr lagging", "OLS Fit", "Hinge Fit", "11yr offset (Trewin)", "15yr-linear (SR1.5)","LOWESS-linear ??cut??","Butterworth", "Autoregression (AR6)", "State Space", "Remove El Nino (Foster)"]

# Loop over each subplot (method) to configure the layout and base plot
for i, ax in enumerate(axes.flat):
    # Scatter plot of the base temperature data
    ax.scatter(years, temperature, color='lightgrey', label='Temperature Data',s=2)

    # Hide the x-axis labels but keep the axis visible
    ax.set_xticks([])
    
    # Only label the y-axis on the left-most column
    if i % npanels2 == 0:
        ax.set_ylabel('GMST Anomaly (°C)')
    else:
        ax.set_yticks([])  # Hide y-axis labels for non-leftmost panes
    
    # Set a title for each method (to be customized later)
    ax.set_title(impl_methods[i])



ax = axes[0, 0]  # First pane
lag_len=20
for i in range(lag_len, len(years)):
    lag_mean = np.mean(temperature[i-lag_len:i])
    color = 'mediumorchid' #if (i - 5) % 5 == 1 else 'thistle'  # Darker every 5th line
    if (i - lag_len) % lag_len == 3:
        ax.plot([years[i-lag_len], years[i]], [lag_mean, lag_mean], color=color)
        ax.plot(years[i], lag_mean, 's', color=color, markersize=6)
    else:
        ax.plot(years[i], lag_mean, 's', color=color, markersize=1)

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
breakyr = 1974-1850
model = LinearRegression().fit(X[breakyr:], temperature[breakyr:])
trendline = model.predict(X[breakyr:])

# Plot the OLS trendline
ax.plot(years[breakyr:], trendline, color='blue', label='OLS Trendline')

model = LinearRegression().fit(X[0:breakyr:], temperature[0:breakyr])
trendline = model.predict(X[0:breakyr])
ax.plot(years[0:breakyr], model.predict(X[0:breakyr]), color='darkgrey', label='OLS Trendline')


ax = axes[0, 3]  # Fourth pane - 11yr offset
breakyr = 1970-1850
offset = 0.11
for i in range(breakyr, len(years)):
    lag_mean = np.mean(temperature[i-11:i])
    color = 'orange' #if (i - 5) % 5 == 1 else 'thistle'  # Darker every 5th line
    if (i) % 11 == 5:
        ax.plot([years[i-11], years[i]], [lag_mean, lag_mean], color=color)
        ax.plot([years[i], years[i]], [lag_mean, lag_mean+offset], color=color)
        ax.plot(years[i], lag_mean+offset, 's', color=color, markersize=6)
    else:
        ax.plot(years[i], lag_mean+offset, 's', color=color, markersize=1)


ax = axes[0, 4]  # Fifth pane - SR1.5 endpoint of 15 yr trend
for i in range(15, len(years)+1):
    # Fit the OLS trendline for each 15-year chunk
    X_chunk = years[i-15:i].reshape(-1, 1)
    y_chunk = temperature[i-15:i]
    model = LinearRegression().fit(X_chunk, y_chunk)
    trendline_chunk = model.predict(X_chunk)
    color='magenta'
    if (i) % 15 == 5:
    # Plot the trendline
        ax.plot(years[i-15:i], trendline_chunk, color='magenta')
        #ax.plot(years[i-1], trendline_chunk[-1], '^', color='magenta', markersize=6)
        
        # Plot the arrow using the first and last points of the chunk
        ax.annotate('', xy=(years[i-1], trendline_chunk[-1]), xytext=(years[i-15], trendline_chunk[0]), 
                arrowprops=dict(facecolor=color, edgecolor=color, arrowstyle='-|>', lw=2,shrinkA=0, shrinkB=0))
    else:
        ax.plot(years[i-1], trendline_chunk[-1], '^', color='magenta', markersize=1)
    


from scipy.signal import butter, filtfilt

# Method: Butterworth Smoother
ax = axes[1, 1] 

# Define the Butterworth filter

def lowpass(indata, frequency, iconb, icone):
    # Define filter order and normalized frequency
    ipts = 10  # 10-point Butterworth filter
    fn = frequency * 2  # Normalized frequency (cycles per time unit)

    # Length of input data
    nn = len(indata)

    # Calculate the amount of padding
    npad = 3 * round(1 / fn)

    # Initialize padded array with zeros
    padded = np.zeros(nn + 2 * npad)

    # Fill the padded array with the input data
    padded[npad:npad + nn] = indata
    padded[npad + nn:] = indata[-1]
    padded[:npad] = indata[0]

    # Implement the boundary constraints for the right-hand side
    for j in range(nn + npad, nn + 2 * npad):
        ipad = j - nn - npad
        if icone == 0:
            apad = np.mean(indata)  # Minimum norm
        elif icone == 1:
            apad = indata[nn - ipad - 1]  # Minimum slope
        else:
            apad = 2 * indata[-1] - indata[nn - ipad - 1]  # Minimum roughness
        padded[j] = apad

    # Implement the boundary constraints for the left-hand side
    for j in range(npad):
        ipad = j
        if iconb == 0:
            apad = np.mean(indata)  # Minimum norm
        elif iconb == 1:
            apad = indata[npad - ipad - 1]  # Minimum slope
        else:
            apad = 2 * indata[0] - indata[npad - ipad - 1]  # Minimum roughness
        padded[j] = apad

    # Apply the Butterworth lowpass filter
    b, a = butter(ipts, fn, btype='low')
    smoothed0 = filtfilt(b, a, padded)

    # Extract the smoothed portion (removing the padded areas)
    smoothed = smoothed0[npad:npad + nn]

    # Calculate the mean squared error (MSE) relative to the original data
    resid = smoothed - indata
    mse = np.var(resid) / np.var(indata)

    return smoothed, mse

def lowpassmin(indata, frequency):
    mse0 = 999.0
    smoothed0 = None
    icb = 0
    ice = 0
    
    # Iterate over boundary condition combinations
    for iconb in range(3):
        for icone in range(3):
            # Apply the lowpass filter with the current boundary conditions
            smoothed, mse = lowpass(indata, frequency, iconb, icone)
            
            # Check if the current mse is smaller than the best mse
            if mse <= mse0:
                icb = iconb
                ice = icone
                mse0 = mse
                smoothed0 = smoothed

    return smoothed0, icb, ice, mse0

def lowpassadaptive(indata, frequency):
    msebest = 1e+15

    # First, optimize the lower constraint
    smoothedlower, icb, ice, mse0 = lowpassmin(indata, frequency)

    # Now, perform adaptive optimization for the upper constraint
    smoothed0, mse0 = lowpass(indata, frequency, icb, 0)
    smoothed1, mse1 = lowpass(indata, frequency, icb, 1)
    smoothed2, mse2 = lowpass(indata, frequency, icb, 2)

    # Search for the best combination of weights
    for weight0 in np.arange(0, 1.01, 0.01):
        for weight1 in np.arange(0, 1.01 - weight0, 0.01):
            weight2 = 1 - weight0 - weight1
            smoothed = weight0 * smoothed0 + weight1 * smoothed1 + weight2 * smoothed2
            mse = np.var(smoothed - indata) / np.var(indata)
            if mse <= msebest:
                w0 = weight0
                w1 = weight1
                w2 = weight2
                msebest = mse
                smoothedbest = smoothed

    return smoothedbest, w0, w1, w2, msebest



# Create the Butterworth filter
order = 3    # Filter order (can be adjusted)
cutoff = 1/30  # Cutoff frequency for smoothing
b, a = butter(order, cutoff, btype='low', analog=False)

# Apply the filter to the temperature data
smoothed_temp1 = filtfilt(b, a, temperature)

smoothed_temp2,mse = lowpass(temperature, 1/40, 1,1)
smoothed_temp3,mse = lowpass(temperature, 1/40, 2,2)


# Plot the smoothed temperature data
#ax.plot(years, smoothed_temp1, color='darkorange', label='Butterworth Smoother')
#ax.plot(years, smoothed_temp3, color='coral', label='Butterworth Smoother')

for i in range(100, len(years)+1):
    smoothed_temp4,_,_,_,_ = lowpassadaptive(temperature[0:i], 1/40)
    if i%40==14:
        ax.plot(years[0:i], smoothed_temp4, color='orange')
        ax.plot(years[i-1], smoothed_temp4[-1],'o', color='orange',markersize=6)
    else:
        ax.plot(years[i-1], smoothed_temp4[-1],'o', color='orange',markersize=1)




# Add a legend outside the panes (customize legend entries later)
handles, labels = ax.get_legend_handles_labels()
#fig.legend(handles, labels, loc='upper right')

# Adjust layout to minimize spacing between panes
plt.tight_layout(rect=[0, 0, 0.9, 1])
#plt.savefig("sec42.pdf")
# Display the figure
plt.show()
