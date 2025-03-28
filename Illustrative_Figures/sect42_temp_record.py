import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# Create a figure and a set of subplots
npanels2=5
fig, axes = plt.subplots(nrows=2, ncols=npanels2, figsize=(12, 7))#expand wi
# Generate fictitious temperature data (for demonstration)
data = np.genfromtxt(open("../Common_Data/toyKFmodelData8c.csv", "rb"),dtype=float, delimiter=',')
data2 = pd.read_csv("../Common_Data/HadCRUT.5.0.2.0.analysis.summary_series.global.annual.csv")
years= data2["Time"].to_numpy()
#years[0]=1850
nyrs = len(years)
temperature = data2["Anomaly (deg C)"].to_numpy()
offset_t = np.mean(temperature[0:51]) #1850-1900
temperature = temperature - offset_t
all_method_names = ["5yr lagging", "10yr lagging",
                    "OLS Fit", "Theil Sen", "Hinge Fit","Quartic","Bayes CP",
                    "11yr offset", "15yr fit end","LOWESS-linear", "Butterworth"
                    "Autoregression",
                    "Remove El Nino", "Remove Greens Functions"]

impl_methods = ["30yr (20yr) centred", "10yr (20yr) lagging","OLS Fit", "Hinge Fit", "11yr offset (Trewin)", "15yr-linear (SR1.5)","Butterworth Adapt. (Mann)",
                "GAM AR1 (Stephenson)", "Kalman: Random Walk", "Remove ENSO (Foster)"]

# Loop over each subplot (method) to configure the layout and base plot
for i, ax in enumerate(axes.flat):
    # Scatter plot of the base temperature data
    ax.scatter(years, temperature, color='lightgrey', label='HadCRUT5, annual',s=2)

    # Hide the x-axis labels but keep the axis visible
    if i >= npanels2:
        ax.set_xlabel('Year')
    else:
        ax.set_xticklabels([]) 
    
    # Only label the y-axis on the left-most column
    if i % npanels2 == 0:
        ax.set_ylabel('GMST Anomaly (°C)\nRelative to 1850-1900')
    else:
        ax.set_yticklabels([])  # Hide y-axis labels for non-leftmost panes
    
    # Set a title for each method (to be customized later)
    ax.set_title(impl_methods[i])


ax = axes[0, 0]  # First pane
avg_len_l=15
avg_len_u=15 #actually 15 yrs below, that year, 14 yrs after
for i in range(avg_len_l, len(years) - avg_len_u):
    lag_mean = np.mean(temperature[i-avg_len_l:i+avg_len_u])
    color = 'black' 
    if (i ) % (avg_len_l+avg_len_u) == 9:
        ax.plot([years[i-avg_len_l], years[i+avg_len_u-1]], [lag_mean, lag_mean], color=color)
        ax.plot(years[i], lag_mean, '|', color=color, markersize=7)
    else:
        ax.plot(years[i], lag_mean, '|', color=color, markersize=3)



avg_len_l=10
avg_len_u=10 #actually 10 yrs below, that year, 9 yrs after
for i in range(avg_len_l, len(years) - avg_len_u):
    chunka=temperature[i-avg_len_l:i+avg_len_u]
    chunkb=temperature[i-avg_len_l+1:i+avg_len_u+1]
    lag_mean = np.mean([chunka,chunkb])
    color = 'black'
    if (i ) % (avg_len_l+avg_len_u) == 4:
        ax.plot([years[i-avg_len_l], years[i+avg_len_u]], [lag_mean, lag_mean], color=color, alpha=0.5)
        ax.plot(years[i], lag_mean, '|', color=color, markersize=6, alpha=0.5)
    else:
        ax.plot(years[i], lag_mean, '|', color=color, markersize=2, alpha=0.5)

axlim =ax.get_xlim()
ax.plot(0, lag_mean, '|', color=color, markersize=6, label="30yr mean")
ax.plot(0, lag_mean, '|', color=color, markersize=4,alpha=0.5, label="20yr mean (*)")
ax.legend()
ax.set_xlim(axlim)

ax = axes[0, 1]  
lag_len=10
for i in range(lag_len, len(years)):
    lag_mean = np.mean(temperature[i-lag_len:i])
    color = 'violet' #if (i - 5) % 5 == 1 else 'thistle'  # Darker every 5th line
    if (i - lag_len) % lag_len == 4:
        ax.plot([years[i-lag_len], years[i]], [lag_mean, lag_mean], color=color)
        ax.plot(years[i], lag_mean, 's', color=color, markersize=6)
    else:
        ax.plot(years[i], lag_mean, 's', color=color, markersize=1)

lag_len=20
for i in range(lag_len, len(years)):
    lag_mean = np.mean(temperature[i-lag_len:i])
    color = 'violet' #if (i - 5) % 5 == 1 else 'thistle'  # Darker every 5th line
    if (i - lag_len) % lag_len == 14:
        ax.plot([years[i-lag_len], years[i]], [lag_mean, lag_mean], color=color,alpha=0.3)
        ax.plot(years[i], lag_mean, 's', color=color, markersize=4,alpha=0.3)
    else:
        ax.plot(years[i], lag_mean, 's', color=color, markersize=1,alpha=0.3)

axlim =ax.get_xlim()
ax.plot(0, lag_mean, 's', color=color, markersize=6, label="10yr lagging")
ax.plot(0, lag_mean, 's', color=color, markersize=4,alpha=0.3, label="20yr lagging")
ax.legend()
handles, labels = ax.get_legend().legend_handles, [text.get_text() for text in ax.get_legend().get_texts()]
del handles[0]
del labels[0]
ax.get_legend()._legend_box = None
ax.legend(handles, labels)
ax.set_xlim(axlim)

from sklearn.linear_model import LinearRegression

# Method 2: Ordinary Least Squares (OLS) Trendline
ax = axes[0, 2]  # Second pane

# Fit the OLS trendline
X = years[:-1].reshape(-1, 1)
model = LinearRegression().fit(X, temperature[:-1])
trendline = model.predict(X)

# Plot the OLS trendline
ax.plot(years[:-1], trendline, color='deepskyblue', label='OLS Trendline')


# Method 3: Ordinary Least Squares (OLS) Trendline
ax = axes[0, 3]  # Second pane

# Fit the OLS trendline
X = years.reshape(-1, 1)
breakyr = 1974-1850
model = LinearRegression().fit(X[breakyr:], temperature[breakyr:])
trendline = model.predict(X[breakyr:])

# Plot the OLS trendline
ax.plot(years[breakyr:], trendline, color='blue', label='OLS Trendline')
ax.plot(years[i], trendline[-1], '*', color='blue', markersize=5)
bk2 = 1984-1850
for i in range(bk2, len(years)):
    # Fit the OLS trendline for each 15-year chunk
    X_chunk = years[breakyr:i].reshape(-1, 1)
    y_chunk = temperature[breakyr:i]
    model = LinearRegression().fit(X_chunk, y_chunk)
    trendline_end = model.predict(X_chunk)
    ax.plot(years[i], trendline_end[-1], '*', color='blue', markersize=2)

model = LinearRegression().fit(X[0:breakyr:], temperature[0:breakyr])
trendline = model.predict(X[0:breakyr])
ax.plot(years[0:breakyr], model.predict(X[0:breakyr]), color='darkgrey', label='OLS Trendline')


ax = axes[0, 4]  # Fourth pane - 11yr offset
breakyr = 1970-1850
offset = 0.11
for i in range(breakyr, len(years)):
    lag_mean = np.mean(temperature[i-11:i])
    color = 'goldenrod' #if (i - 5) % 5 == 1 else 'thistle'  # Darker every 5th line
    if (i) % 11 == 9:
        ax.plot([years[i-11], years[i]], [lag_mean, lag_mean], color=color)
        ax.plot([years[i], years[i]], [lag_mean, lag_mean+offset], color=color)
        ax.plot(years[i], lag_mean+offset, 's', color=color, markersize=6)
    else:
        ax.plot(years[i], lag_mean+offset, 's', color=color, markersize=1)


ax = axes[1, 0]  # Fifth pane - SR1.5 endpoint of 15 yr trend
for i in range(15, len(years)+1):
    # Fit the OLS trendline for each 15-year chunk
    X_chunk = years[i-15:i].reshape(-1, 1)
    y_chunk = temperature[i-15:i]
    model = LinearRegression().fit(X_chunk, y_chunk)
    trendline_chunk = model.predict(X_chunk)
    color='magenta'
    if (i) % 15 == 10:
    # Plot the trendline
        ax.plot(years[i-15:i], trendline_chunk, color='magenta')
        #ax.plot(years[i-1], trendline_chunk[-1], '^', color='magenta', markersize=6)
        
        # Plot the arrow using the first and last points of the chunk
        ax.annotate('', xy=(years[i-1], trendline_chunk[-1]), xytext=(years[i-15], trendline_chunk[0]), 
                arrowprops=dict(facecolor=color, edgecolor=color, arrowstyle='-|>', lw=2,shrinkA=0, shrinkB=0))
    else:
        ax.plot(years[i-1], trendline_chunk[-1], '^', color='magenta', markersize=1)
    


from scipy.signal import butter, filtfilt, welch
import os
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


def lowpassadaptive_w(indata, frequency,weight0,weight1,weight2):
    msebest = 1e+15

    # First, optimize the lower constraint
    smoothedlower, icb, ice, mse0 = lowpassmin(indata, frequency)

    # Now, perform adaptive optimization for the upper constraint
    smoothed0, mse0 = lowpass(indata, frequency, icb, 0)
    smoothed1, mse1 = lowpass(indata, frequency, icb, 1)
    smoothed2, mse2 = lowpass(indata, frequency, icb, 2)
    smoothed = weight0 * smoothed0 + weight1 * smoothed1 + weight2 * smoothed2
    mse = np.var(smoothed - indata) / np.var(indata)
    return smoothed, mse


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
weight_file = "../Methods/42_Temp_Alone/2_LT_Fits/butterworth/butterworth_weights.npy"
weights_exist = False
st_date =100
curcol='darkorange'
if os.path.exists(weight_file):
        # Load the weights from file
    weights = np.load(weight_file)
    weights_exist=True
else:
    weights = np.empty((len(years)+1-st_date, 3))
        
for i in range(st_date, len(years)+1):
    if (weights_exist):
        smoothed_temp4,_ = lowpassadaptive_w(temperature[0:i], 1/40,weights[i-st_date,0],weights[i-st_date,1],weights[i-st_date,2])
    else:
        smoothed_temp4,weights[i-st_date,0],weights[i-st_date,1],weights[i-st_date,2],_ = lowpassadaptive(temperature[0:i], 1/40)
    if i%40==15:
        ax.plot(years[0:i], smoothed_temp4, color=curcol)
        ax.plot(years[i-1], smoothed_temp4[-1],'o', color=curcol,markersize=6)
    else:
        ax.plot(years[i-1], smoothed_temp4[-1],'o', color=curcol,markersize=1)

if not(weights_exist):
    np.save(weight_file, weights)

ylim=ax.get_ylim()

# Inset plot for frequency domain (Power Spectral Density)
ax_inset = ax.inset_axes([0.15, 0.55, 0.5, 0.3])
fs=1
# Frequency domain representation (using Welch's method to estimate power spectral density)
f_original, Pxx_original = welch(temperature, fs, nperseg=20)
f_filtered, Pxx_filtered = welch(smoothed_temp4, fs, nperseg=20)

# Plot original and filtered signal power in frequency domain
ax_inset.plot(f_original, Pxx_original, label="Orig.", color='gray', alpha=0.6)
ax_inset.plot(f_filtered, Pxx_filtered, label="Filter", color=curcol, lw=2)
ax_inset.set_title("Frequency", fontsize=10)
#ax_inset.set_xlabel("Frequency (Hz)", fontsize=8)
#ax_inset.set_ylabel("Power", fontsize=8)
#ax_inset.set_yticklabels([])
ax_inset.tick_params(axis='both',labelsize=7)
ax_inset.legend(loc="upper right", fontsize=7)

# Annotate the inset plot as "Removing High-Frequency Noise"
#ax_inset.annotate("Removing High-Frequency Noise", xy=(3, 0.01), xytext=(5, 0.02),
#                  arrowprops=dict(facecolor='black', arrowstyle='->'), fontsize=8)



# Method: AR1
ax = axes[1, 2]
fits_df = pd.read_csv("../Methods/42_Temp_Alone/4_GAM_AR1/GAM_AR_Stephenson/gamAR1_fits_historical.csv")
st_date = 30
for i in range(0, len(years)-st_date+1):
    if (i+st_date)%40==15:
        ax.plot(years[0:i+st_date], fits_df.iloc[i,0:(st_date+i)]-offset_t, color='olive')
        ax.plot(years[i+st_date-1], fits_df.iloc[i,(st_date+i-1)]-offset_t,'d', color='olive',markersize=6)
    else:
        ax.plot(years[i+st_date-1],fits_df.iloc[i,(st_date+i-1)]-offset_t,'d', color='olive',markersize=1)


from statsmodels.tsa.stattools import acf
lastfitar1 = fits_df.iloc[-1,:]
# Load the CSV files with fitted values and basis functions
lastfitar0= pd.read_csv("../Methods/42_Temp_Alone/4_GAM_AR1/GAM_AR_Stephenson/gam_AR0_last_fit.csv")
basis_functions = pd.read_csv("../Methods/42_Temp_Alone/4_GAM_AR1/GAM_AR_Stephenson/gam_modelmatrix.csv")
ax_inset = ax.inset_axes([0.15, 0.55, 0.5, 0.3])
acf_resid0 = acf(lastfitar0.values[:,0]-temperature)
acf_resid1 = acf(lastfitar1.values-temperature)
bar_width = 0.35
x_pos = np.arange(1,7)
ax_inset.bar(x_pos - bar_width/2, acf_resid0[1:7], bar_width, color="limegreen", label='AR0')
ax_inset.bar(x_pos + bar_width/2, acf_resid1[1:7], bar_width, color="olive", label='AR1')
ax_inset.legend(loc="upper right", fontsize=7)
ax_inset.set_title("Residual Autocorr.", fontsize=10)
ax_inset.tick_params(axis='both',labelsize=7)
ax_inset.set_xticks(x_pos)

ax.plot(years,basis_functions/16-.15,alpha=0.3)
ax.text(1950,-0.15,"Basis Funct. (10)",size=7,alpha=0.5)
ax.set_ylim(ylim)

# Method: State Space Kalman Filter
ax = axes[1, 3]
import sys
sys.path.append("/Users/JohnMatthew/Downloads/Thorne_15_codefigurestats/Methods/42_Temp_Alone/5_Kalman/")
import Kalman_EM_linRW as KF
ax.scatter(years[0:KF.n_iter], KF.OceanRec-offset_t, color='skyblue', label='Ocean SST',s=2,alpha=0.4)
ax.text(1950,0,"Ocean Surf.\n(HADSST4)",size=7,color='skyblue',alpha=0.8)
curcol = "mediumseagreen"
st_date=0
for i in range(0, len(years)-st_date-1):
    if (i+st_date)==len(years)-2:
        ax.plot(years[0:i+st_date], KF.xhathat[0:(st_date+i),0,0]-offset_t, color=curcol)
        ax.plot(years[i+st_date], KF.xhat[(st_date+i),0,0]-offset_t,'P', color=curcol,markersize=6)
    else:
        ax.plot(years[i+st_date],KF.xhat[(st_date+i),0,0]-offset_t,'+', color=curcol,markersize=3)
        

# Data for the pie chart
#xhat[k] = xhatminus[k]+np.matmul(K[k],eta[k])
lnd_weight = np.mean(KF.K[5:,0,0])
sea_weight = np.mean(KF.K[5:,0,1]) #*0.9/0.7
sizes = [lnd_weight, sea_weight, 1-lnd_weight-sea_weight]
labels = ['Land', 'SST', 'Persist']
colors = ['green', 'skyblue', 'mediumseagreen']

# Create a figure with subplots
ax_inset = ax.inset_axes([0.15, 0.55, 0.5, 0.3])
# Create a pie chart in the subplot
ax_inset.pie([lnd_weight+sea_weight,sizes[2]], colors=['lightgrey',colors[2]], startangle=180,radius=1.2)
wedges, texts = ax_inset.pie(sizes, colors=colors, startangle=180)#,wedgeprops=dict(width=0.4))
for i, wedge in enumerate(wedges):
    angle = (wedge.theta2 + wedge.theta1) / 2
    x = wedge.r * 0.6 * np.cos(np.radians(angle))
    y = wedge.r * 0.6 * np.sin(np.radians(angle))
    ax_inset.text(x, y, labels[i], horizontalalignment='center', verticalalignment='center', fontsize=7)
ax_inset.text(-0.8, -0.8, "GMST", horizontalalignment='center', verticalalignment='center', fontsize=7,color="grey")

# Equal aspect ratio ensures that pie is drawn as a circle
ax_inset.set_title("Forward Weight", fontsize=10)
ax_inset.axis('equal')


# Method: Remove ENSO
ax = axes[1, 4]
import statsmodels.api as sm
enso = pd.read_csv("../Common_Data/meiv2.csv")
omega = 2 * np.pi / 1 
cos_2yr = np.cos(0.5 * omega * years)
X = pd.DataFrame({
    'MEI': enso['Avg'][108:-1], #1971 start
    'AOD': data[129:, 3],
    'TSI': data[129:, 8],
    'cos2': cos_2yr[129:-2],
    'time': years[129:-2]-1971,
})
# Add a constant term for the intercept
X = sm.add_constant(X)
y = temperature[129:-2]

model = sm.OLS(y, X).fit()
print(model.params)

X2 = pd.DataFrame({
    'MEI': enso['Avg'][1:-1],
    'AOD': data[22:, 3],
    'TSI': data[22:, 8],
    'cos2': cos_2yr[22:-2],
    'time': np.zeros(len(years)-22-2)
})
X2 = sm.add_constant(X2)
pred2=model.predict(X2)
offset = np.mean(temperature[22:]) - np.mean(pred2)
coefMEI=model.params['MEI'] #0.1

##X0 = pd.DataFrame({
##    'MEI': enso['AVG'][1:],
##    #'AOD': data[22:, 3],
##    #'TSI': data[22:, 8],
##    'cos2': cos_2yr[22:],
##    'time': years[22:]-1971
##})
##X_expand = pd.concat([X0,enso.iloc[1:,1:13]], axis=1)
##from mlxtend.feature_selection import SequentialFeatureSelector
##from sklearn import linear_model
##sfs = SequentialFeatureSelector(linear_model.LinearRegression(),
##                                k_features=4,
##                                forward=True,
##                                scoring='neg_mean_squared_error',
##                                cv=None)
##modelexp =sfs.fit( X_expand,temperature[22:])
##print(modelexp.k_feature_names_)

ax.plot(years[22:-1],temperature[22:-1]-enso['Avg'][1:]*coefMEI, color = "red") #

yenso= enso['Avg'][1:]*coefMEI+0.75-offset_t
ax.fill_between(years[22:-1], yenso, 0.75-offset_t, where=(enso['Avg'][1:] > 0), facecolor='red', alpha=0.6, interpolate=True)  # Red above y=0
ax.fill_between(years[22:-1], yenso, 0.75-offset_t, where=(enso['Avg'][1:] < 0), facecolor='blue', alpha=0.6, interpolate=True)  # Blue below y=0
ax.text(1850,1-offset_t,"Multivariate ENSO Index",verticalalignment='center',size=10)
ax.set_ylim(ylim)

#ax.plot(years[22:],(X2['TSI']-np.mean(X2['TSI']))*model.params['TSI']+0.9, color='goldenrod', alpha=0.6)
#ax.plot(years[22:],X2['AOD']*model.params['AOD']+1, color='teal', alpha=0.6) 
print(model.params)
#ax.plot(years[22:],temperature[22:]-pred2-offset, color = "pink",lw=1) 
# Add a legend outside the panes (customize legend entries later)
handles, labels = ax.get_legend_handles_labels()
#fig.legend(handles, labels, loc='upper right')



# Adjust layout to minimize spacing between panes
plt.tight_layout(rect=[0, 0, 0.9, 1])
plt.savefig("sec42.pdf")
# Display the figure
plt.show()
