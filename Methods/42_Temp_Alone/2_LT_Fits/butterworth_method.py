import numpy as np
from scipy.signal import butter, filtfilt, welch
import os

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




# Plot the smoothed temperature data
#ax.plot(years, smoothed_temp1, color='darkorange', label='Butterworth Smoother')
#ax.plot(years, smoothed_temp3, color='coral', label='Butterworth Smoother')



def run_method(years, temperature, uncert, model_run, experiment_type):


    # Create the Butterworth filter
    order = 3    # Filter order (can be adjusted)
    cutoff = 1/30  # Cutoff frequency for smoothing
    b, a = butter(order, cutoff, btype='low', analog=False)



    dir_path = os.path.dirname(os.path.realpath(__file__))
    weight_file =  dir_path+"/butterworth/butterworth_weights_"+experiment_type+str(model_run)+".npy"
    weights_exist = False
    empser  = np.full(np.shape(years),np.nan)
    means = empser.copy()
    ses = empser.copy()
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs
    st_date =100
    
    if os.path.exists(weight_file):
            # Load the weights from file
        weights = np.load(weight_file)
        weights_exist=True
    else:
        weights = np.empty((len(years)+1-st_date, 3))
    
    min_uncert = np.std(temperature[0:st_date])/np.sqrt(st_date)
    for i in range(st_date, len(years) +1):
            # Apply the filter to the temperature data
        
        smoothed_temp1 = filtfilt(b, a, temperature[0:i])
        smoothed_temp0 = smoothed_temp1+ min_uncert #np.concatenate((smoothed_temp1[0:8],filtfilt(b, a, temperature[8:i])))
        smoothed_temp2,mse = lowpass(temperature[0:i], 1/40, 1,1)
        smoothed_temp3,mse = lowpass(temperature[0:i], 1/40, 2,2)
        
        if (weights_exist):
            smoothed_temp4,_ = lowpassadaptive_w(temperature[0:i], 1/40,weights[i-st_date,0],weights[i-st_date,1],weights[i-st_date,2])
        else:
            smoothed_temp4,weights[i-st_date,0],weights[i-st_date,1],weights[i-st_date,2],_ = lowpassadaptive(temperature[0:i], 1/40)

        means[i-1] = smoothed_temp4[-1]
        ses[i-1] = (max(smoothed_temp4[-1], smoothed_temp3[-1],smoothed_temp2[-1],smoothed_temp1[-1], smoothed_temp0[-1]) -
            min(smoothed_temp4[-1], smoothed_temp3[-1],smoothed_temp2[-1],smoothed_temp1[-1], smoothed_temp0[-1]))

    if not(weights_exist):
        np.save(weight_file, weights)

    smsses = (np.array([smoothed_temp4, smoothed_temp3,smoothed_temp2,smoothed_temp1,smoothed_temp0]).max(axis=0) -
           np.array([smoothed_temp4, smoothed_temp3,smoothed_temp2,smoothed_temp1,smoothed_temp0]).min(axis=0))
    
    return means, ses, smoothed_temp4, smsses

