import numpy as np
from scipy import stats
#import pdb;

def confidence_interval(x, y_pred, x_mean, s_e, t_crit,n, factor=1):
    return t_crit * s_e * np.sqrt(1/n + (x - x_mean)**2 / np.sum((x - x_mean)**2)) * factor
    
def run_method(years, temperature, uncert, model_run, experiment_type):

    data = np.genfromtxt(open("./Common_Data/co2_mlo_icedome.csv", "rb"),dtype=float, delimiter=',')
    datesmlo=data[:,0]
    co2mlo=data[:,1]
    datesicedome=data[:,2]
    co2icedome=data[:,3]
    #create full record Jarvis 2024 paper
    Co2=np.full(np.shape(years),np.nan)
    Co2[int(datesmlo[0]-1850):]=co2mlo
    Co2[:int(datesmlo[0]-1850)]= np.interp(np.arange(1850,datesmlo[0]),datesicedome,co2icedome)
    
    means = np.full(np.shape(years),np.nan)
    ses = np.full(np.shape(years),np.nan)
    #use OLS unlike Jarvis who used weighted least squares

    temp_base = np.mean(temperature[0:30])

    for endd in range(1850,1900): #use half Arrhenius's 5.5°C/doubling estimate, with now 1.5-4.5, so ~ 2°C ses uncertainty
        means[endd-1850] = temp_base + 5.5 /2 * (Co2[endd-1850]-Co2[0])/Co2[0]
        ses[endd-1850] =  3/4  * (Co2[endd-1850]-Co2[0])/Co2[0]
    
    for endd in range(1900,years[-1]+1):
        
        regX = Co2[0: (endd+1 -1850)]
        n=len(regX)
        regY = temperature[0:(endd+1-1850)]
        slope, intercept, _, _, _ = stats.linregress(regX, regY) #can also give stdrr of slope, intercept

        y_pred = slope * regX + intercept
        residuals = regY - y_pred
        # Calculate standard error of the estimate (s_e)
        s_e = np.sqrt(np.sum(residuals**2) / (n - 2))

        # Critical t-value for 1 se confidence level
        t_crit = stats.t.ppf(.5+0.3413, df=n-2) #one standard error

        means[endd-1850] = y_pred[-1]
        cis = confidence_interval(regX, y_pred, np.mean(regX), s_e, t_crit, n )
        ses[endd-1850] = cis[-1]
    print("Fitted line")
    print(slope,intercept)
    means2= y_pred
    ses2= cis
    
    return means, ses, means2,ses2

