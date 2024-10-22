#from Livezey et. al. 2007 https://doi.org/10.1175/2007JAMC1666.1

import numpy as np

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar



# Dummy data for `dates` and `temps`
dates = np.arange(1850, 2021)  # Replace with actual date values
temps = np.random.normal(loc=15, scale=2, size=len(dates))  # Replace with actual temperature values

def OCNerr(N, g, b, t):
    # function for optimization
    nerr= (1+g)/(1+g+(N-1)*(1-g));
    nerr = nerr + (b*((N-1)/2+t))**2;
    return nerr


def run_method(dates, temps, uncert, model_run, experiment_type):

    empser = np.full(np.shape(dates),np.nan)
    means = empser.copy()
    ses = empser.copy()
    meansf = empser.copy()
    sesf = empser.copy()
    (temps_CIl, temps_CIu) = uncert
    temps_1std = (temps_CIu - temps_CIl) / 4 #pm2 stdevs
    
    # Initialization
    et50 = 1850
    computing = 2
    Nguess1 = 30
    sp = 0  # start
    ep = Nguess1  # end
    tau = 0
    Noptimal = 30
    Nguess = Nguess1
    std_guess = np.std(temps[0:60])
    plotting = False
    while computing >= 1:
        #print(sp)
        # Generate trend, sigma, beta, g for this guess interval Nguess=30
        p = np.polyfit(dates[sp:ep], temps[sp:ep], 1)
        model = p[0] * dates[sp:ep] + p[1]
        sigma = np.sqrt(np.sum((model - temps[sp:ep]) ** 2) / (len(model) - 2))  # regression std error
        beta = p[0] / sigma
        g = np.dot(-(model[:-1] - temps[sp:ep-1]), -(model[1:] - temps[sp+1:ep])) / np.sum((model - temps[sp:ep]) ** 2)
        
        # Plot the model
        if(plotting):
            plt.plot(dates[sp:ep], model, linestyle='-', linewidth=0.5, color='grey')
        
        # Find Noptimal using fminbnd equivalent
        def funct(N):
            return OCNerr(N, g, beta, tau + 0.1)
        
        result = minimize_scalar(funct, bounds=(1, 50), method='bounded')
        Noptimal = result.x
        Noptimal = max(round(Noptimal), 2)

        if round(Noptimal) + sp > len(temps):
            Noptimal = len(temps) - sp
            computing = 0  # extended too far

        if Noptimal <= Nguess1:
            Nguess = round((Noptimal + Nguess1) / 2)
        else:
            Nguess = Noptimal
            ep = sp + Nguess
            continue
        
        # Compute stats and plot for sp through sp+Noptimal-1
        ep = sp + round(Noptimal) 
        if ep <= len(temps):
            avg_now = np.mean(temps[sp:ep])
            std_now = np.std(temps[sp:ep])
        else:
            real_data = len(temps) - sp
            proj_data = ep - len(temps)
            temps2 = np.concatenate([temps[sp:], p[0] * (np.arange(dates[-1]+1, dates[-1] + proj_data+1)) + p[1]])
            avg_now = np.mean(temps2[sp:ep])
            std_now = np.std(temps2[sp:ep])
        
        # Add rectangle and line plot
        if(plotting):
            plt.gca().add_patch(plt.Rectangle((sp + et50, avg_now - std_now / np.sqrt(Noptimal)),
                                              Noptimal, 2 * std_now / np.sqrt(Noptimal),
                                              linewidth=0, facecolor='blue', alpha=0.3))
            plt.plot([sp + et50, ep + 1 + et50], [avg_now, avg_now], linewidth=1, color='green')
            plt.errorbar(sp + Noptimal / 2 + et50, avg_now, yerr=2 * std_now, linewidth=2, color='red', capsize=5)
        
        # Store results in OCNM
        meansf[sp:ep] = avg_now
        sesf[sp:ep] = std_now / np.sqrt(Noptimal)
        for i in range(sp,ep):
            
            if (i-sp)<=1:
                ses[i]= std_guess / np.sqrt(2)
                means[i] = temps[i]
            else:
                ses[i] = np.std(temps[sp:i])/np.sqrt(i-sp -1) #sample std
                means[i] = np.mean(temps[sp:i])
        
        
        # Reset interval
        sp = sp + round(Noptimal)
        if sp + 1 >= len(temps):
            computing = 0
            break
        
        if sp + Nguess > len(temps):
            ep = len(temps)
            computing = 1
        else:
            ep = sp + Nguess

    # Display plot
   # plt.show()

    
    return means, ses, meansf, sesf

