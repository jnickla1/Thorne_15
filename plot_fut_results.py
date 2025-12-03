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


#sel_methods = ['CGWL10y_sfUKCP', 'EBMKF_ta4','FaIR_comb_unB', 'FaIR_nonat_unB','GAM_AR1','GWI_anthro_SR15','Kal_flexLin','lowess1dt36wnc']
#already sorted
#sel_methods_colors = [ "#CC3311", "#999933","#44AA99" ,"#44AA99","#117733","#AA4499","#88CCEE","#DDCC77"]
sel_methods = ['CGWL10y_sfUKCP', 'EBMKF_ta4','FaIR_comb_unB', 'GAM_AR1','GWI_tot_CGWL','Kal_flexLin','lowess1dt36wnc']
sel_methods_colors = [ "#CC3311", "#999933","#44AA99" ,"#117733","#AA4499","#88CCEE","#DDCC77"]
sel_methodsA = ['FaIR_nonat_unB','GWI_anthro_CGWL']
sel_methods_colorsA = ["#44AA99","#AA4499"]
# Define the file pattern and loop through all matching files
file_pattern = "/Results/current_fut_statistics_fut_ESM1-2-LR_SSP{ssp_id}_constVolc{run}.csv"
rmse_cols = ['RMS', '100RMS','75RMS', 'e100RMS','e75RMS']

if __name__ == '__main__':

    for ssp_id in [126, 245, 370]:
        df_list = []
        for run in range(50):
            file_name = os.getcwd() + file_pattern.format(ssp_id=ssp_id, run=run)
            if os.path.exists(file_name):  # Check if the file exists
                # Read the CSV file
                df = pd.read_csv(file_name)
                
                # Filter rows with the desired method_name
                filtered = df[df['method_name'].isin(sel_methods)]
                df_sorted = filtered.sort_values(by='method_name')
                filteredA = df[df['method_name'].isin(sel_methodsA)]
                df_sortedA = filteredA.sort_values(by='method_name')
                
                df_list.append(df)
                #print(run)

                # Extract Edyrs15 values and add to the corresponding list
                #breakpoint()
                if ssp_id == 126:
                    ssp_126.append(df_sorted['Fdyrs15'].values)
                    ssp_126e.append(df_sortedA['eFdyrs15'].values)
                    #if (abs(df_sorted.loc[df_sorted['method_name'] == 'EBMKF_ta2', 'Fdyrs15'].values[0]) > 10):
                    #    print(file_name)
                elif ssp_id == 245:
                    ssp_245.append(df_sorted['Fdyrs15'].values)
                    ssp_245e.append(df_sortedA['eFdyrs15'].values)
                elif ssp_id == 370:
                    ssp_370.append(df_sorted['Fdyrs15'].values)
                    ssp_370e.append(df_sortedA['eFdyrs15'].values)
        if df_list:
            # Concatenate all filtered/sorted dataframes
            all_data = pd.concat(df_list)
            # Compute numeric means excluding RMSE columns
            numeric_cols = all_data.select_dtypes(include=[np.number]).columns.difference(rmse_cols)
            averaged_numeric = all_data.groupby('method_name')[numeric_cols].mean().reset_index()
            # Compute RMSE-style means for the specified columns
            rmse_processed = (
                all_data
                .groupby('method_name')[rmse_cols]
                .apply(lambda x: np.sqrt((x**2).mean()))
                .reset_index())
            # Extract constant text columns
            text_columns = all_data[['method_name', 'method_class', 'c/r']].drop_duplicates('method_name')
            # Merge all parts together
            final_result = (
                text_columns
                .merge(averaged_numeric, on='method_name', how='inner')
                .merge(rmse_processed, on='method_name', how='inner')
            )
            # Save to CSV
            final_result.to_csv(f"averaged_runs{ssp_id}.csv", index=False)
    try:
        np.array(ssp_126)
    except:
        breakpoint()

    file_pattern2 = "/Results/current_fut_statistics_fut_NorESM_RCP45_{volc_id}{run}.csv"

    v_45 = []
    v_45e = []
    cv_45 = []
    cv_45e = []

    for volc_id in ["Volc","VolcConst"]:
        df_list=[]
        for run in range(60):
            file_name = os.getcwd() + file_pattern2.format(volc_id=volc_id, run=run)
            if os.path.exists(file_name):  # Check if the file exists
                # Read the CSV file
                df = pd.read_csv(file_name)
                # Filter rows with the desired method_name
                filtered = df[df['method_name'].isin(sel_methods)]
                df_sorted = filtered.sort_values(by='method_name')
                filteredA = df[df['method_name'].isin(sel_methodsA)]
                df_sortedA = filteredA.sort_values(by='method_name')
                
                df_list.append(df)
                # Extract Edyrs15 values and add to the corresponding list
                if volc_id == "Volc":
                    v_45.append(df_sorted['Fdyrs15'].values)
                    v_45e.append(df_sortedA['eFdyrs15'].values)
                    
                elif volc_id == "VolcConst":
                    cv_45.append(df_sorted['Fdyrs15'].values)
                    cv_45e.append(df_sortedA['eFdyrs15'].values)
                    #if (abs(df_sorted.loc[df_sorted['method_name'] == 'EBMKF_ta2', 'Fdyrs15'].values[0]) > 10):
                    #    print(file_name)
        if df_list:
            # Concatenate all filtered/sorted dataframes
            all_data = pd.concat(df_list)
            # Compute numeric means excluding RMSE columns
            numeric_cols = all_data.select_dtypes(include=[np.number]).columns.difference(rmse_cols)
            averaged_numeric = all_data.groupby('method_name')[numeric_cols].mean().reset_index()
            # Compute RMSE-style means for the specified columns
            rmse_processed = (
                all_data
                .groupby('method_name')[rmse_cols]
                .apply(lambda x: np.sqrt((x**2).mean()))
                .reset_index())
            # Extract constant text columns
            text_columns = all_data[['method_name', 'method_class', 'c/r']].drop_duplicates('method_name')
            # Merge all parts together
            final_result = (
                text_columns
                .merge(averaged_numeric, on='method_name', how='inner')
                .merge(rmse_processed, on='method_name', how='inner')
            )
            # Save to CSV
            final_result.to_csv(f"averaged_runs{volc_id}.csv", index=False)

    # Create a figure with 3 subplots for the smoothed histograms
    fig, axs = plt.subplots(5, 2, figsize=(10, 10), sharex=True,gridspec_kw={'height_ratios': [1,1,1,1,1]})
    histstack=True
    # Define a helper function to plot soothed histograms
    def plot_smoothed_histogram(ax, data0, title , leg=False, switch = 0,anthro=False):
        
        if( not histstack):
            for i in range(np.shape(data0)[1]):
                data = data0[:,i] #don't need to flip it
                kde = gaussian_kde(data)
                x_range = np.linspace(min(data)-1.5, max(data)+1.5, 1000)

                
                #if(i//4 !=switch): switch to plot only the label in half of them
                #    templabel=None
                
                if(not anthro):
                    templabel = sel_methods[i]
                    ax.plot(x_range, kde(x_range), color=sel_methods_colors[i], label = templabel, lw=2)

                else:
                    templabel = sel_methodsA[i]
                    ax.plot(x_range, kde(x_range), color=sel_methods_colorsA[i], label = templabel, lw=2)
                #else:
                #    ax.plot(x_range, kde(x_range), color=sel_methods_colors[i], label = templabel, lw=2,linestyle="--")
                #ax.hist(data,density=True,color=sel_methods_colors[i],alpha=0.3)
                #ax.fill_between(x_range, kde(x_range), color=sel_methods_colors[i], alpha=0.2)

        else:
            bins = np.arange(int(np.min(data0)-1.5), int(np.max(data0)+1.5), 1)
            if (not anthro):
                ax.hist(data0,bins=bins,density=True,stacked=True,
                    color=sel_methods_colors,alpha=0.5,edgecolor="none",label=sel_methods)
            else:
                ax.hist(data0,bins=bins,density=True,stacked=True,
                    color=sel_methods_colorsA,alpha=0.5,edgecolor="none",label=sel_methodsA)
        
        ax.set_title(title, fontsize=14)
        ax.grid(True)
        ax.set_ylim([0,0.35])
        ax.tick_params(axis='x', labelbottom=True)
        if (leg and anthro):
            ax.legend(loc='best',fontsize=8.3)
        elif (leg and not anthro):
            ax.legend(loc='upper right',fontsize=8.3,ncol=3,bbox_to_anchor=(1.1, 1))

    # Plot smoothed histograms for each SSP
    plot_smoothed_histogram(axs[2,0], np.array(ssp_126), "MPI ESM1-2-LR SSP126",leg=True )
    plot_smoothed_histogram(axs[1,0], np.array(ssp_245), "MPI ESM1-2-LR SSP245", )
    plot_smoothed_histogram(axs[0,0], np.array(ssp_370), "Realized warming:\n comparison within one member\nMPI ESM1-2-LR SSP370")
    plot_smoothed_histogram(axs[3,0], np.array(cv_45), "NorESM RCP45 VolcConst",)
    plot_smoothed_histogram(axs[4,0], np.array(v_45), "NorESM RCP45 Volc")

    plot_smoothed_histogram(axs[2,1], np.array(ssp_126e), "MPI ESM1-2-LR SSP126",leg=True,switch=1 ,anthro=True) 
    plot_smoothed_histogram(axs[1,1], np.array(ssp_245e), "MPI ESM1-2-LR SSP245", anthro=True)
    plot_smoothed_histogram(axs[0,1], np.array(ssp_370e), "Attributable warming:\n comparison to ensemble average\nMPI ESM1-2-LR SSP370",anthro=True)
    plot_smoothed_histogram(axs[3,1], np.array(cv_45e), "NorESM RCP45 VolcConst",anthro=True)
    plot_smoothed_histogram(axs[4,1], np.array(v_45e), "NorESM RCP45 Volc",anthro=True)

    axs[4,0].set_xlim(-10,10)
    # Set x-axis label
    axs[4,0].set_xlabel("Error in Crossing Instant (years)", fontsize=12)
    axs[4,1].set_xlabel("Error in Crossing Instant (years)", fontsize=12)


    cai = [0,0]
    caiht = np.array([.2,.2])
    wid_txt = 2.5

    for i,ax in enumerate([axs[1,0],axs[1,1]]):
        ax.arrow(cai[i]-2*wid_txt, caiht[i], -wid_txt, 0,head_width=caiht[i]*3/12,head_length=wid_txt/2,facecolor='white',lw=2)
        ax.text(cai[i]-2.2*wid_txt, caiht[i]*.75,"Too Early", fontsize=12,horizontalalignment='center',color='black')
        ax.arrow(cai[i]+2*wid_txt, caiht[i], wid_txt, 0,head_width=caiht[i]*3/12,head_length=wid_txt/2,facecolor='white',lw=2)
        ax.text(cai[i]+2.2*wid_txt, caiht[i]*.75,"Too Late", fontsize=12,horizontalalignment='center',color='black')


    enscross_time = [2033.125, 2035.458333,  2042.125, 2033.791666, 2034.208333 ]
    for i in range(5):
        ax = axs[i,1]
        ax.text(0, 0.3,'{0:.1f}'.format(enscross_time[i]), fontsize=12,horizontalalignment='center',color='black',fontweight='bold')

    fig.suptitle("Time of 1.5°C Threshold Crossing: \n(Method's) First Reported — 20-yr Mean of ...", fontsize=16)
    # Adjust layout and show the plot
    plt.tight_layout()

    plt.savefig("combined7histogram"+str(histstack)+".png",dpi=500,bbox_inches='tight')
#plt.show()
