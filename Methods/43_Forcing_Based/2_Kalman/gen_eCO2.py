import os
import numpy as np
import pandas as pd
import xarray as xr
#from netCDF4 import Dataset
#import argparse
import pdb


def average_every_n(lst, n):
    """Calculates the average of every n elements in a list."""
    return np.array([np.mean(lst[i:i + n]) for i in range(0, len(lst), n)])



def read_co2_values(file_path):
    """Reads modeled CO2 values in kg from a file."""
    #dataset = Dataset(filepath, 'r')
    #variable = dataset.variables['co2mass']
    with xr.open_dataset(file_path) as ds:
        vararray = ds['co2mass'].values
    return average_every_n(vararray[0,:], 12)

def convert_co2_to_ppm(co2_values):
    """Converts CO2 values from kg to ppm using the given conversion factor."""
    conversion_factor =  4.27e-10 / 3341
    return co2_values * conversion_factor

def read_n2o_fraction(ssp, n2o_dir):
    """Reads the N2O fraction from a matching file in the directory."""
    files = os.listdir(n2o_dir)
    split_files = [f.split("-") for f in files]
    n2o_file = False
    for inner_array in split_files:
        # Ensure the inner array has enough elements to check indices
        if len(inner_array) > 8 :
            if (inner_array[8] == ssp or inner_array[9] == ssp)  and inner_array[3] == "nitrous":
                n2o_file = "-".join(inner_array)
                break
    
    if not n2o_file:
        raise FileNotFoundError(f"No matching N2O file found for SSP {ssp}.")
    n2o_path = os.path.join(n2o_dir, n2o_file)
    with xr.open_dataset(n2o_path) as ds:
        n2o_fraction = ds["mole_fraction_of_nitrous_oxide_in_air"].values

    #1750-2500, we want 1850-2100
    return n2o_fraction[100:(2100-1750+1),0]  #0th is global average

def calculate_relative_forcing(co2_ppm, n2o_ppb):
    """
    Calculate the relative forcing (RF_CO2) using the Oslo line-by-line model equations.
    https://doi.org/10.5194/gmd-13-3571-2020, table 3

    Parameters:
    co2_ppm (float): CO2 concentration in ppm.
    n2o_ppb (float): N2O concentration in ppb.

    Returns:
    float: The relative forcing RF_CO2 in W/m2
    """
    # Constants
    a1 = -2.4785e-7  # W/m2 ppm-2
    b1 = 0.00075906  # W/m2 ppm-1
    c1 = -0.0021492  # W/m2 ppb0.5
    d1 = 5.2488      # W/m2
    C0 = 277.15      # ppm
    Ca_max = C0 - (b1 / (2 * a1))  # Approx. 1808 ppm
    
    # Calculate α'
    if co2_ppm > Ca_max:
        alpha_prime = d1 - (b1**2 / (4 * a1))
    elif C0 < co2_ppm <= Ca_max:
        alpha_prime = d1 + (a1 * (co2_ppm - C0)**2) + (b1 * (co2_ppm - C0))
    else:
        alpha_prime = d1
    
    # Calculate α_N2O
    alpha_n2o = c1 * np.sqrt(n2o_ppb)
    
    # Calculate RF_CO2
    rf_co2 = (alpha_prime + alpha_n2o) * np.log(co2_ppm / C0)
    
    return rf_co2

def read_erf_data(erf_file):
    """Reads anthropogenic GHG forcing from the ERF file."""
    erf_data = pd.read_csv(erf_file)
    ghg_columns = [ #"co2",
         "ch4", "n2o","other_wmghg", "o3","h2o_stratospheric","contrails",
        "land_use","bc_on_snow"
    ]
    erf_data["total_anthro_ghg_forcing_mco2"] = erf_data[ghg_columns].sum(axis=1)
    return erf_data["total_anthro_ghg_forcing_mco2"].values

def recompute_volc_erf_data(erf_file):
    """Reads natural GHG forcing from the ERF file."""
    erf_data = pd.read_csv(erf_file)
    volc = erf_data["volcanic"].values[100:(2100-1750+1)]
    volc[(2014-1850):] = np.mean(volc[:(2014-1850)])
    return volc, erf_data["solar"].values[100:(2100-1750+1)]

def calculate_equivalent_co2(forcing_values):
    """Converts forcing values into equivalent CO2 concentration in ppm."""
    equivalent_co2 = 10**((forcing_values + 31.149 )/ 12.7433 )
    return equivalent_co2

    

def recompute_eCO2(ssp, cmip_model):
    # File paths
    co2_file = os.path.expanduser(f"~/climate_data/{cmip_model}/co2mass/{ssp}_co2mass.nc")
    n2o_dir = os.path.expanduser("~/climate_data/SSP_inputdata/UoM1_2_0")
    erf_file = os.path.expanduser(f"~/climate_data/SSP_inputdata/ERFs-Smith-ar6/ERF_{ssp}_1750-2500.csv")
    output_dir = os.path.expanduser("~/climate_data/SSP_inputdata")

    # Step 1: Read CO2 values and convert to ppm
    co2_hist = os.path.expanduser(f"~/climate_data/{cmip_model}/co2mass/historical_co2mass.nc")
    co2_values_hist = read_co2_values(co2_hist)
    
    co2_values = read_co2_values(co2_file)
    co2_ppm = convert_co2_to_ppm(np.concatenate((co2_values_hist,co2_values)))

    # Step 2: Read N2O fraction
    n2o_fraction = read_n2o_fraction(ssp, n2o_dir)

    # Step 3: Calculate relative forcing
    erf_co2 = np.array([calculate_relative_forcing(co2,n2o) for co2,n2o in zip(co2_ppm, n2o_fraction)])

    # Step 4: Read ERF data and calculate total forcing
    total_forcing = read_erf_data(erf_file)

    # Step 5: Calculate equivalent CO2
    equivalent_co2 = calculate_equivalent_co2(total_forcing[100:(2100-1750+1)] + erf_co2)
    print(equivalent_co2)
    #breakpoint()
    # Step 6: Save to CSV
    years = np.arange(1850, 2101)  # Generate years from 1850 to 2100
    df = pd.DataFrame({"year": years,"eCO2": equivalent_co2})

    # Write the DataFrame to a CSV file
    output_path = os.path.join(output_dir, f"eCO2_{cmip_model}_{ssp}.csv")
    #df.to_csv(output_path, index=False)
    print(f"Saved equivalent CO2 avalues to {output_path}")

    volc,solar = recompute_volc_erf_data(erf_file)
    # Write the DataFrame to a CSV file
    output_path2 = os.path.join(output_dir, f"ERF_{cmip_model}_{ssp}.csv")
    df2 = pd.DataFrame({"year": years,"CO2": co2_ppm, "ERF_CO2": erf_co2, "ERF_anthro": total_forcing[100:(2100-1750+1)] + erf_co2,
                        "ERF_volc": volc, "ERF_natural":volc+solar})
    df2.to_csv(output_path2, index=False)

if __name__ == "__main__":
    #parser = argparse.ArgumentParser(description="Generate equivalent CO2 values for a given SSP and CMIP model.")
    #parser.add_argument("ssp", type=str, help="The SSP scenario (e.g., ssp370).")
    #parser.add_argument("cmip_model", type=str, help="The CMIP model (e.g., ESM1-2-LR).")
    #args = parser.parse_args()

    recompute_eCO2("ssp245","ESM1-2-LR")
