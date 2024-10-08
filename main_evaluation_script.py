import os
import importlib
import numpy as np
import pandas as pd

def run_methods(years, avg_temperatures, model_run, experiment_type, methods_folder='Methods'):
    # Find all *_method.py files in the folder and subfolders
    method_files = []
    for root, _, files in os.walk(methods_folder):
        for file in files:
            if file.endswith('_method.py'):
                method_files.append(os.path.join(root, file))

    # Initialize result storage
    results = {}

    # Send data to each method and retrieve results
    for method_path in method_files:
        # Extract method class (folder structure relative to methods_folder)
        method_class = os.path.relpath(os.path.dirname(method_path), methods_folder)

        # Dynamically import the module
        #print(method_path)
        module_name = method_path.replace('/', '.').replace('.py', '')
        #print(module_name)
        method_module = importlib.import_module(module_name, package = "Thorne_15_codefigurestats")

        # Call the method function (assumed to be named "run_method" in each *_method.py)
        result = method_module.run_method(years, avg_temperatures, model_run, experiment_type)

        # Store result along with method class (folder structure)
        method_name = os.path.basename(method_path).replace('_method.py', '')
        results[method_name] = {
            'method_class': method_class,
            'LT_trend': result
        }
    return results

if __name__ == '__main__':
    # First evaluation
    data = pd.read_csv("./Common_Data/HadCRUT5.csv")
    temps_obs = data.loc[:,"Anomaly"].to_numpy()
    years=data.loc[:,"Time"].to_numpy()
    nyrs = len(years)
    model_run = 'hadcrut5'
    experiment_type = 'historical_obs_20yr_mean'

    # Run the methods
    results = run_methods(years, temps_obs, model_run, experiment_type)

    # Example of handling results
    for method_name, method_data in results.items():
        print(f"Results from {method_name} (Method Class: {method_data['method_class']}):")
        result = method_data['LT_trend']
        if isinstance(result, tuple):
            # Central estimate and SE
            central_estimate, se = result
            print(f"  Central Estimate: {central_estimate}")
            print(f"  Standard Error: {se}")
        else:
            # 100 percentiles
            percentiles = result
            print(f"  Percentiles: {percentiles}")
