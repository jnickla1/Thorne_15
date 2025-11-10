# Code written by Bruce T. T. Calvert - July 2025
# Obtains a 100 representative member subset of a 1000 member radiative forcing timeseries ensemble using balanced k-means algorithm

import numpy as np
from sklearn.cluster import KMeans
from scipy.optimize import linear_sum_assignment
from pathlib import Path
import netCDF4
import pandas as pd
from openpyxl import Workbook

def balanced_kmeans(ensemble,number_of_clusters,weights=["null"]):
	ensemble = np.array(ensemble)
	if round(number_of_clusters,0) != number_of_clusters or number_of_clusters < 1:
		raise TypeError("number of clusters must be a positive integer")
	if ensemble.ndim == 1:
		ensemble = ensemble.reshape((-1,1))
	if ensemble.ndim != 2:
		raise TypeError("ensemble must be a 2 dimensional array")
	if ensemble.shape[1] % number_of_clusters != 0:
		raise TypeError("number of ensemble members must be divisible by number of clusters")
	if weights[0] != "null":
		weights = np.array(weights)
		weights = weights.reshape(weights.size)
	if weights[0] == "null":
		weights = np.ones(ensemble.shape[0])
	if weights.shape[0] != ensemble.shape[0]:
		raise TypeError("weights dimensions must agree with number of time periods of hierarchy")
	if weights.shape[0] != weights.size:
		raise TypeError("weights must be a 1 dimensional vector")
	cluster_size = int(ensemble.shape[1]/number_of_clusters)
	ensemble = np.transpose(ensemble)
	# Assign weights
	ensemble_mean = np.mean(ensemble,axis=0)
	ensemble = np.multiply(ensemble-ensemble_mean,np.sqrt(weights)) + ensemble_mean
	# Initialize clusters using K means algorithm
	kmeans = KMeans(n_clusters=number_of_clusters,random_state=np.random.randint((1<<31)-1)).fit(ensemble)
	cluster_centers = kmeans.cluster_centers_
	# Iteratively assign equal numbers of ensemble members to each cluster using Hungarian-like algorithm and recalculate cluster centers
	for iterations in range(1, 100):
		distance_matrix = np.zeros((ensemble.shape[0],number_of_clusters))
		for i in range(ensemble.shape[0]):
			for j in range(number_of_clusters):
				distance_matrix[i,j] = np.linalg.norm(ensemble[i]-cluster_centers[j])
		row_index, column_index = linear_sum_assignment(np.kron(np.ones((1,int(cluster_size))),distance_matrix))
		assigned_cluster = column_index % number_of_clusters;
		for i in range(number_of_clusters):
			cluster_centers[i] = ensemble[assigned_cluster==i].mean(axis=0)
	return assigned_cluster

if __name__ == '__main__':
	# Set Seed of Random Number Generator
	np.random.default_rng(seed=0)
	# download ERF data from https://zenodo.org/records/15630666 and put it in the same folder as this py file before running the code
	variable_names = ['total','co2','ch4','n2o','halogen','o3','aerosol-radiation_interactions','aerosol-cloud_interactions','contrails','land_use','bc_snow','h2o_strat','solar','volcanic']
	wb = Workbook()
	(wb.active).title = 'total'
	wb.save('representative_ERF_ensemble_1750_to_present.xlsx')
	writer = pd.ExcelWriter('representative_ERF_ensemble_1750_to_present.xlsx',engine='openpyxl',mode='a',if_sheet_exists='replace')
	for i in range(len(variable_names)):
		if i == 0:
			data_file = netCDF4.Dataset(Path(__file__).with_name('ERF_DAMIP_1000.nc'))
		else:
			data_file = netCDF4.Dataset(Path(__file__).with_name('ERF_DAMIP_1000_full.nc'))
		data = np.ma.getdata(data_file.variables[variable_names[i]]).data
		if i == 0:
			clusters = balanced_kmeans(data,100)
			selected_members = np.array([])
			for j in range(100):
				selected_members = np.append(selected_members,np.arange(1000)[clusters==j][np.random.permutation(np.arange(sum(clusters==j)))[0]])
		output_ensemble = np.ma.getdata(data_file.variables['time']).data.reshape(-1,1)
		# Obtain one ensemble from each cluster
		for j in range(100):
			output_ensemble = np.append(output_ensemble,data[:,int(selected_members[j])].reshape((-1,1)),axis=1)
		with pd.ExcelWriter('representative_ERF_ensemble_1750_to_present.xlsx',engine='openpyxl',mode='a',if_sheet_exists='replace') as writer:
			df = pd.DataFrame(output_ensemble,columns=['time']+['ensemble'+str(j) for j in range(1, 101)])
			df.to_excel(writer,sheet_name=str(variable_names[i]),index=False)
