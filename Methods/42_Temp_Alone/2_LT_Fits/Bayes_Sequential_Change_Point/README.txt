%************************************************************************/
%* The Sequential Bayesian Change Point algorithm - A program to        */
%* caluclate the posterior probability of a change point in a time      */
%* series.                                                              */
%*                                                                      */
%* Please acknowledge the program author on any publication of          */
%* scientific results based in part on use of the program and           */
%* cite the following article in which the program was described.       */
%*                                                                      */
%* E. Ruggieri and M. Antonellis.  "An exact approach to Bayesian       */
%* sequential change point detection"                                   */
%*                                                                      */
%* Program Author: Eric Ruggieri                                        */
%* College of the Holy Cross                                            */
%* Worcester, MA 01610                                                  */
%* Email:  eruggier@holycross.edu                                       */
%*                                                                      */
%* Copyright (C) 2016  College of the Holy Cross                        */
%*									*/
%* Acknowledgements: This work was supported by a grant from the        */
%* National Science Foundation, DMS-1407670 (E. Ruggieri, PI)           */
%*                                                                      */
%* The Sequential Bayesian Change Point algorithn is free software:     */
%* you can redistribute it and/or modify it under the terms of the GNU  */
%* General Public License as published by the Free Software Foundation, */
%* either version 3 of the License, or (at your option) any later       */
%* version.                                                             */
%*                                                                      */
%* The Sequential Bayesian Change Point algorithm is distributed in the */
%* hope that it will be useful, but WITHOUT ANY WARRANTY; without even  */ 
%* the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR  */
%* PURPOSE.  See the GNU General Public License for more details.       */
%*                                                                      */
%* You should have received a copy of the GNU General Public License    */
%* along with the Sequential Bayesian Change Point algorithm.  If not,  */
%* see <http://www.gnu.org/licenses/> or write to the                   */
%* Free Software Foundation, Inc.                                       */
%* 51 Franklin Street Fifth Floor                                       */
%* Boston, MA 02110-1301, USA.                                          */
%************************************************************************/


Main Script: Bayes_changepoint_sequential.m
By default, the algirhtm is coded to run the Global Surface Temperature Anomaly 
time series [Temperature_Anomaly2013.txt] analyzed in the manuscript. To run a different 
data set, modifications are needed only where the data is loaded and where 
figures are created (plot_results.m).

 Outline of the Sequential Bayesian Change Point algorithm:
 1) Load the data  
 2) Define the parameter values
 3) For the initial data set:
   a) Calculate the probability density of the data for each sub-interval
   b) Forward Recursion [Dynamic Programming]
   c) If a change point is detected:
       Stochastic Backtrace via Bayes rule
       i) Sample a number of change points
       ii) Sample the location of the change points
       iii) Sample the regression parameters between adjacent change points
       iv) Plot the results
 4) (As needed) Add additional data points one at a time to existing data set
 5) Repeat steps 3a-3c for new observation
 6) Plot the results
________________________________________________________________________
Data Sets:
Temperature_Anomaly2013.txt - Global surface temperature anomalies from 1880-2013.  
	Source: http://www.ncdc.noaa.gov/monitoring-references/faq/anomalies.php
well_log_data.txt - The original well-log data set of:
	O’Ruanaidh, J. and Fitzgerald, W.J. (1996) Numerical Bayesian methods applied to signal processing. 
	New York: Springer-Verlag.
	Source: http://onlinelibrary.wiley.com/journal/10.1111/(ISSN)1467-9868/homepage/65_4.htm
well_log_data_no_outliers.txt - The well-log data set with outliers manually removed.  
	This was the actual data set used in the analysis of Section 5.3

_______________________________________________________________________
Associated functions (See below for more complete definition of variables):
addone.m - Function performs steps 3a and 3b above for each new observation added 
	to the data set.  Specifically, function calculates the probability of any
	segment of the data which ends with the new observation (Equation 2 from 
	Ruggieri 2013) and then performs the Forward Recursion step (Equation 2) 
	of the Sequential Bayesian Change Point algorithm.
	Input: User specified parameters (see below), data set (Y), model (X),
		current versions of matrices Py and P
	Output: Updated versions of matrices Py and P
addone_well_log.m - The same as the above function but contains additional lines of 
	code that address the potential underflow issues with very large data sets.
	Input: User specified parameters (see below), data set (Y), model (X),
		current versions of matrices Py and P
	Output: Updated versions of matrices Py and P
Bayes_changepoint_sequential_well_log.m - The same as the main script 
	(Bayes_changepoint_sequential.m), except that the well-log data set is used 
	instead of temperature anomalies.  In addition,	code exists to remove outliers 
	from the data set and the appropriate well-log function calls are used instead 
	of the basic function calls.  The difference is due to potential underflow 
	issues associated with long data sets.
	Input: None
	Output: None
find_chgpts_simulation.m - For the multiple change point simulation of Section 5.2, 
	this function draws samples from the exact posterior distribution of the
	number and location of change points for each new change point that has been
	detected.  Function returns a vector chgpt_loc that contains the sampled 
	change point positions.
	Input: User specified parameters (see below), matrices Py and P, vector k
		containing posterior distribution of number of change points
	Output: Vector chgpt_loc containing the sampled change point locations
nCk.m - A replacement for Matlab's nchoosek function when 'N' is large
	Input: Values for N and k
	Output: The value of the combination N choose K
partition_fn.m - Performs the Forward Recursion step of the Sequential Bayesian 
	Change Point algorithm [Equation (2)].
	Input: Matrix Py, kmax (maximum number of change points), N (length of data set)
	Output: Matrix P
pick_k1.m - General function that samples from a discrete probability distribution. 
	Input: Vector containing a discrete probability distribution
	Output: A random sample from that vector.
plot_results.m - General function that draws samples of the number of change points,
	their locations, and the parameters of the regression model from the exact
	posterior distribution of these quantities.  Axis labels for the resulting
	figure are specific to the temperature anomalies data set, but can easily 
	be modified to suit a particular need.
	Input: User specified parameters (see below), data set (Y), model (X),
		time points (x), vector k (posterior distribution of # change points)
		current versions of matrices Py and P
	Output: Plot of the data set, inferred model, and posterior probability of 
		a change point at each location
plot_results_well_log.m - Same as the plot_results function above, except that it 
	contains additional lines of code to address underflow in large data sets.
	Input: User specified parameters (see below), data set (Y), model (X),
		time points (x), vector k (posterior distribution of # change points)
		current versions of matrices Py and P
	Output: Plot of the data set, inferred model, and posterior probability of 
		a change point at each location
process_chgpts.m - For the multiple change point simulation of Section 5.2, after 
	the change point locations have been sampled from the posterior distribution
	using the function find_chgpts_simulation, this function takes the output
	vector chgpt_loc and processes the vector to obtain detection speeds, 
	detection rates, bias, and variance of the location of each of the three 
	change points for each of the three detection criteria (mean, median, and
	mode of the posterior distribution of the number of change points)
	Input: vector chgpt_loc, number of sampled solutions, actual location of 
		change points, current versions of matrices locations, variances,
		found, and detect.
	Output: Updated versions of matrices locations (gives bias of change point
		locations), variances (standard deviation of poserior distribution), 
		found (detection rate), and detect (detection speed)
Sequential_simulation.m - Main script that creates and processes the simulated 
	data sets as described in Section 5.1. Specifically, the script generates
	a constant model containing no change points and a constant linear trend 
	containing no change points.
	Input: None
	Output: None
Sequential_Simulation_chgpts.m - Script to generate the multiple change point 
	simulation of Section 5.2.  Calls the functions find_chgpts_simulation.m
	and process_chgpts.m
	Input: None
	Output: None

________________________________________________________________________
Description of Variables:


User Specified Variables:
Y_complete = The complete data set / time series		N x 1 vector
x_complete = The complete set of time points, which may or 
		may not be equally spaced			N x 1 vector

Y = A subset of the data set / time series.			N x 1 vector
	Subset so that you can add additional observations 
	one at a time (sequentially)   		 
x = Subset of the time points which may or may not be equally spaced
X = Predictor variables.     					N x m matrix 
    In the case of the global surface temperature anomalies time series, a linear model

N = Number of data points in data set being analyzed
m = Number of regressors

d_min = Minimum distance between adjacent change points
k_0 = Hyperparameter for the prior on the regression coefficients (betas)
v_0, sig_0 = Hyperparameter for scaled inverse chi-square prior on the variance
k_max = Maximum number of change points allowed
num_samp = Number of sampled solutions

parameters = The six parameters listed above in the order given

beta0 = Mean of multivariate normal prior on regression coefficients.		m x 1 vector 
	Note: beta0 is also defined at the top of the addone function.

Variables Generated by Algorithm:
Py = Probability density of the data for every possible substring of the data.  N x N matrix  
     (Equation 2 from Ruggieri 2013).  

P = Probability density of the data containing k change points.			kmax x N matrix
    Equation (2).    
k = Posterior distribtuion on the number of change points.  			k_max x 1 vector
    Equations (4) and (7)
KK = matrix where each column is the vector k at a given point in time		kmax x N matrix
     i.e. You can use KK to view how the posterior distribution on the number
     of change points changes through time.

M = Current number of change points according to the posterior distribution
    using the specified detection criteria
current_num = Number of change points that have already been detected.
		(so if M > current_num, a new change point has been detected)

____________________________________________________________________________
Potential variables of interest within plot_results function:

chgpt_loc = Posterior probability of a change point at each data point		1 x N vector 
samp_holder = Contains the 'num_samp' individual change point solutions		num_samp x k_max matrix
	      Equations (5) and (8)

BETA = The average regression coefficient for each of the predictors at each data point
       (Equation 9 from Ruggieri 2013)						m X N matrix
model = The average model produced by the Sequential Bayesian Change Point algorithm	
										1 x N vector 


____________________________________________________________________________
Potential variables of interest in simulations (not already listed above):

*************** NO CHANGE POINTS (Section 5.1) *************************
num_sim = Number of simulated data sets to generate
num_data = Length of each simulated data set

intercept / sign / trend = Values used to create the constant linear model
			   which contains no change points
NOISE = magnitude of Gaussian white noise added to each simulated data set

M1/M2/M3 = Current number of change points according to the posterior distribution
    	using the specified detection criteria (1 = Mode, 2 = Median, 3 = Mean)
current_num1 / current_num2 / current_num3 = Number of change points that have 
		already been detected by the three detection criteria.
		(so if M > current_num, a new change point has been detected)

count_mean = Total # of falsely detected change points using the mean as detection criteria
count_median = Total # of falsely detected change points using the median as detection criteria
count_mode = Total # of falsely detected change points using the mode as detection criteria


************* MULTIPLE CHANGE POINTS (Section 5.2) ********************
found_mode / found_median / found_mean / found_total 				num_sim x 3 matrices
	Indicates whether or not a particular change point has been detected 
	for each of the three detection criteria in a given simulation.  
	Total represents a batch analysis of the data set

locations_mode / locations_median / locations_mean / locations_total		num_sim x 3 matrices 
	The difference between the inferred and true change point location for 
	each of the three change points and each of the three detection criteria 
	in each simulation. 
locations = Average of the differences recorded in the above matrices			4x3 matrix
stdev_loc = Standard deviation of the differences recorded in the above matrices	4x3 matrix

variances_mode / variances_median / variances_mean / variances_total 		num_sim x 3 matrices
	Standard deviation of the posterior distribution for each of the change
	points identified by the three detection criteria.
variances = Average standard deviation of the posterior distribution for 		4x3 matrix
	each change point	

detect_mode / detect_median / detect_mean 					num_sim x 3 matrices
	Specifies the detection speed for each of the three detection criteria
	for each simulation.  Matrix entry is the number of data points after the
	true change point location at which detection occurred.
detect = Avearge detection speed for each change point and each detection criteria	4x3 matrix

pct_mode / pct_median / pct_mean / pct_total						1 x 3 vectors
	The percentage of times that each of the three change points was detected
	by each of the three detection criteria.  
	Total represents a batch analysis of the data set. 

CHGPTS = Actual change point locations for each simulation
chgpt_loc = Posterior probability of a change point at each data point			1 x N vector
_________________________________________________________________________________

Additional variables of interest for Well-Log data set (not already listed above):

******************** Section 5.3 - Well-Log Data Set **************************
outliers = vector indicating the positions of outliers within the data set.  

