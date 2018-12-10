# VL_spatial_model

This repository houses code developed for the analysis presented in the scientific publication "Spatio-temporal modelling of Leishmania infantum infection among domestic dogs: a simulation study and sensitivity analysis applied to rural Brazil" by Elizabeth Buckingham-Jeffery, Edward M. Hill, Samik Datta, Erin Dilger, and Orin Courtenay.

Note, this code will not run "out of the box" as it additionally requires datasets which are available from the authors on reasonable request and for use in the context of the research study. The corresponding author can be contacted at e.buckingham-jeffery@manchester.ac.uk.

Please find below an explainer of the directory structure within this repository.

## FlyDistributionAnalysis
Smoothing seasonal sandfly data (figures S2, S3, and S4).

## HostDistributionAnalysis
Fitting Poisson and negative binomial distributions to household host data (figures 2 and S1).
Computation of AIC value for each fitted distribution (table S1).

## ModelSimnCode 
Files to carry out spatial simulations of VL infection dynamics. The call structure of the simulation files is outlined below:

### First level

- VL_model_SimnCallScript.m
	-- Parameter values specified. Passed as inputs to VL_model_function.m

### Second level

- VL_model_function.m
	-- Accept as inputs the biological parameters, settlement configuration, simulation variables.

### Third level 

The following files are called within VL_model_function.m

- assign_sandfly_abundance.m
	--  Allocate sandfly population to each household 
	--  ASSUMPTION: These are the baseline sandfly population at each household. Scaled later by seasonality.

- host_pref_allocation.m
	-- Computes preference towards each host type (household specific)

- hostDist.m
	-- Sample from the distribution of the specific host in the specified setting 	
 
- IP_setup.m
	-- Compute infectious pressure surface for current timestep

- run_epidemic.m
	-- Update disease status variables on current timestep

- set_household_locations.m
	-- Specify if households locations should be data informed or use synthetically generated setup

## SensitivityAnalysis
Rank sensitivity of biological parameters (figures 6 and 7) based on method described in: Damiani et al (2013)
Parameter sensitivity analysis of stochastic models: Application to catalytic reaction networks
*Computational biology and chemistry* **42**: 5-17.

## SimnOutputFiles
Simuation outputs that were used to assess biological parameter sensitivity (figure 5).
