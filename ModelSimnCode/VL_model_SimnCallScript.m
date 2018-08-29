%% LEISHMANIASIS CODE - SCRIPT TO CALL VL_model_function
%SENSITIVITY ANALYSIS OF BIOLOGICAL PARAMETERS

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER VALUES TO BE TESTED     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set up to allow for different number of test values per parameter

BiolParam_to_test = 15; %Number of parameters to test

biol_param = cell(BiolParam_to_test,1);
%Row vectors of values to test for biologic parameters of
%interest.
%Note, first value for each biological parameter is the baseline value
biol_param{1} = [0.3,0.02,0.7,2]; %radius of effect (in km) around house within which an infected dog can infect sandflies
biol_param{2} = [0.55,0.14,0.28,0.42]; %Propn. of infected dogs that are never infectious 
biol_param{3} = [0.37,0.25,0.6,0.8]; %Propn. of dogs that are highly infectious
biol_param{4} = [0.13,0.0064,0.29,0.43]; %Propn. of dogs that are infected when "born" (through immigrations of adult dogs, not puppies) 
biol_param{5} = [182.5,153,212,240]; %Average latent period length (days)
biol_param{6} = [0.046,0.036,0.070,0.095]; %Exposed/Never infectious dog mortality rate (per month)
biol_param{7} = [0.060,0.036,0.080,0.095]; %Low infectious dog mortality rate (per month)
biol_param{8} = [0.064,0.036,0.080,0.095]; %Highly infectious dog mortality rate (per month)
biol_param{9} = [0.038,0.0083,0.020,0.030]; %All other dog class mortality rate (per month)
biol_param{10} = [121,0,243,578]; %Dog replacement time (average length of time, in days, between dog death and new dog being brought into household)
biol_param{11} = [1/3,0.25,0.4,0.5]; %Sandfly bite rate (per day)
biol_param{12} = [0.01,0.002,0.1,0.26]; %Baseline proportion of sandflies infected
biol_param{13} = [0.321,0.1,0.2,0.5]; %Prob. a susceptible dog becomes infected when fed on by an infected sandfly
biol_param{14} = [0.275,0.023,0.15,0.45]; %Propn. of sandflies infected by infectious dogs
biol_param{15} = [0.9,0.75,0.8,0.85]; %Propn. of female sandfly population not observed in CDC traps 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% JOB SUBMISSION DETAILS - BEGIN     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set job number (interger between 1-46 inclusive)
job_num = 1;
fprintf('job index is %d \n', job_num);

%Seed random number generator
new_seed = job_num;
fprintf('new seed is %d \n', new_seed);
rng(new_seed,'twister');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Construct array of parameter index values for all simulations
%Rows - job number being used; Columns - Biological parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param_sets_to_test =  sum(cellfun('length',biol_param)) + 1 - BiolParam_to_test; 
param_idx_array = ones(param_sets_to_test,BiolParam_to_test); 
end_idx = 1;
for i=1:BiolParam_to_test
    NonBaseline_values = numel(biol_param{i}) - 1; %Number of non-baseline values parameter i has 
    start_idx = end_idx + 1; %Update row index to begin edits
    end_idx = end_idx + NonBaseline_values; %Update row to end edits
    param_idx_array(start_idx:end_idx,i) = 2:1:NonBaseline_values+1;
end

%Get param index numbers for current run
param_idxs = param_idx_array(job_num,:);

%%% param_idxs vector - correspond to following biological parameters %%%%%%
%In total, 15 parameters
%param_idxs(1) - radius of effect around house within which an infected dog can infect sandflies
%param_idxs(2) - Propn. of infected dogs that are never infectious 
%param_idxs(3) - Propn. of dogs that are highly infectious
%param_idxs(4) - Propn. of dogs that are infected when "born" (i.e. immigrations of adult dogs, not puppies)
%param_idxs(5) - Average latent period length
%param_idxs(6) - Exposed/Never infectious dog mortality rate (per month)
%param_idxs(7) - Low infectious dog mortality rate (per month)
%param_idxs(8) - Highly infectious dog mortality rate (per month)
%param_idxs(9) - All other dog class mortality rate (per month)
%param_idxs(10) - Dog replacement time (average length of time between dog death and new dog being brought into household)
%param_idxs(11) - Sandfly bite rate (per day)
%param_idxs(12) - Baseline proportion of sandflies infected
%param_idxs(13) - Prob. a susceptible dog becomes infected when fed on by an infected sandfly
%param_idxs(14) - Propn. of sandflies infected by infectious dogs
%param_idxs(15) - Propn. of female sandfly population not observed in CDC traps 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% JOB SUBMISSION DETAILS - END      %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER VALUE ASSIGNMENT - BEGIN %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Construct vector of chosen biological param values
chosen_biol_param = zeros(BiolParam_to_test,1);
for i=1:BiolParam_to_test
    chosen_biol_param(i) = biol_param{i}(param_idxs(i));
end

%Use param_idxs to get value for current run to be used as an input to
%function call at end of this script
dog_rad_of_effect = biol_param{1}(param_idxs(1));
neverinf_dog_propn = biol_param{2}(param_idxs(2));
propn_infectious_dogs_highinf = biol_param{3}(param_idxs(3));
propn_newborn_dog_inf = biol_param{4}(param_idxs(4));
dog_latent_period = biol_param{5}(param_idxs(5));
dog_Exp_NeverInf_MortRate = biol_param{6}(param_idxs(6));
dog_LowInf_MortRate = biol_param{7}(param_idxs(7));
dog_HighInf_MortRate = biol_param{8}(param_idxs(8));
dog_Sus_MortRate = biol_param{9}(param_idxs(9));
dog_Replacement_Time = biol_param{10}(param_idxs(10));
Sandfly_BiteRate = biol_param{11}(param_idxs(11));
Baseline_Sandfly_InfPropn = biol_param{12}(param_idxs(12));
Prob_SusDog_AcquireInf = biol_param{13}(param_idxs(13));
Propn_Sandflies_InfbyInfDogs = biol_param{14}(param_idxs(14));
Propn_flies_NotObvs = biol_param{15}(param_idxs(15));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER VALUE ASSIGNMENT - END   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETTLEMENT TYPE DECLARATION - BEGIN  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify settlement spatial config. and type and household num, grid size
% Entry definitions
% 1: Configuration type (1 - real village data, 2 - random, 3 - regular, 4 - large cluster, 5 - small cluster)
% 2: Settlement type  (0 for rural, 1 for semi urban, 2 for urban)
% 3: household_num
% 4: grid size
settlement_variables = [1,0,200,10];

%If using regular or clustered household configuration, check requested
%file for requested number of households is available
if settlement_variables(1) == 3 || settlement_variables(1) == 4 || settlement_variables(1) == 5
   if  settlement_variables(3) ~= 100 && settlement_variables(3) ~= 200 && settlement_variables(3) ~= 500 && settlement_variables(3) ~= 1000
       error('Regular and clustered spatial arrangements only available for 100, 200, 500 and 1000 households. Current number of households specified is %f', settlement_variables(3))
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETTLEMENT TYPE DECLARATION - END    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BIOLOGICAL PARAMETERS BEGIN %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify relationship between infected dog density and impact on IP surface
%Have two regions, one in local neighbourhood of max effect, one of redcued
%effect further away

% inf_dog_IP_effect_range Entry definitions
% 1: Radius of full effect region (0.3 matches flight range of sandfly!)
% 2: Size of reduced effect region (whihc is a ring around full effect region)
inf_dog_IP_effect_range = [dog_rad_of_effect,0];


%% Flag variables to determine how sandfly abundance is allocated, 
sandfly_alloc_flag = 1; %1 - sample uniformly from data, 2 - sample from distribution fitted to data


%% alter_sheds_popns_everysimn_flag entry definitions
%% Indicator variables: 0 - fixed across all simulation runs; 1 - resampled each run

% 1: human/dog/chicken distributions
% 2: sandfly distribution
alter_popns_everysimn_flag = [1 1];

%% Specify how baseline host preference is calculated 

host_pref_var = 'simple';

if strcmp(host_pref_var,'simple') == 0 && strcmp(host_pref_var,'complex') == 0
    error('Accepted string inputs for host preference function are simple and complex');
end

%% DOG INFECTION STATUS RELATED PARAMETERS %%
% e = latent, i = highly infectious, a = low infectious, n = never infectious, s = susceptible, d = dead
% 1: dog e -> i, a or n
% 2: dog a -> d,  (lowly infectious class death rate)
% 3: dog i -> d (highly infectious death rate)
% 4: dog e -> d, n -> d (exposed and never infectious death rate)
% 5: dog s -> d (standard death rate)
% 6: dog d -> s. So rate of getting a new dog once a dog has died. 
% 7: Proportion of immigrant dogs seropositive (move directly to one of the three infectious classes) 
d_rates = [1/(dog_latent_period) (12*dog_LowInf_MortRate)/365 (12*dog_HighInf_MortRate)/365 (12*dog_Exp_NeverInf_MortRate)/365 ...
    (12*dog_Sus_MortRate)/365 1/(dog_Replacement_Time) propn_newborn_dog_inf];
 
%% INFECTED DOGS - INFECTIOUSNESS CLASS PROPORTIONS %%
% 1: highly infectious
% 2: low infectious
% 3: never infectious

%Input parameters used here: neverinf_dog_propn, propn_infectious_dogs_highinf 

%Calculate proportion of dogs in highinf and lowinf classes, assign to
%vector
highinf_dog_propn = (1 - neverinf_dog_propn)*propn_infectious_dogs_highinf;
lowinf_dog_propn = 1 - neverinf_dog_propn - highinf_dog_propn;
infectious_dog_class_propns = [highinf_dog_propn, lowinf_dog_propn, neverinf_dog_propn];

%%  RATE OF INFECTED SANDFLY BITE PARAMETER VALUES %%%
bite_rate = Sandfly_BiteRate;
baseline_propn_inf_flies = Baseline_Sandfly_InfPropn;

%%  HOST DEPENDENT PROBABILITY OF INFECTION WHEN BITTEN BY INFECTED SANDFLY %%%

% prob. of host acquiring infection from infected sandfly bite, entry defs
% 1: dogs
% 2: adult
% 3: child
prob_host_acquiring_inf = [Prob_SusDog_AcquireInf,0.1605,0.1605];

%% Parameters/functions to alter proportion of infected sandflies in regions with infected dogs 

%Probability of susceptible fly becoming infected when biting an infectious
%dog 
%NOTE, THIS IS AVERAGED OVER BOTH HIGH AND LOW INFECTIOUS DOGS
% Given by parameter: Propn_Sandflies_InfbyInfDogs;

relative_propn_dog_infectious_class = infectious_dog_class_propns(1:2) ./ sum(infectious_dog_class_propns(1:2));

%Calculate thresholds for highly infectious and low infectious 
%Thresholds correpond to max attainable proportion of sandflies infected, if 
%that type of infectious dog type was the only host present in region of influence
%NOTE, ASSUME 80% OF TRANSMISSION IS DUE TO HIGHLY INFECTIOUS DOGS
highinf_dog_sandfly_inf_prob = 0.8*Propn_Sandflies_InfbyInfDogs /relative_propn_dog_infectious_class(1);
lowinf_dog_sandfly_inf_prob = 0.2*Propn_Sandflies_InfbyInfDogs /relative_propn_dog_infectious_class(2);

%Assign to variable (for function input)
%Entry 1 - threshold for highly infectious dogs
%Entry 2 - threshold for low infectious dogs
sandfly_inf_propn_thresholds = [highinf_dog_sandfly_inf_prob,lowinf_dog_sandfly_inf_prob];

if any(sandfly_inf_propn_thresholds > 1)
    error('infectious_dog_class_propns, Propn_Sandflies_InfbyInfDogs values give sandfly_inf_propn_thresholds above 1, so are not realistic!')
end

%% Sandfly population scaling calculation
abundance_scaling_param = 1/(1 - Propn_flies_NotObvs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BIOLOGICAL PARAMETERS END %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulation type variables
run_num = 5; %Number of runs (for published study, set to 1000 per parameter set)
maxt = 10*365; %Time period to be simulated per run
burn_in_time = 2*365; %Specify time when data is not recorded
timestep = 1;
simulation_start_date = 1; %Allow simulation to begin at a specified time of the year!
RunNum_SaveStartIdx = 100; %Index used when saving MAT File (to avoid partially overwriting files from other runs, should use multiples of 50)

%% Function call
VL_model_function(settlement_variables,inf_dog_IP_effect_range,sandfly_alloc_flag,...
alter_popns_everysimn_flag,host_pref_var,abundance_scaling_param,...
d_rates,infectious_dog_class_propns,bite_rate,baseline_propn_inf_flies,prob_host_acquiring_inf,...
sandfly_inf_propn_thresholds,...
run_num,maxt,burn_in_time,timestep,simulation_start_date,RunNum_SaveStartIdx,job_num);
