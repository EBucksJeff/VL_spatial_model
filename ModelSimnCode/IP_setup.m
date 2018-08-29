function [household_level_dog_ip_vector,yest_abundance] = IP_setup(numHouses, bite_rate, baseline_propn_inf_flies,...
    prob_host_acquiring_inf, yest_abundance,...
    simulation_day_of_year, sandfly_household_base_abundance, hostpref_household, d,... 
    dog_stat, inf_dog_IP_effect_range,...
    sandfly_inf_propn_thresholds,overall_mass_per_household)
%Input explainer:
%
% prob_host_acquiring_inf:(from infected sandfly bite)
%                           Entry 1 -> dogs
%                           Entry 2 -> adult
%                           Entry 3 -> child
%
% yest_abundance: Sandfly abundances per household at previous timestep 
%
% sandfly_household_base_abundance: Baseline sandfly abundance per household, later scaled by seasonality
%
% hostpref_household: array of household level host preference. 
%                       Row i - corresponds to household i ,
%                      Column - host type  [Column 1 - Dogs, Column 2 - Adults, Column 3 - Children]
%
% d: array giving distances between households
%
% dog_stat: matrix is from the previous
%          day i-1. So in here the status matrices are just 2D not 3D
%
% inf_dog_IP_effect_range: Relationship between infected dog density and impact on IP surface
%                            % Entry definitions
%                            % Entry 1 -> Radius of full effect region 
%                            % Entry 2 -> Size of reduced effect region (whihc is a ring around full effect region)
%
% sandfly_inf_propn_thresholds: Max attainable proportion of sandflies infected, if 
%                               that type of infectious dog type was the only host present 
%                               in region of influence.
%                               Entry 1 - threshold for highly infectious dogs
%                               Entry 2 - threshold for low infectious dogs
% overall_mass_per_household: vector of biomass amount at each household 


%Outputs:
% household_level_dog_ip_vector: vectors of infectious pressure (vectors of length 5)
% yest_abundance: Output the abundance of sandflies for tomorrow

%%% METHOD DESCRIPTION: SETTING UP INFECTIOUS PRESSURE %%%
% (i) Compute maximum sandfly abundance per house from seasonal data and
% grow the abundance at any houses which are not at maximum using
% yesterdays abundance
% (ii) Update proportion of infected sandflies based on the presence of
% infected dogs
% (iii) Any abundances that are too low put back up to maximum


%% (i) Compute maximum sandfly abundance per house from seasonal data and
% grow the abundance at any houses which are not at maximum using
% yesterdays abundance%Get possible maximum abundance for yesterday from seasonality and data
load('Host_dist_data/seasonality365_start01012011.mat');

if(simulation_day_of_year == 1)
    yesterday_day_of_year = 365;
else
    yesterday_day_of_year = simulation_day_of_year - 1;
end
seasonality_scaling_yest = seasonalityVec(yesterday_day_of_year)/max(seasonalityVec)+0.5;
max_abundance_yest = sandfly_household_base_abundance*seasonality_scaling_yest;

%Get possible maximum abundance today from seasonality and data
seasonality_scaling = seasonalityVec(simulation_day_of_year)/max(seasonalityVec) + 0.5;
max_abundance = sandfly_household_base_abundance*seasonality_scaling;

today_abundance = max_abundance;
for i=1:numHouses
    if(yest_abundance(i) < max_abundance_yest(i))
        L = max_abundance(i); %Maximum to grow to (today's maximum)
        y = yest_abundance(i); %Current level
        logistic_day = -10*log((L-y)/y); %Function this! %Find the day using the inverse 
        
        %Add one to the day and compute the new value using logistic growth
        logistic_day = logistic_day+1;
        potential_abundance = L/(1+exp(-0.1*logistic_day));
        if(potential_abundance < max_abundance(i))
            today_abundance(i) = potential_abundance;
        end
    end
end


%% (ii) Update proportion of infected sandflies based on the presence of infected dogs
% Assume additive effect from multiple dogs 

%Find households with infected dogs
dog_highly_inf_idx = find(dog_stat(:,5) > 0); 
dog_low_inf_idx = find(dog_stat(:,4) > 0);

%For each household, find all other households in sandfly travel distance range of effect 
house_in_range_indicator = d <= inf_dog_IP_effect_range(1);

%Calculate weighting factor for those households within sandfly range
%ASSUMING LINEAR DECREASE WITH DISTANCE
%At household of interest, d = 0,  so distance_weighting will be 1. 
%For d>effect_range, house_inf_range_indicator will be zero, giving distance_weighting of 0
distance_weighting = ((inf_dog_IP_effect_range(1) - d) ./inf_dog_IP_effect_range(1)) .*house_in_range_indicator;  

if any(distance_weighting > 1)  %Check for implausible values
   error('Sandfly infected proportion distance weighting exceeds 1!') 
end

if any(distance_weighting < 0) %Check for implausible values
   error('Sandfly infected proportion distance weighting has gone below 0!') 
end

dog_HighInfStatus_replicate = ones(numHouses,1)*dog_stat(:,5)';
dog_LowInfStatus_replicate = ones(numHouses,1)*dog_stat(:,4)';
biomass_replicate = ones(numHouses,1)*overall_mass_per_household';

%For each household, number of highly infectious/low infectious dogs within range that impact that household IP.
num_HighInf_dog_within_range = sum(dog_HighInfStatus_replicate(:,dog_highly_inf_idx).*distance_weighting(:,dog_highly_inf_idx),2); 
num_LowInf_dog_within_range = sum(dog_LowInfStatus_replicate(:,dog_low_inf_idx).*distance_weighting(:,dog_low_inf_idx),2); 

%Get total biomass of houses within range, check no household report zero biomass
biomass_within_range = sum(biomass_replicate.*distance_weighting,2);

if any(biomass_within_range == 0)
   error('Biomass of zero found, but every household should have at least one person!') 
end

%Compute relative proportion of biomass that is high inf/low inf dogs per
%household
highinf_dog_biomass_propn = (2*num_HighInf_dog_within_range)./biomass_within_range;
lowinf_dog_biomass_propn = (2*num_LowInf_dog_within_range)./biomass_within_range;

%Input relative biomass proportions and raw numbers of high inf/low inf
%dogs into a function to give updated proportion of infected sandflies

%Get updated proportion of infected sandlfies at each household (entry
%i corresponds to household i)
updated_sandfly_inf_propn = sandfly_inf_propn_func_calc(highinf_dog_biomass_propn,lowinf_dog_biomass_propn,...
    baseline_propn_inf_flies,sandfly_inf_propn_thresholds); 

%% (iii) Fix any abundances that are too low. 

for i=1:numHouses
    if(today_abundance(i) < 0.1*max_abundance(i))
        today_abundance(i) = 0.1*max_abundance(i);
    end
end

%% Infectious pressure vectors for outputting
% Row i corresponds to household i
% Column 1: abundance (at household)
% Column 2: bite rate (fixed across households)
% Column 3: proportion of infected sandflies (at household)
% Column 4: probability of transmission on bite (fixed across households)
% Column 5: host preference scaling (contained in [0, 1])

household_level_dog_ip_vector = [today_abundance repmat(bite_rate, numHouses, 1) updated_sandfly_inf_propn ...
    repmat(prob_host_acquiring_inf(1), numHouses,1) hostpref_household(:, 1)];


yest_abundance = today_abundance; 

end

function modified_sandfly_propn = sandfly_inf_propn_func_calc(highinf_dog_biomass_propn,lowinf_dog_biomass_propn,...
    baseline_fly_inf_propn,sandfly_inf_propn_thresholds)

%PLACEHOLDER FUNCTIONS FOR THE TIME BEING!
%If set to 1, only the proportion of total biomass that is in an infectious dog class matters, absolute number of infectious dogs does not!
highinf_dog_mod_factor =  1;
lowinf_dog_mod_factor =  1;


FlyInfPropn_increase_HighInfDogs = highinf_dog_biomass_propn.*highinf_dog_mod_factor.*(sandfly_inf_propn_thresholds(1) - baseline_fly_inf_propn);
FlyInfPropn_increase_LowInfDogs = lowinf_dog_biomass_propn.*lowinf_dog_mod_factor.*(sandfly_inf_propn_thresholds(2) - baseline_fly_inf_propn);

%Proportion of sandflies infected given by baseline added to increase given
%by highly infectious dogs, increase given by lowly infectious dogs
modified_sandfly_propn = baseline_fly_inf_propn + FlyInfPropn_increase_HighInfDogs + FlyInfPropn_increase_LowInfDogs;

%Error checks
if any(modified_sandfly_propn > sandfly_inf_propn_thresholds(1))
    error('Proportion of flies infected exceeds threshold value assuming only highly infectious dogs are present')
end

end



