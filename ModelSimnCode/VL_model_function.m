%% LEISHMANIASIS CODE - function call version

%Accept as inputs the settlement config, biological parameters
%Input definitions provided in VL_model_SimnCallScript

function VL_model_function(settlement_variables,inf_dog_IP_effect_range,sandfly_alloc_flag,...
alter_popns_everysimn_flag,host_pref_var,abundance_scaling_param,...
d_rates,infectious_dog_class_propns,bite_rate,baseline_propn_inf_flies,prob_host_acquiring_inf,...
sandfly_inf_propn_thresholds,...
num_runs,maxt,burn_in_time,timestep,simulation_start_date,RunNum_SaveStartIdx,job_num)


%% HOUSEHOLD LOCATIONS %%
% Set household locations using set_household_locations function
%Outputs: MAT_filename for saving the results with the current parameter values
%and d the matrix of distances between all the households
[household_num, d] = set_household_locations(settlement_variables);
settlement_type = settlement_variables(2);

%Specify MAT filename for output file
MAT_filename = (['OutputMATfile_#',num2str(RunNum_SaveStartIdx + job_num),'.mat']);

%% HOUSEHOLD HOST NUMBERS - - FIXED ACROSS ALL RUNS %%
% Set host locations using hostDist function
% Outputs: vectors of the number of humans (split into two age categories), dogs, chickens per household
% If don't want distributions informed by data, use placeholder vectors
% with constant number of hosts per household. 

if alter_popns_everysimn_flag(1) == 0
    
    [adults_per_household] = hostDist(household_num, 0, settlement_type);
    %adults_per_household = 0*ones(household_num,1);  %Placeholder vector: everything has 0 humans
    
    [children_per_household] = hostDist(household_num, 1, settlement_type);
    %children_per_household = 0*ones(household_num,1);  %Placeholder vector: everything has 0 humans
    
    [dogs_per_household] = hostDist(household_num, 2, settlement_type);
    % dogs_per_household = 2*ones(household_num,1);  %Placeholder vector: everything has 2 dogs
    
    [chickens_per_household] = hostDist(household_num, 3, settlement_type);
    % chickens_per_household = 10*ones(household_num,1); %Placeholder vector: everything has 10 chickens
end

%% SANDFLY ABUNDANCE NUMBERS 

%Generate new sandfly abundance at start of each run if required
%sandfly_alloc_flag - determines method by which abundance is allocated 
if alter_popns_everysimn_flag(2) == 0
    [sandfly_household_abundance] = assign_sandfly_abundance(household_num,sandfly_alloc_flag,abundance_scaling_param);
end

%% Timing parameters
t = 1:timestep:maxt; % time vector, accounts for length of timestep
burn_in_timestep_num = floor(burn_in_time/timestep); %Get number of timesteps required to leave burn in period
num_timesteps_tracked = numel(t) - burn_in_timestep_num; %Used to size prevalence count vectors. Number of timesteps where info. is recorded


%% MATRICES NEEDED FOR THE SIMULATION

dog_num_inf = zeros(3, num_runs); % recording number of dog infection events per simulation run: Col 1 - low inf, Col 2 - high inf, Col 3 - Never inf 
dog_import_num_inf = zeros(3, num_runs); % recording number of dog importations that are infected per simulation run: Col 1 - low inf, Col 2 - high inf, Col 3 - Never inf 

%Vectors for number of at-risk dogs 
unique_dog_count = zeros(1,num_runs); %Across entire run
AtRiskDogsPerYear = zeros(num_runs,floor(num_timesteps_tracked/365)); %Yearly at-risk dog count

%Initialise storage arrays
dog_sus_prev = zeros(num_runs,num_timesteps_tracked); %Total number of susceptible dogs per day (prevalence count)  
dog_exposed_prev = zeros(num_runs,num_timesteps_tracked); %Total number of exposed dogs per day (prevalence count)  
dog_neverinf_prev = zeros(num_runs,num_timesteps_tracked); %Total number of never infectious dogs per day (prevalence count)  
dog_lowinf_prev = zeros(num_runs,num_timesteps_tracked); %Total number of low infectious dogs per day (prevalence count)  
dog_highinf_prev = zeros(num_runs,num_timesteps_tracked); %Total number of highly infectious dogs per day (prevalence count) 

%Note, the next three arrays track infected cases that are not importations
dog_low_inf_per_day = zeros(num_runs,num_timesteps_tracked); %Number of new low infectious dogs (summed over all households) per day (incidence count) 
dog_high_inf_per_day = zeros(num_runs,num_timesteps_tracked); %Number of new highly infectious dogs (summed over all households) per day (incidence count) 
dog_never_inf_per_day = zeros(num_runs,num_timesteps_tracked); %Number of new never infectious dogs (summed over all households) per day (incidence count) 

dog_imports_low_inf_per_day = zeros(num_runs,num_timesteps_tracked); %Number of new imported low infectious dogs (summed over all households) per day (incidence count) 
dog_imports_high_inf_per_day = zeros(num_runs,num_timesteps_tracked); %Number of new imported highly infectious dogs (summed over all households) per day (incidence count) 
dog_imports_never_inf_per_day = zeros(num_runs,num_timesteps_tracked); %Number of new imported never infectious dogs (summed over all households) per day (incidence count) 

total_species_popn = zeros(3,num_runs); %Number of dogs/adults/children in each simulation (used for calibration purposes)

%% START LOOPS
%fprintf('delta = %3.1f, radius= %3.1fkm, longevity = %d days, %d sheds.\n', delta_s, pherom_dist_var(1), pherom_longev_var(1), num_sheds_s); % output info
%fprintf('%d sheds with pheromone\n',  num_sheds_s); % output info

for loop_num = 1:num_runs% loop through simulations
    tic
    
    %Generate new host distribution at start of each run if required
    if alter_popns_everysimn_flag(1) == 1
        
        [adults_per_household] = hostDist(household_num, 0, settlement_type);
        %adults_per_household = 0*ones(household_num,1);  %Placeholder vector: everything has 0 humans
        
        [children_per_household] = hostDist(household_num, 1, settlement_type);
        %children_per_household = 0*ones(household_num,1);  %Placeholder vector: everything has 0 humans
        
        [dogs_per_household] = hostDist(household_num, 2, settlement_type);
        % dogs_per_household = 2*ones(household_num,1);  %Placeholder vector: everything has 2 dogs
        
        [chickens_per_household] = hostDist(household_num, 3, settlement_type);
        % chickens_per_household = 10*ones(household_num,1); %Placeholder vector: everything has 10 chickens   
        
    end
    
    %Generate new sandfly abundance at start of each run if required
    %sandfly_alloc_flag - determines method by which abundance is allocated 
    if alter_popns_everysimn_flag(2) == 1
        [sandfly_household_abundance] = assign_sandfly_abundance(household_num,sandfly_alloc_flag,abundance_scaling_param);
    end
    
    %Assign populations for current run 
    total_species_popn(1,loop_num) = sum(dogs_per_household);
    total_species_popn(2,loop_num) = sum(adults_per_household);
    total_species_popn(3,loop_num) = sum(children_per_household);

    %% Define storage for this loop
    dog_highly_inf_counter = zeros(household_num,1); %Vector to count number of highly infectious dog infection events in each household
    dog_low_inf_counter = zeros(household_num,1); %Vector to count number of low infection dog infection events in each household
    dog_never_inf_counter = zeros(household_num,1); %Vector to count number of never infectious dog infection events in each household
    
    %Vectors to count number of imported infected dogs of each type 
    dog_import_highly_inf_counter = zeros(household_num,1); 
    dog_import_low_inf_counter = zeros(household_num,1); 
    dog_import_never_inf_counter = zeros(household_num,1); 
    
    dog_stat = zeros(household_num, 6, numel(t)); % dog infection status matrix per house per timestep: Col. 1 - Suscep, Col.2 - NeverInf, Col.3 - Exposed, Col.4 - Low Infectious, Col. 5 - Highly Infectious, Col. 6 - dead.

    %History arrays to use in visual simulation!
    ip_history = zeros(household_num, numel(t)); %IP at each house per timestep
    
    
    %% Initialise number of suceptibles at initial time
    dog_stat(:,1,1) = dogs_per_household; %All dogs start susceptible 

    %% SIMULATE EPIDEMIC
    simulation_day_of_year = simulation_start_date; 
    record_count = 0;
    year_idx = 1; %Vector indexing term used with AtRiskDogsPerYear, track year post burn-in
    for i=2:numel(t) %Cycle through all days from 2 to the end (day one is where everything is susceptible)
        %% Event update storage vectors: keep track of infection status changes that happen this day per household
                     
        %%% METHOD DESCRIPTION: ON EACH TIMESTEP %%%
        % (i) Set up IP surface   
        % (ii) Run epidemic

        %% (i) Set up IP surface      
        
        % First compute baseline host preference values at each household
        % based on dog, chicken, human populations.
        % This has to be done each day because some dogs can die. 
    
        alive_dogs_per_household = sum(dog_stat(:, 1:5, i-1),2); %Sum each row, ignore dead dogs (column 6)   
        [hostpref_household,overall_mass_per_household] = host_pref_allocation(host_pref_var,...
            chickens_per_household,alive_dogs_per_household,...
            adults_per_household,children_per_household,household_num);

        %On the first day, define yesterdays abundance as the max sandfly abundance 
        %possible yesterday for each house, otherwise we have it as an
        %output of the IP surface function
        if(i == 2)
            load('Host_dist_data/seasonality365_start01012011.mat');
            
            if(simulation_day_of_year == 1)
                yesterday_day_of_year = 365;
            else
                yesterday_day_of_year = simulation_day_of_year - 1;
            end
            
            seasonality_scaling_yest = seasonalityVec(yesterday_day_of_year)/max(seasonalityVec)+0.5;
            yest_abundance = sandfly_household_abundance*seasonality_scaling_yest; 
        end
        
        [household_level_dog_ip_vector,yest_abundance] = IP_setup(household_num, bite_rate, baseline_propn_inf_flies,...
            prob_host_acquiring_inf, yest_abundance,...
            simulation_day_of_year, sandfly_household_abundance, hostpref_household, d,...
            dog_stat(:, :, i-1), inf_dog_IP_effect_range,...
            sandfly_inf_propn_thresholds,overall_mass_per_household);
        
        %Multiply the 5 elements of the ip vector together for each host
        household_level_dog_ip = household_level_dog_ip_vector(:, 1).*...
            household_level_dog_ip_vector(:, 2).*household_level_dog_ip_vector(:, 3).*...
            household_level_dog_ip_vector(:, 4).*household_level_dog_ip_vector(:, 5);
        
        %Split infectious pressure by number of dogs at each household
        household_level_dog_ip(alive_dogs_per_household > 0) = household_level_dog_ip(alive_dogs_per_household > 0)  ./ alive_dogs_per_household(alive_dogs_per_household > 0);
        household_level_dog_ip(alive_dogs_per_household == 0) = 0; %No dogs present, so nothing to infect
        dog_ip = household_level_dog_ip;        

        %% (ii) Run epidemic and update matrices
     
        %Run Epidemic
        [dog_event_update,...
            d_exp_to_low_inf_num, d_exp_to_high_inf_num,d_exp_to_never_inf_num,...
            d_dead_to_low_inf_num, d_dead_to_high_inf_num,d_dead_to_never_inf_num,...
            d_dead_to_reborn_num]= run_epidemic(dog_stat(:,:,i-1),dog_ip,...
            d_rates, infectious_dog_class_propns, household_num, timestep);
      
        %Update infection status matrices with this timesteps changes
        dog_stat(:,:,i) = dog_stat(:,:,i-1) + dog_event_update;
        
        %If burn_in_time exceeded, update count vectors
        if t(i) > burn_in_time
            record_count = record_count + 1;
            
            %Update unique dog counter
            %If first iteration after burn_in_time, initialise with number of dogs susceptible at
            %end of previous timestep
            if record_count == 1
                unique_dog_count(1,loop_num) = sum(dog_stat(:,1,i-1));
                unique_dog_count(1,loop_num) = unique_dog_count(1,loop_num) + sum(d_dead_to_reborn_num);
            elseif record_count ~= num_timesteps_tracked %Newborns from final timestep have no chance of being infected, so not included in calculation
                unique_dog_count(1,loop_num) = unique_dog_count(1,loop_num) + sum(d_dead_to_reborn_num);
            end
            
            %Update at-risk population per year counts
            if mod(record_count,365) == 1
                AtRiskDogsPerYear(loop_num,year_idx) = sum(dog_stat(:,1,i-1)); %Count number of susceptibles at start of year
                AtRiskDogsPerYear(loop_num,year_idx) = AtRiskDogsPerYear(loop_num,year_idx) + sum(d_dead_to_reborn_num); %Newborns can not be infected until next timestep!
            elseif mod(record_count,365) == 0
                %Newborns at end of year can be first infected at start of
                %next year. This info is contained in susceptible count at year start
                year_idx = year_idx + 1;
            else
                AtRiskDogsPerYear(loop_num,year_idx) = AtRiskDogsPerYear(loop_num,year_idx) + sum(d_dead_to_reborn_num);
            end
            
            %Update counters of number of infection events over whole
            %simulation
            dog_low_inf_counter = dog_low_inf_counter + d_exp_to_low_inf_num;
            dog_highly_inf_counter = dog_highly_inf_counter + d_exp_to_high_inf_num;           
            dog_never_inf_counter = dog_never_inf_counter + d_exp_to_never_inf_num;           

            %Update infection importation counters
            dog_import_low_inf_counter = dog_import_low_inf_counter + d_dead_to_low_inf_num;
            dog_import_highly_inf_counter = dog_import_highly_inf_counter + d_dead_to_high_inf_num;           
            dog_import_never_inf_counter = dog_import_never_inf_counter + d_dead_to_never_inf_num;  
            
            %Update prevalence counters
            dog_sus_prev(loop_num,record_count) = sum(dog_stat(:,1,i));
            dog_neverinf_prev(loop_num,record_count) = sum(dog_stat(:,2,i));
            dog_exposed_prev(loop_num,record_count) = sum(dog_stat(:,3,i));
            dog_lowinf_prev(loop_num,record_count) = sum(dog_stat(:,4,i));
            dog_highinf_prev(loop_num,record_count) = sum(dog_stat(:,5,i));

            %Update number of infection events per day per simulation loop (for
            %all houses)
            dog_low_inf_per_day(loop_num,record_count) =  sum(d_exp_to_low_inf_num);
            dog_high_inf_per_day(loop_num,record_count) =  sum(d_exp_to_high_inf_num);
            dog_never_inf_per_day(loop_num,record_count) =  sum(d_exp_to_never_inf_num);       
            
            dog_imports_low_inf_per_day(loop_num,record_count) =  sum(d_dead_to_low_inf_num);
            dog_imports_high_inf_per_day(loop_num,record_count) =   sum(d_dead_to_high_inf_num);
            dog_imports_never_inf_per_day(loop_num,record_count) =   sum(d_dead_to_never_inf_num);                  

        end
        
        %ERROR CHECK
        if nnz(dog_stat(:,:,i) < 0) > 0
            error('Number of dogs has gone negative!');
        end
        
        %ERROR CHECK
        if nnz(sum(dog_stat(:,:,i),2) ~= dogs_per_household) > 0
            error('Number of hosts a households has changed from original value!');
        end
                
        %Update arrays to use in visual simulation
        ip_history(:,i) = dog_ip;
        
        %Update day of year tracker
        simulation_day_of_year = mod(simulation_day_of_year + timestep,365);
        if simulation_day_of_year == 0
            simulation_day_of_year = 365; %Deal with case of last day of year resetting to zero!
        end
        
    end
    
    %Record total number of low, high and never infectious dogs 
    %NOTE - NOT IMPORTED INFECTIONS    
    dog_num_inf(1, loop_num) = sum(dog_low_inf_counter); 
    dog_num_inf(2, loop_num) = sum(dog_highly_inf_counter);
    dog_num_inf(3, loop_num) = sum(dog_never_inf_counter);
    
    %Record total number of low, high and never infectious dogs
    %IMPORTED INFECTED CASES
    dog_import_num_inf(1, loop_num) = sum(dog_import_low_inf_counter);
    dog_import_num_inf(2, loop_num) = sum(dog_import_highly_inf_counter);
    dog_import_num_inf(3, loop_num) = sum(dog_import_never_inf_counter);
    
    
    %% Display simulation number every 100 simulations
    if mod(loop_num, 10) == 0
        fprintf('%d loops done.\n', loop_num);
    end
    toc
end

%Save num_inf array to file
save(MAT_filename,'dog_num_inf','dog_import_num_inf',...
    'dog_sus_prev','dog_exposed_prev','dog_lowinf_prev','dog_highinf_prev','dog_neverinf_prev',...
    'dog_low_inf_per_day','dog_high_inf_per_day','dog_never_inf_per_day',...
    'dog_imports_low_inf_per_day','dog_imports_high_inf_per_day','dog_imports_never_inf_per_day',...
    'total_species_popn','unique_dog_count','AtRiskDogsPerYear');
end
