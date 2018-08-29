function [dog_event_update,...
    d_exp_to_low_inf_num, d_exp_to_highly_inf_num,d_exp_to_never_inf_num,...
    d_dead_to_low_inf_num, d_dead_to_highly_inf_num,d_dead_to_never_inf_num,...
    d_dead_to_reborn_num] =...
    run_epidemic(dog_status, dog_ip,...
    d_rates, infectious_dog_class_propns, household_num, timestep)

%Inputs:
% dog_status: matrix from the previous day, so i-1. So in here the status matrices are just 2D not 3D
% dog_ip: Infectious pressure values against dogs for each household.
% d_rates: Dog-related epidmiological the transition rates. 
%      e = latent, i = highly infectious, a = low infectious, n = never infectious, 
%      s = susceptible, d = dead
%       1: dog e -> i, a or n
%       2: dog a -> d,  (lowly infectious class death rate)
%       3: dog i -> d (highly infectious death rate)
%       4: dog e -> d, n -> d (exposed and never infectious death rate)
%       5: dog s -> d (standard death rate)
%       6: dog d -> s (a, i, or n) rate of getting a new dog once a dog has died. 
%       7: Proportion of immigrant dogs seropositive (move directly to one of the three infectious classes) 
%
% infectious_dog_class_propns: Entry 1 -> highly infectious
%                              Entry 2 -> low infectious
%                              Entry 3 -> never infectious
% household_num, timestep: 

%Outputs:
%  dog_event_update: 2D array. Row per household. Store overall change in 
%                              number of dogs per compartment
%  d_exp_to_low_inf_num etc: Per houshold, number of labelled transition
%                             event to occur on curret timestep

dog_event_update = zeros(household_num, 6); %Col. 1 - Suscep, Col.2 - NeverInf, Col.3 - Exposed, Col.4 - low infectious, Col. 5 - highly infectious, Col. 6 - dead

%% dog s -> d  (death of susceptible)
[d_sus_death_num] = transition(d_rates(5), dog_status(:, 1), household_num, timestep);
dog_event_update(:, 6) = dog_event_update(:, 6) + d_sus_death_num; %moving into dead
dog_event_update(:, 1) = dog_event_update(:, 1) - d_sus_death_num; %moving out of s

%% dog s -> e 
sus_dog_households = find(dog_status(:,1) > 0);

for j = 1:numel(sus_dog_households)
    sus_dog_num = dog_status(sus_dog_households(j),1) - d_sus_death_num(sus_dog_households(j)); %get number of s dogs because some have died
    sus_to_exp_success_num = sum(1 - exp(-dog_ip(sus_dog_households(j),1)*timestep) > rand(sus_dog_num, 1)); % randomly determines what has high enough prob
    
    dog_event_update(sus_dog_households(j),3) = dog_event_update(sus_dog_households(j),3) + sus_to_exp_success_num; % moving into exp. class
    dog_event_update(sus_dog_households(j),1) = dog_event_update(sus_dog_households(j),1) - sus_to_exp_success_num; % moving out of sus. class
end

%% dog e -> d (death and replacement of latent dog)
[d_exp_death_num] = transition(d_rates(4), dog_status(:, 3), household_num, timestep);

dog_event_update(:, 6) = dog_event_update(:, 6) + d_exp_death_num ; % moving into dead
dog_event_update(:, 3) = dog_event_update(:, 3) - d_exp_death_num; % move out of latent

%% dogs e -> i, e -> a, e -> n 
updated_exp_dog_status = dog_status(:, 3) - d_exp_death_num; %To take into account that some have died. 
[d_exp_out_num] = transition(d_rates(1), updated_exp_dog_status, household_num, timestep);

d_exp_to_highly_inf_num = zeros(household_num, 1);
d_exp_to_low_inf_num = zeros(household_num, 1);
moving_households = find(d_exp_out_num > 0);
for j = 1:numel(moving_households)
    r = rand(d_exp_out_num(moving_households(j)), 1);
    d_exp_to_highly_inf_num(moving_households(j)) = sum(infectious_dog_class_propns(1) >= r);
    d_exp_to_low_inf_num(moving_households(j)) = sum(infectious_dog_class_propns(1) < r & (infectious_dog_class_propns(1) + infectious_dog_class_propns(2) >= r));
end
d_exp_to_never_inf_num = d_exp_out_num - d_exp_to_low_inf_num - d_exp_to_highly_inf_num;

dog_event_update(:, 4) = dog_event_update(:, 4) + d_exp_to_low_inf_num; % moving into low infectious
dog_event_update(:, 5) = dog_event_update(:, 5) + d_exp_to_highly_inf_num; % moving into high infectious
dog_event_update(:, 2) = dog_event_update(:, 2) + d_exp_to_never_inf_num; % moving into never infectious
dog_event_update(:, 3) = dog_event_update(:, 3) - d_exp_out_num; % moving out of e

%% dogs a -> d
[d_low_inf_to_dead_num] = transition(d_rates(2), dog_status(:, 4), household_num, timestep);

dog_event_update(:, 6) = dog_event_update(:, 6) + d_low_inf_to_dead_num; % moving into dead
dog_event_update(:, 4) = dog_event_update(:, 4) - d_low_inf_to_dead_num; % moving out of low infectious

%% dogs i -> d
[d_highly_inf_to_dead_num] = transition(d_rates(3), dog_status(:, 5), household_num, timestep);  

dog_event_update(:,6) = dog_event_update(:,6) + d_highly_inf_to_dead_num; % moving into dead
dog_event_update(:,5) = dog_event_update(:,5) - d_highly_inf_to_dead_num; % moving out of high infectious

%% dogs n -> d
[d_neverinf_to_dead_num] = transition(d_rates(4), dog_status(:, 2), household_num, timestep);

dog_event_update(:, 6) = dog_event_update(:, 6) + d_neverinf_to_dead_num; % moving into dead
dog_event_update(:, 2) = dog_event_update(:, 2) - d_neverinf_to_dead_num; % move out of neverinf


%% dogs d -> s, d -> a, d -> i, d -> n
%Find number of new dogs per household introduced this timestep
[d_dead_to_reborn_num] = transition(d_rates(6), dog_status(:, 6), household_num, timestep);

%Initialise vectors - count transitions into each class per household 
d_dead_to_sus_num = d_dead_to_reborn_num;
d_dead_to_highly_inf_num = zeros(household_num, 1);
d_dead_to_low_inf_num = zeros(household_num, 1);

%Determine households with new dogs 
moving_households = find(d_dead_to_reborn_num > 0);
for j = 1:numel(moving_households)
    r = rand(d_dead_to_reborn_num(moving_households(j)), 1);
    
    %%Find number of new dogs (imports) that are already infected
    d_dead_to_inf_num = sum(d_rates(7) >= r);

    if d_dead_to_inf_num > 0 
        %Reduce number of new susceptible dogs, subtract number already infected
        d_dead_to_sus_num(moving_households(j)) = d_dead_to_sus_num(moving_households(j)) - d_dead_to_inf_num;
        
        %Generate new random numbers, used to assign infected dogs to a class
        new_r = rand(d_dead_to_inf_num, 1); 
        d_dead_to_highly_inf_num(moving_households(j)) = sum(infectious_dog_class_propns(1) >= new_r);
        d_dead_to_low_inf_num(moving_households(j)) = sum(infectious_dog_class_propns(1) < new_r & (infectious_dog_class_propns(1) + infectious_dog_class_propns(2) >= new_r));      
    end
end
d_dead_to_never_inf_num = d_dead_to_reborn_num - d_dead_to_sus_num - d_dead_to_low_inf_num - d_dead_to_highly_inf_num;

dog_event_update(:, 6) = dog_event_update(:, 6) - d_dead_to_reborn_num; %moving out of dead
dog_event_update(:, 4) = dog_event_update(:, 4) + d_dead_to_low_inf_num; % moving into low infectious
dog_event_update(:, 5) = dog_event_update(:, 5) + d_dead_to_highly_inf_num; % moving into high infectious
dog_event_update(:, 2) = dog_event_update(:, 2) + d_dead_to_never_inf_num; % moving into never infectious
dog_event_update(:, 1) = dog_event_update(:, 1) + d_dead_to_sus_num; %moving into s.

end

function [to_move] = transition(rate, status_vector, household_num, timestep)
% This takes the number of hosts in each house of the type we're interested in(status_vector)
% the rate of movement between the classes we're interested in
% and the number of households.

% It computes the number of hosts to move from each household

%Find which households we are interested in
households = find(status_vector > 0);
to_move = zeros(household_num, 1);

for j = 1:numel(households)
    number_hosts = status_vector(households(j)); %get number of hosts at this household
    transfer_number = sum(1 - exp(-rate*timestep) > rand(number_hosts, 1)); %use rate to determine how many transfer classes
    
    to_move(households(j)) = transfer_number; % numbers moving per household
end
end
