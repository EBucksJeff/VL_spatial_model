%Script to return vector of values for histrogram for parameter i, value j
%all others at baseline
function [av_prev_per_sim] = return_histogram_values_prev(i, j)

    %%% Import data files for parameter combination of interest
    if(j == 1)
        simn_idx = 1; %Baseline parameter values
    else
        simn_idx = (i-1)*3+j;
    end
    
    filename = (['../SimnOutputFiles/VLmodel_BiolParamSensitivitySimn_Run#' num2str(simn_idx) '.mat']);
    load(filename); 


    %The following is adjusted from the plot code script.
    %Gives the prevalence averaged over the last year for each simulation.

    % indices for the last year
    prev_indices = max(1, size(dog_highinf_prev, 2) - 364):size(dog_highinf_prev, 2); 

    % prevalence per day for the last year
    av_prev_per_sim = mean((dog_highinf_prev(:, prev_indices) + dog_lowinf_prev(:, prev_indices)...
        + dog_neverinf_prev(:, prev_indices) + dog_exposed_prev(:, prev_indices))...
    ./(dog_highinf_prev(:, prev_indices) + dog_lowinf_prev(:, prev_indices) ...
    + dog_neverinf_prev(:, prev_indices) + dog_sus_prev(:, prev_indices)...
    + dog_exposed_prev(:, prev_indices)), 2);
    
end
