function  [sandfly_household_abundance] = assign_sandfly_abundance(household_num,sandfly_alloc_flag,abundance_scaling_param)
%  Allocate sandfly population to each household 

%Inputs:
%   household_num: Number of households. 
%   sandfly_alloc_flag: Flag variable (1 - sample from data;
%                                       2 - fit gamma distribution, draw samples from fitted distribution). 
%   abundance_scaling_param: To account for proportion of female sandflies unobserved. 

%Outputs: 
%   sandfly_household_abundance: Baseline sandfly populations ar each household. 
%   ASSUMPTION: THESE ARE THE BASELINE SANDFLY POPULATIONS AT EACH
%   HOUSEHOLD. SCALED LATER BY SEASONALITY

data_houseF = dlmread('Sandfly_data/SandflyData.txt');
                             
houseF = data_houseF*abundance_scaling_param;

%Use raw data values sampled uniformly
if sandfly_alloc_flag == 1
    r = randi(numel(houseF),household_num,1);
    
    sandfly_household_abundance = houseF(r); %Assign value drawn from data to each household
    
elseif sandfy_alloc_flag == 2
    
    %Fit gamma distribution to data, draw samples from that
    data = houseF;
    pd = fitdist(data, 'Gamma');
    sandfly_household_abundance = round(random(pd, household_num, 1));
end