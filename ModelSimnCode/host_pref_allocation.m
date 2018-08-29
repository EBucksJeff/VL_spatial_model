%Computes preference towards each host type (household specific)

function [hostpref_household,overall_mass_per_household] =...
    host_pref_allocation(host_pref_baseline_choice,chickens_per_household...
    ,dogs_per_household,adults_per_household,children_per_household,household_num)
%Inputs:
%   host_pref_baseline_choice: Choose how to calculate host preferece. Either using 'simple' or 'complex' biomass dependence function.
%   Other inputs are host numbers at each household, the number of households.

%Outputs: 
%   hostpref_household: array of household level host preference. 
%   (Note, Row i - corresponds to household i ,Column - host type  [%Column 1 - Dogs, Column 2 - Adults, Column 3 - Children])
%   overall_mass_per_household: vector of biomass amount at each household 

if strcmp(host_pref_baseline_choice,'simple') == 1 %"Simple" biomass dependence function 
    %Biomass values (units are equivalent num. of chickens)
    child_mass = 5;
    adult_mass = 10;
    dog_mass = 2;

    %Get total biomass of each host type at every household
    child_mass_per_household = children_per_household*child_mass;
    adult_mass_per_household = adults_per_household*adult_mass;
    dog_mass_per_household = dogs_per_household*dog_mass;
    
    %Get total biomass at every household & index for households with
    %no hosts present (possible if have no numans in model)
    overall_mass_per_household = child_mass_per_household + adult_mass_per_household + dog_mass_per_household + chickens_per_household;
    mass_present = find(overall_mass_per_household > 0);
    
    %Calculate host preference for each household
    hostpref_household = zeros(household_num,3); %Column 1 - Dogs, Column 2 - Adults, Column 3 - Children
    hostpref_household(mass_present,1) = dog_mass_per_household(mass_present)./overall_mass_per_household(mass_present);
    hostpref_household(mass_present,2) = adult_mass_per_household(mass_present)./overall_mass_per_household(mass_present);
    hostpref_household(mass_present,3) = child_mass_per_household(mass_present)./overall_mass_per_household(mass_present);
end

if strcmp(host_pref_baseline_choice,'complex') == 1 %"Complex" biomass dependence function 
    
    chicken_to_child_conversion = @(x) 0.2*x;
    
    %Biomass values (units are equivalent num. of chickens)
    dog_mass = 2;
    
    %Get prefrence towards each host type at every household (in units of
    %children)
    child_mass_per_household = children_per_household;
    adult_mass_per_household = adults_per_household*2; %Treat adults as twice the mass of a child
    chicken_mass_per_household = chicken_to_child_conversion(chickens_per_household);
    dog_mass_per_household =(chicken_to_child_conversion(dogs_per_household))*dog_mass;
    
    %Get total biomass at every household & index for households with
    %no hosts present (possible if have no numans in model)
    overall_mass_per_household = child_mass_per_household + adult_mass_per_household + dog_mass_per_household + chicken_mass_per_household;
    mass_present = find(overall_mass_per_household > 0);
    
    %Calculate host preference for each household
    hostpref_household = zeros(household_num,3); %Column 1 - Dogs, Column 2 - Adults, Column 3 - Children
    hostpref_household(mass_present,1) = dog_mass_per_household(mass_present)./overall_mass_per_household(mass_present);
    hostpref_household(mass_present,2) = adult_mass_per_household(mass_present)./overall_mass_per_household(mass_present);
    hostpref_household(mass_present,3) = child_mass_per_household(mass_present)./overall_mass_per_household(mass_present);
end

end
