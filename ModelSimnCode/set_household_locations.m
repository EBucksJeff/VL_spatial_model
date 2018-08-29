% This sets the locations of the households for the simulation

function [household_num, d] = set_household_locations(settlement_variables)
%Outputs: 
%   household_num 
%   d - array giving distances between households 
%Inputs 
%   settlement variables -  as assigned below
settlement_config = settlement_variables(1); 
household_num = settlement_variables(3); 
grid_size = settlement_variables(4); 

if settlement_config == 1
    %Set up household locations from data
    
    %Import lat/long data from file
    filename_rural_Brazil = 'Spatial_configuration_files/Brazil_location_data/Brazil_rural_village_Calderao_235_households.mat';
    load(filename_rural_Brazil)
    x = village_longitude_CL;
    y = village_latitude_CL;
       
    household_num = length(x); %Number of households
    
elseif settlement_config == 2
    % Set up household locations not from data: random placement on a grid
    x = grid_size*rand(household_num, 1); y = grid_size*rand(household_num, 1);
     
elseif settlement_config == 3
    %Set up household locations not from data: regular placement on a grid
    
    %Load household locations from MAT file
    % Options are currently limited to 100, 200, 500, 10000 houses
    % If want more need to run code: in dropbox
    filename_regular = (['Spatial_configuration_files/Regular_locations_N_',num2str(household_num),'_Xlength_10_Ylength_10.mat']);
    load(filename_regular)
    
    %Houses generated on 10x10 grid, scale to match grid size being used here
    scale_factor_adjustment = 10/grid_size;
    x = x_list/scale_factor_adjustment; y = y_list/scale_factor_adjustment;
        
elseif settlement_config == 4
    %Set up household locations not from data: large clusters
    
    %Load household locations from MAT file.
    % Options are currently limited to 100, 200, 500, 1000 houses
    % If want more need to run code: in dropbox
    filename_large_cluster = (['Spatial_configuration_files/Spatialarrangement_N_',num2str(household_num),'_B_0.4_ratio_50.mat']);
    load(filename_large_cluster)
    
    %Cluster configuration generated on 50x50 grid, scale to match grid size being used here
    scale_factor_adjustment = 50/grid_size;
    x = x_location/scale_factor_adjustment; y = y_location/scale_factor_adjustment;
        
elseif settlement_config == 5
    % Set up household locations not from data: small clusters
    
    %Load household locations from MAT file.
    % Options are currently limited to 100, 200, 500, 1000 houses
    % If want more need to run code: in dropbox
    filename_small_cluster = (['Spatial_configuration_files/Spatialarrangement_N_',num2str(household_num),'_B_1_ratio_50.mat']);
    load(filename_small_cluster)
    
    %Cluster configuration generated on 50x50 grid, scale to match grid size being used here
    scale_factor_adjustment = 50/grid_size; 
    x = x_location/scale_factor_adjustment; y = y_location/scale_factor_adjustment;

    
else
    error('No spatial configuration has been selected!')
end



%% Calculate distances between households
% d is matrix (size household_num x household_num) giving the distance
% between any two households

% If using real data then compute the distance using the latitude and longitude
if settlement_config == 1
    x_array_1 = x*ones(1, household_num); %Row i gives longitude of household i
    x_array_2 = ones(household_num,1)*x'; %Column j gives longitude of household j
    y_array_1 = y*ones(1, household_num); %Row i gives latitude of household i
    y_array_2 = ones(household_num,1)*y'; %Column j gives latitude of household j
    d = great_circle_distance_calc(y_array_1,y_array_2,x_array_1,x_array_2); 
    
    %Error check
    if sum(isnan(d(:)) == 1)
       error('NaN values in distance matrix, location data is missing!'); 
    end
    
else %if using synthetic household locations
    d = sqrt((x*ones(1, household_num) - ones(household_num, 1)*x').^2 + (y*ones(1, household_num) - ones(household_num, 1)*y').^2);
end


end

%Calculate the great circle distance (premises to premises)
function dist = great_circle_distance_calc(lat1,lat2,lon1,lon2)

radius=6371;            %Earths Radius

lat1=lat1*pi/180;
lat2=lat2*pi/180;
lon1=lon1*pi/180;
lon2=lon2*pi/180;
deltaLat=lat2-lat1;
deltaLon=lon2-lon1;
a=sin((deltaLat)/2).^2 + cos(lat1).*cos(lat2).* sin(deltaLon/2).^2;
c=2*atan2(sqrt(a),sqrt(1-a));
%c=2*asin(min(ones(length(a),1),sqrt(a)));

dist = radius.*c;

end
