function [rtn] = hostDist(n, hostType, location)
% Assign population count for each host type to households

%Inputs:
%
%   n: The number of households.
%   hostType: 0 for adults, 1 for children, 2 for dogs, 3 for chickens. 
%   location: 0 for rural, 1 for semi urban, 2 for urban.
%
%   Data used for this:
%       Adult rural: 
%       Adult semi urban:
%       Adult urban:
%       Children rural: 
%       Children semi urban:
%       Children urban:
%       Dog rural:
%       Dog semi urban:
%       Dog urban:
%       Chicken rural:
%       Chicken semi urban:
%       Chicken urban:
%

%Outputs: 
%   rtn: a vector of length n. Each entry is a sample from the distribution 
%           of the specific host in the specific location 

%% Read the data
load('Host_dist_data/HostDistData.mat')

if(hostType == 0) % Adults
    if(location == 0) % Rural (Marajo)
        data = adult_rural;  %Poiss
        pd = poissfit(data);
        rtn = poissrnd(pd, n, 1);
		
		%Ensure at least one adult per household!
        rtn(rtn == 0) = 1;
    elseif(location == 1) % Semi urban (Campo Grande)
        data = adult_semi_urban;  %Poiss
        pd = poissfit(data);
        rtn = poissrnd(pd, n, 1);
		
		%Ensure at least one adult per household!
        rtn(rtn == 0) = 1;
    else % Urban (Aracatuba)
        data = adult_urban;    %Poiss
        pd = poissfit(data);
        rtn = poissrnd(pd, n, 1);
		
		%Ensure at least one adult per household!
        rtn(rtn == 0) = 1;	
    end
    
elseif(hostType == 1) % Children
    if(location == 0) % Rural (Marajo)
        data = children_rural;  %NB 
        pd = fitdist(data, 'NegativeBinomial');        
        rtn = random(pd, n, 1);
    elseif(location == 1) % Semi urban (Campo Grande)
        data = children_semi_urban;  %NB
        pd = fitdist(data, 'NegativeBinomial');       
        rtn = random(pd, n, 1);
    else % Urban (Aracatuba)
        data = children_urban;    %NB
        pd = fitdist(data, 'NegativeBinomial');
        rtn = random(pd, n, 1);
    end
    
elseif(hostType == 2) % Dogs
    if(location == 0) % Rural (Marajo)
        data = dog_rural;  %NB 
        pd = fitdist(data, 'NegativeBinomial');
        rtn = random(pd, n, 1);
    elseif(location == 1) % Semi urban (Campo Grande) 
        data = dog_semi_urban;  %poiss
        pd = poissfit(data);
        rtn = poissrnd(pd, n, 1);
    else % Urban (Aracatuba)
        data = dog_urban;    %NB 
        pd = fitdist(data, 'NegativeBinomial');
        rtn = random(pd, n, 1);
    end
    
    
else % Chickens
    if(location == 0) % Rural (Marajo)
        data = chicken_rural;  %NB 
        pd = fitdist(data, 'NegativeBinomial');        
        rtn = random(pd, n, 1);
    elseif(location == 1) % Semi urban (Campo Grande)
        data = chicken_semi_urban;  %NB
        pd = fitdist(data, 'NegativeBinomial');       
        rtn = random(pd, n, 1);
    else % Urban (Aracatuba)
        data = chicken_urban;    %NB
        pd = fitdist(data, 'NegativeBinomial');
        rtn = random(pd, n, 1);
    end
end
    
    
end

