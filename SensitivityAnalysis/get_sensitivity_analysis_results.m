%Use this script to get results for the sensitivity analysis 

%NOTE: Run in MATLAB 2016b

%The function "normalised_sensitivity_analysis" takes as inputs 
%a vector of parameter values and the number (1 - 15) of the parameter being tested. 
%It outputs the sensitivity as described in the Damiani paper. 
%We can then rank the parameters according to these sensitivities. 

clear 

%Parameter values we consider
biol_param{1} = [0.3,0.02,0.7,2];
biol_param{2} = [0.55,0.14,0.28,0.42];
biol_param{3} = [0.37,0.25,0.6,0.8];
biol_param{4} = [0.13,0.0064,0.29,0.43]; 
biol_param{5} = [182.5,153,212,240];
biol_param{6} = [0.046,0.036,0.070,0.095]; 
biol_param{7} = [0.060,0.036,0.080,0.095]; 
biol_param{8} = [0.064,0.036,0.080,0.095]; 
biol_param{9} = [0.038,0.0083,0.020,0.030];
biol_param{10} = [121,0,243,578]; 
biol_param{11} = [1/3,0.25,0.4,0.5]; 
biol_param{12} = [0.01,0.002,0.1,0.26];
biol_param{13} = [0.321,0.1,0.2,0.5];
biol_param{14} = [0.275,0.023,0.15,0.45];
biol_param{15} = [0.9,0.75,0.8,0.85];

S_normalised = zeros(1, 15);
%To obtain normalised results
for(i=1:15)
    S_normalised(i) = normalised_sensitivity_analysis(biol_param{i}, i);
    display(['Parameter set ' num2str(i) ' of 15 done (normalised).'])
end

%Save S_normalised for plotting