%% Code to produce Figure 6

%Purpose:
%Plot median & 95% prediction interval of VL prevalence in dogs  
%for various biological parameter sets

%Figure specifics:
%Panel per parameter set tested
%Violin plot per parameter value, plotted in ascending order

%NOTE: Run in MATLAB 2016b

%--------------------------------------------------------------------------
clear variables

%--------------------------------------------------------------------------
% ADD PATH DEPENDENCIES
%--------------------------------------------------------------------------
addpath('Violinplot-Matlab') 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PARAMETER VALUES TESTED     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set up to allow for different number of test values per parameter

BiolParam_to_test = 15; %Number of parameters to test

biol_param = cell(BiolParam_to_test,1);
%Row vectors of values to test for biologic parameters of
%interest. (Again, order matches the list in Parameter_table.docx))
%Note, first value for each biological parameter is the baseline value
biol_param{1} = [0.3,0.02,0.7,2];
biol_param{2} = [0.55,0.14,0.28,0.42];
biol_param{3} = [0.37,0.25,0.6,0.8];
biol_param{4} = [0.13,0.0064,0.29,0.43]; 
biol_param{5} = [182.5,153,212,240];
biol_param{6} = [0.046,0.036,0.070,0.095]; 
biol_param{7} = [0.060,0.036,0.080,0.095]; 
biol_param{8} = [0.064,0.036,0.080,0.095]; 
biol_param{9} = [0.038,0.032,0.034,0.036];
biol_param{10} = [121,0,243,578]; 
biol_param{11} = [1/3,0.25,0.4,0.5]; 
biol_param{12} = [0.01,0.002,0.1,0.26];
biol_param{13} = [0.321,0.1,0.2,0.5]; %VALUES TO BE CONFIRMED!
biol_param{14} = [0.275,0.023,0.15,0.45];
biol_param{15} = [0.9,0.75,0.8,0.85];

%List name of paramters to be used as x-axis labels in figures
BiolParam_xlabels = cell(BiolParam_to_test,1);

BiolParam_xlabels{1} = 'Dog range of effect (km)';
BiolParam_xlabels{2} = 'Propn. of dogs never infectious';
BiolParam_xlabels{3} = 'Propn. of infectious dogs that are highly infectious';
BiolParam_xlabels{4} = 'Prob. of newly introduced dog being infected';
BiolParam_xlabels{5} = 'Average latent period (days)';
BiolParam_xlabels{6} = 'Average mortality rate of never infectious dog (per month)';
BiolParam_xlabels{7} = 'Average mortality rate of low infectious dog (per month)';
BiolParam_xlabels{8} = 'Average mortality rate of highly infectious dog (per month)';
BiolParam_xlabels{9} = 'Average mortality rate of susceptible dog (per month)';
BiolParam_xlabels{10} = 'Replacement time (days)';
BiolParam_xlabels{11} = 'Sand fly bite rate (per day)';
BiolParam_xlabels{12} = 'Background propn. of sand flies infected';
BiolParam_xlabels{13} = 'Prob. suscept. dog infected when bitten by infected sand fly';
BiolParam_xlabels{14} = 'Prob. suscept. sand fly infected when biting infectious dog';
BiolParam_xlabels{15} = 'Propn. of female sand flies unobserved';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Specify number of households used
n = 200;

%Specify number of runs used
run_num = 1000;

%Specify number of biological parameters varied
BiolParam = 15;

%%% Import data files for number of infection events - BASELINE BIOLOGICAL PARAM. VALUES %%%
%NOTE: This will be used in all plots
%Brazil village
filename = ('../SimnOutputFiles/VLmodel_BiolParamSensitivitySimn_Run#1.mat'); 
load(filename); 

% work out average prevalence in last year
prev_indices = max(1, size(dog_highinf_prev, 2) - 364):size(dog_highinf_prev, 2); % indices of interest
prev_baseline = mean((dog_highinf_prev(:, prev_indices) + dog_lowinf_prev(:, prev_indices) + dog_neverinf_prev(:, prev_indices) + dog_exposed_prev(:, prev_indices))...
    ./(dog_highinf_prev(:, prev_indices) + dog_lowinf_prev(:, prev_indices) + dog_neverinf_prev(:, prev_indices)...
    + dog_sus_prev(:, prev_indices) + dog_exposed_prev(:, prev_indices)), 2);



%--------------------------------------------------------------------------
%Iterate through each run 

% Initialise cell to store prevalence quantities
prev_BiolSens = cell(BiolParam,1);

simn_idx = 2; %Set up simulation number counter, used to load files 
for i = 1:BiolParam
    non_baseline_values = numel(biol_param{i}) - 1; %Get number of simulations carried out for current biological parameter of interest   

    % prevalence
    prev_BiolSens{i,1} = zeros(non_baseline_values,run_num);

%     % incidence
%     incidence_BiolSens{i,1} = zeros(non_baseline_values,run_num);

    j = 1;
    while j<= non_baseline_values
        %Load file, number simn_idx
        filename = (['../SimnOutputFiles/VLmodel_BiolParamSensitivitySimn_Run#' num2str(simn_idx) '.mat']);
        load(filename);

        % prevalence
        prev_BiolSens_temp = mean((dog_highinf_prev(:, prev_indices) + dog_lowinf_prev(:, prev_indices) + dog_neverinf_prev(:, prev_indices) + dog_exposed_prev(:, prev_indices))...
            ./(dog_highinf_prev(:, prev_indices) + dog_lowinf_prev(:, prev_indices) + dog_neverinf_prev(:, prev_indices)...
            + dog_sus_prev(:, prev_indices) + dog_exposed_prev(:, prev_indices)), 2); % average prevalence
        prev_BiolSens{i,1}(j,:) = prev_BiolSens_temp;

%         
        j = j + 1;
        simn_idx  = simn_idx + 1;
    end
end

%%
%--------------------------------------------------------------------------
%%% Amalgamate baseline values and non_basline_value outputs into single arrays

%Initialise storage cells
Combined_prev_BiolSens = cell(BiolParam,1);


%Iterate through each biological parameter set
%Assign baseline value related output to first row
%Non-baseline related value outputs to subsequent rows
for ii = 1:BiolParam
    %Prevalence
    Combined_prev_BiolSens{ii,1} = [prev_baseline'; prev_BiolSens{ii}];

end


%--------------------------------------------------------------------------
% Iterate through each parameter set, reorder tested values & outputs so in
% ascending parameter value order 

%Initialise storage cells for ordered biological parameter values and array 
%index locations
biol_param_AscendOrder = cell(BiolParam,1);
BiolParam_AscendOrderIdx = cell(BiolParam,1);

%Initialise storage cells fo arrays combining baseline and non-baseline
%outputs
Combined_prev_BiolSens_AscendOrder = cell(BiolParam,1);


for ii = 1:BiolParam

    %Reorder biol_param values and simulation_output
    %biol_param into ascending order (from simulation index order). 
    %Column Indx of simulation_output match ascending order paramter values
    [biol_param_AscendOrder{ii,1},BiolParam_AscendOrderIdx{ii,1}] = sort(biol_param{ii});
    
    Combined_prev_BiolSens_AscendOrder{ii,1} = Combined_prev_BiolSens{ii,1}(BiolParam_AscendOrderIdx{ii,1},:);
    
end


%For parameter set 11, alter tested values so 1/3 is printed on graph as
%0.333
biol_param_AscendOrder{11} = [0.25,0.333,0.4,0.5]; 


%% Construct one overall plot with a subplot for each parameter tested (prevalence)
scrsz = get(0, 'ScreenSize');

figure('Color',[1 1 1]);
set(gcf, 'Position', [0 0 round(scrsz(3)) round(scrsz(4))]);
for ii = 1:BiolParam
    subplot(4,4,ii)
    hold on
    
    %Transpose array of avg. prevalence values per simulation
    %Now column per parameter value tested, Row per simulation replicate
    PrevDataTransposed = Combined_prev_BiolSens_AscendOrder{ii}'; 
    violins = violinplot(PrevDataTransposed, biol_param_AscendOrder{ii},'ShowData',false);
    
    TestParamNum = numel(biol_param{ii}); %Get number of simulations carried out for current biological parameter of interest
    jj = 1;
    while jj <= TestParamNum         
        %If baseline parameter value, plot violins in black, with red symbol for median
        %Otherwise, plot violins in blue, with black symbol denoting median
        if BiolParam_AscendOrderIdx{ii,1}(jj) == 1
            
            %Set violin plot region properties
            violins(jj).ViolinColor = [0.5 0.5 0.5];
            violins(jj).ViolinAlpha = 1.0; %shading transparency
            violins(jj).EdgeColor = [0.8 0 0];
            
            %Set median marker properties
            violins(jj).MedianColor = [0.8 0 0];
            violins(jj).MedianPlot.Marker = 's';
            violins(jj).MedianPlot.MarkerEdgeColor = [0.8 0 0];
            violins(jj).MedianPlot.LineWidth = 1;
            
            %Set whisker line properties
            violins(jj).BoxColor = [0.8 0 0];
        else
            %Set violin plot region properties
            violins(jj).ViolinColor = [0 0 0.8];
            violins(jj).ViolinAlpha = 0.5; %shading transparency
            violins(jj).EdgeColor = [0 0 0];
            
            %Set median marker properties
            violins(jj).MedianPlot.MarkerFaceColor = [1 1 1];
            violins(jj).MedianPlot.MarkerEdgeColor = [0 0 0];
            violins(jj).MedianPlot.LineWidth = 1;
            
            %Set whisker line properties
            violins(jj).BoxColor = [0 0 0];
        end             
        jj = jj + 1;
    end
    
    FirstColPlotIdx = [1 5 9 13]; %Pick out plots that will be given y-axis label  
    if sum(ii == FirstColPlotIdx) == 1
        ylabel 'Average prevalence (%)'
    end
    
    %Get xlabel from cell array
    xlabel(BiolParam_xlabels{ii}); 
    
    %Set y-axis properties
    yticks([0.3 0.5 0.7 0.9])
    yticklabels({'30','50','70','90'})
    ylim([0.22 0.98])
    
    ax = gca;
    set(gca, 'FontSize', 10);
    
    %Add ID number to plot
    x1 = -0.86;
    y1 = 1.05;
    txt1 = ['(', num2str(ii) ,')'];
    text(x1,y1,txt1,'FontWeight','bold','FontSize',12)
    
    %Add surrounding box
    box on
    
    hold off
end
