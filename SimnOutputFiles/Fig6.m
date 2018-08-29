%% Code to produce Figure 6

clear

%--------------------------------------------------------------------------
%%% LOAD DATA (using baseline biological parameters)
filename = ('VLmodel_BiolParamSensitivitySimn_Run#1.mat'); 
load(filename);

%%
%--------------------------------------------------------------------------
%% PREVALENCE CALCULATION (All four infected classes in numerator)
%--------------------------------------------------------------------------

%Get number of iterations in data sample
run_num = size(total_species_popn,2);

%Get number of timesteps (number of columns in prevalence count arrays)
timestep_num = size(dog_lowinf_prev,2);

%Get number of alive dogs at each timestep
alive_dog_count = dog_sus_prev + dog_exposed_prev + dog_neverinf_prev + dog_lowinf_prev + dog_highinf_prev;

%Calculate normalsied prevalence (proportion of dogs exposed or infected)
normalised_prev = 1 - (dog_sus_prev./alive_dog_count);

%Get 2.5th, 50th, 97.5th percentiles for prevalence (at each time step)
NormalisedPrevSummStat = prctile(normalised_prev,[2.5 50 97.5]);
%   --> outputs a 2d array, row for each requested prctile value, 
%                           col for each timestep

%--------------------------------------------------------------------------
% PRODUCE FIGURE
figure()
clf
position = [100, 100, 550, 450];
set(0, 'DefaultFigurePosition', position);
hold on

%Set column idx to access 
StartIdx = (3*365) + 1;

%Set x axis values
x = 0:1:(5*365)-1;

%Shade prediction interval region grey
%Done by creating 2D polygon
x2 = [x, fliplr(x)]; 
inBetween = [NormalisedPrevSummStat(1,StartIdx:end), fliplr(NormalisedPrevSummStat(3,StartIdx:end))];
f1 = fill(x2, inBetween,[0.5 0.5 0.5],'Linestyle','none','DisplayName','95% prediction interval');

%Plot an example simulation replicate as thin, blue dashed line
p1 = plot(normalised_prev(100,StartIdx:end)',':','Color',[0 0 0.8],'Linewidth',1,'DisplayName','Example individual runs');
plot(normalised_prev(199,StartIdx:end)',':','Color',[0 0 0.8],'Linewidth',1);

%Plot median as thicker, red line
p2 = plot(NormalisedPrevSummStat(2,StartIdx:end),'--','Color',[0.8 0 0],'Linewidth',2.5,'DisplayName','Median'); 

%Add a legend
legend([f1,p2,p1])

%Set x-axis properties
xlabel('Time elapsed (years)');
xticks([0*365 1*365 2*365 3*365 4*365 5*365])
xticklabels({'5','6','7','8','9','10'})
xlim([0 5*365])

%Set y-axis properties
ylabel('VL prevalence (%)');
%ylabel('\it L. infantum\rm prevalence (%)')
yticks([0.45 0.50 0.55 0.60 0.65 0.70 0.75])
yticklabels({'45','50','55','60','65','70','75'})
ylim([0.45 0.76])

%Set figure properties
set(gca,'FontSize',16);
set(gca,'LineWidth',1);
box on
