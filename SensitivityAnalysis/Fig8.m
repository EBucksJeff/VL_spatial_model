%% Code to produce Figure 8

%Purpose:
%Plot sensitivity parameter ranking
%Place in ranking order

%--------------------------------------------------------------------------
clear variables

%--------------------------------------------------------------------------
% IMPORT SENSITIVITY COEFFICIENT DATA
%--------------------------------------------------------------------------
load('normalised_sensitivities.mat')

%--------------------------------------------------------------------------
% ORDER SENSITIVITY COEFFICIENTS IN ASCENDING ORDER
%--------------------------------------------------------------------------
[Ranked_S_normalised,SortIdx] = sort(S_normalised,'ascend');

SandflyParamFlag = SortIdx >=11; %Indicator variable, pick out params associated with sandflies
YVal = 1:1:15; %Y values to be passed to plot function.

%%
%--------------------------------------------------------------------------------------
% PRODUCE FIGURE (Left Axis - Parameter Ranking; Paramter set ID number placed next to data point)
%---------------------------------------------------------------------------------------
fig = figure();
clf
position = [100, 100, 550, 1.25*450];
set(0, 'DefaultFigurePosition', position);
hold on


%Add markers to plot
clf
%Plot dog associated parameter sens. coefficients
p1 = plot(Ranked_S_normalised(SandflyParamFlag == 0), YVal(SandflyParamFlag == 0), 'x', 'MarkerSize', 10, 'LineWidth', 3,'DisplayName','Dog associated only');
hold on 
%Plot sandfly associated parameter sens. coefficients
p2 = plot(Ranked_S_normalised(SandflyParamFlag == 1), YVal(SandflyParamFlag == 1), '.', 'MarkerSize', 30, 'LineWidth', 3,'DisplayName','Sand fly associated');

%List parameter set IDs in reverse ranking order, as a string
YAxisLabels = cell(1,15);
for ii=1:15  
    YAxisLabels{ii} = num2str(SortIdx(ii));
end
StrArray_YAxisLabels = string(YAxisLabels);

%Text label offset parameters
MarkerLabelXOffset = 200;
LabelStartX = Ranked_S_normalised + MarkerLabelXOffset;
LabelStartY = YVal;

% Add text description to data points (number only)
LabelText = strcat(StrArray_YAxisLabels);
text(LabelStartX,LabelStartY,LabelText,'Fontsize',12)

%Set left y-axis properties, parameter ranking numbers
ylabel('Parameter ranking')
yticks(1:1:15)
yticklabels({'15','14','13','12','11','10','9','8','7','6','5','4','3','2','1'})
ylim([0.5 15.5])

%Set x-axis properties
xlabel('Sensitivity coefficient (\Upsilon_{p}^{u})')
xlim([0 6500])

%Set up legend
leg = legend([p1;p2],'Location','Southeast');
title(leg,'Parameter type')

%Set figure properties
set(gca,'Fontsize',16)
set(gca,'LineWidth',1)
box on
