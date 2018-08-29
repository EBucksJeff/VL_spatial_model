%% This creates figure S1: seasonality of sandfly abundance per household
clear

%% Load data
load('FigS1_PlotData.mat')

%% Smooth 
for(i=1:8)
    total_by_house{i} = total(site == i);
    date_by_house{i} = date(site == i);
    smoothed{i} = smooth(date_by_house{i}, total_by_house{i}, 0.3, 'lowess');
end

%% Plot 
figure(1)
clf

titles = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'};
x_panel = [(date(95)+date(96))/2 (date(171)+date(172))/2 (date(171)+date(172))/2 (date(95)+date(96))/2]; 
y_panel = [0 0 1200 1200];

for(i=1:8)
    subplot(4, 2, i)
    pl1 = patch(x_panel, y_panel, [0.9 0.9 0.9], 'LineStyle', 'none', 'DisplayName', 'Wet season');
    hold on
    pl2 = plot(date_by_house{i}, total_by_house{i}, '.k',...
        'MarkerSize', 15, 'DisplayName', 'Female sandflies trapped');
    datetick('x', 'mmm')
    ylim([0 1200])
    pl3 = plot(date_by_house{i}, smoothed{i}, 'r', 'LineWidth', 2, 'DisplayName', 'Lowess smoother');
    hold off
    set(gca, 'FontSize', 15)
    title(['House ', titles{i}]) 
    set(gca, 'Layer', 'top') %Bring axis infront of patch object   
    
    if(i==2)
    legend([pl1, pl2, pl3])
    end    
end

