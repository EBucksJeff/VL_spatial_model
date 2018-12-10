%% This creates figure S4: seasonality of sandfly abundance by location
clear

%% Load data
load('FigS2_PlotData.mat')

% uHosts is chicken, dog, house in that order

%% Smooth

smoothed{1} = smooth(uDate, uHosts{1}, 0.3, 'lowess');
smoothed{2} = smooth(uDate, uHosts{2}, 0.3, 'lowess');
smoothed{3} = smooth(uDate, uHosts{3}, 0.3, 'lowess');

%% Plot

figure(1)
clf

titles = {'In the chicken shed', 'At the dog', 'Inside the house'};
x_panel = [(date(95)+date(96))/2 (date(171)+date(172))/2 (date(171)+date(172))/2 (date(95)+date(96))/2]; 
y_panel = [0 0 350 350];

for(i=1:3)
    subplot(3, 1, i)    
    pl1 = patch(x_panel, y_panel, [0.9 0.9 0.9], 'LineStyle', 'none', 'DisplayName', 'Wet season');
    hold on
    pl2 = plot(uDate, uHosts{i}, '.k', 'MarkerSize', 15, 'DisplayName', 'Mean female sand flies trapped');
    pl3 = plot(uDate, smoothed{i}, 'r', 'LineWidth', 2, 'DisplayName', 'Lowess smoother');
    datetick('x', 'mmm')
    ylim([0 350])
    set(gca, 'FontSize', 15)
    title(titles{i})
    set(gca, 'Layer', 'top') %Bring axis infront of patch object 
    if(i==1)
      legend([pl1 pl2 pl3])
    end
end

