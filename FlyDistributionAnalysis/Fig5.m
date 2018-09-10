%% Code to produce Figure 5: fly seasonality 

clear

%Read in the data
load('Fig5_PlotData.mat');

for i = 1:numel(uDate)-1
   uTotal(i) = mean(SandflySeasonalityData(indx(i):indx(i+1)-1));

end
uTotal(numel(uDate)) = mean(SandflySeasonalityData(indx(end):numel(SandflySeasonalityData)));

%% Plot points with smoothing lines

figure()
clf

%Smooth the whole thing
smoothAll = smooth(uDate, uTotal, 0.35, 'lowess');

%Interpolate  the smoothed function to get value for every day 
allDate = uDate(1):uDate(1)+365;
uDateAdj = [uDate; uDate(1)+365];
smoothAllAdj = [smoothAll; smoothAll(1)];
seasonalityVec = interp1(uDateAdj, smoothAllAdj, allDate);

%Set up background shading
xSeason = [734898 734989 734989 734898];
ySeason_2 = [0 0 1200 1200];
p2 = patch(xSeason, ySeason_2, [0.9 0.9 0.9],...
    'DisplayName', 'Wet season', 'LineStyle', 'none');
hold on

%Plot number of sandflies in each household
t2 = plot(date, SandflySeasonalityData, 'k.', 'MarkerSize', 15,...
    'DisplayName', 'Female sand flies trapped');

%Plot mean over trapping sites;
m2 = plot(uDate, uTotal,'--','Color',[0 0 0.8], 'LineWidth', 2,...
    'DisplayName', 'Mean over all sites');

%Plot smoother curve
l2 = plot(allDate, seasonalityVec, 'r', 'Linewidth', 3,'DisplayName','Lowess smoother');

datetick('x', 'mmm')
xlim([date(1) date(1)+365])
set(gca, 'FontSize', 22)
xlabel('Date')
ylabel('Number of sand flies')
ylim([0 1200])
legend([p2, t2, m2, l2]);
set(gca, 'Layer', 'top') %Bring axis infront of patch object 

 

