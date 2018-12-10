%% This script fits to the rural data and produces plots for the paper and the AIC values of best fits.
clear 

%% Read in the data
load('host_dist_data.mat')

%% Figure 2
figure()
clf

%Adults and adolescents
subplot(2, 2, 1)
hold on
set(gca, 'FontSize', 15)
xlabel('Number of adults and adolescents')
ylabel('Probability')
pd_ad = poissfit(adult_rural);
ppdf_ad = poisspdf(0:max(adult_rural), pd_ad);
plot(0:max(adult_rural), ppdf_ad, 'x', 'LineWidth', 3, 'Color', [0 0.4470 0.7410])
box on
xlim([0 10.5])

%Children
subplot(2, 2, 2) 
hold on
set(gca, 'FontSize', 15)
xlabel('Number of children')
ylabel('Probability')
pd_ch_nb = fitdist(children_rural, 'NegativeBinomial');
ppdf_ch_nb = pdf(pd_ch_nb, 0:max(children_rural));
plot(0:max(children_rural), ppdf_ch_nb, 'x', 'LineWidth', 3, 'Color', [0 0.4470 0.7410])
box on
xlim([0 4.5])

%Dogs
subplot(2, 2, 3)
hold on
set(gca, 'FontSize', 15)
xlabel('Number of dogs')
ylabel('Probability')
pd_d_nb = fitdist(dog_rural, 'NegativeBinomial');
ppdf_d_nb = pdf(pd_d_nb, 0:max(dog_rural));
plot(0:max(dog_rural), ppdf_d_nb, 'x', 'LineWidth', 3, 'Color', [0 0.4470 0.7410])
box on
xlim([0 12.5])

%Chickens
subplot(2, 2, 4)
hold on
set(gca, 'FontSize', 15)
xlabel('Number of chickens')
ylabel('Probability')
pd_ck_nb = fitdist(chicken_rural, 'NegativeBinomial');
ppdf_ck_nb = pdf(pd_ck_nb, 0:max(chicken_rural));
plot(0:max(chicken_rural), ppdf_ck_nb, 'x', 'LineWidth', 3, 'Color', [0 0.4470 0.7410])
box on
xlim([0 50.5])



