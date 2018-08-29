%% This script fits to rural host count data using Poisson and negative binomial
% It produces figure 4
% And the values for supplementart table 1
clear 

%% Read in the data
load('host_dist_data.mat')

%% Figure 4
figure(1)
clf

%Adults and adolescents: just fit Poisson
subplot(2, 2, 1)
hold on
histogram(adult_rural, 'Normalization', 'probability',...
    'FaceColor', [0.75 0.75 0.75])
set(gca, 'FontSize', 15)
xlabel('Number of adults and adolescents')
ylabel('Probability')
pd_ad = poissfit(adult_rural);
ppdf_ad = poisspdf(0:max(adult_rural), pd_ad);
plot(0:max(adult_rural), ppdf_ad, 'LineWidth', 3, 'Color', [0 0.4470 0.7410])
xlim([0 10.5])
box on

%Children: fit both
subplot(2, 2, 2) 
hold on
histogram(children_rural, 'Normalization', 'probability',...
    'FaceColor', [0.75 0.75 0.75])
set(gca, 'FontSize', 15)
xlabel('Number of children')
ylabel('Probability')
pd_ch_p = poissfit(children_rural);
ppdf_ch_p = poisspdf(0:max(children_rural), pd_ch_p);
plot(0:max(children_rural), ppdf_ch_p, 'LineWidth', 3,...
    'Color', [0 0.4470 0.7410])

pd_ch_nb = fitdist(children_rural, 'NegativeBinomial');
ppdf_ch_nb = pdf(pd_ch_nb, 0:max(children_rural));
plot(0:max(children_rural), ppdf_ch_nb, 'LineWidth', 3,...
    'Linestyle', '-.', 'Color', [0.85 0.3250 0.098])

xlim([-0.5 4.5])
box on

%Dogs: fit both
subplot(2, 2, 3)
hold on
histogram(dog_rural, 'Normalization', 'probability',...
    'FaceColor', [0.75 0.75 0.75])
set(gca, 'FontSize', 15)
xlabel('Number of dogs')
ylabel('Probability')
pd_d_p = poissfit(dog_rural);
ppdf_d_p = poisspdf(0:max(dog_rural), pd_d_p);
plot(0:max(dog_rural), ppdf_d_p, 'LineWidth', 3,...
    'Color', [0 0.4470 0.7410])

pd_d_nb = fitdist(dog_rural, 'NegativeBinomial');
ppdf_d_nb = pdf(pd_d_nb, 0:max(dog_rural));
plot(0:max(dog_rural), ppdf_d_nb, 'LineWidth', 3,...
    'Linestyle', '-.', 'Color', [0.85 0.3250 0.0980])

xlim([-0.5 12.5])
box on

%Chickens: fit both
subplot(2, 2, 4)
hold on
histogram(chicken_rural, 10, 'Normalization', 'probability',...
    'FaceColor', [0.75 0.75 0.75])
set(gca, 'FontSize', 15)
xlabel('Number of chickens')
ylabel('Probability')
pd_ck_p = poissfit(chicken_rural);
ppdf_ck_p = poisspdf(0:max(chicken_rural), pd_ck_p);
plot(0:max(chicken_rural), ppdf_ck_p, 'LineWidth', 3, 'Color', [0 0.4470 0.7410])

pd_ck_nb = fitdist(chicken_rural, 'NegativeBinomial');
ppdf_ck_nb = pdf(pd_ck_nb, 0:max(chicken_rural));
plot(0:max(chicken_rural), ppdf_ck_nb, 'LineWidth', 3,...
    'Linestyle', '-.', 'Color', [0.85 0.3250 0.098])

xlim([-0.5 50.5])
box on

%% AIC computation for which distribution is best fit
% For children, dogs, and chickens. 

% For Poisson fits:
AIC_ad = 2 - 2*sum(log(poisspdf(adult_rural, pd_ad)));
AIC_ch_p = 2 - 2*sum(log(poisspdf(children_rural, pd_ch_p)));
AIC_d_p = 2 - 2*sum(log(poisspdf(dog_rural, pd_d_p)));
AIC_ck_p = 2 - 2*sum(log(poisspdf(chicken_rural, pd_ck_p)));

% For neg binomial
AIC_ch_nb = 2*2 - 2*sum(log(pdf(pd_ch_nb, children_rural)));
AIC_d_nb = 2*2 - 2*sum(log(pdf(pd_d_nb, dog_rural)));
AIC_ck_nb = 2*2 - 2*sum(log(pdf(pd_ck_nb, chicken_rural)));




