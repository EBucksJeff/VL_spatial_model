%% Code based on method described in: Damiani et al

function [S] = normalised_sensitivity_analysis(param_values_SimOrder, which_param)
%When we consider a parameter value of 0 we cannot compute the logarithm. 
%Affected variable has unit values of days. Therefore, changing 0 values to 1.
param_values_SimOrder(param_values_SimOrder == 0) = 1;

%%
%%%%%%% Step 1: estimation of the probability density functions
%
% Need to estimate one pdf of the output per parameter
%
% Do not impose any analytical form on the distribution, just estimate it
% from the output of the model simulated n times. Smooth the histrogram of
% the simulated values of the output via a kernal density estimation
% process.
%
% They specifically used the Matlab Statistics Toolbox function ksdensity
% and a window parameter that is a function of the number of sampled
% points.
%
% The axis of the pdfs should span from the minimum to the maximum of the
% output observed across all cases
%
% Output: average prevalence over the final year
%

% Obtain average prevalence for the final year per simulation for the
% output of choice
n_values=length(param_values_SimOrder);

for(j=1:n_values)
    simulation_output_SimOrder(:, j) = return_histogram_values_prev(which_param, j);
end

%Reorder param_values and simulation_output
%Param_values into ascending order (from simulation index order). 
%Column Indx of simulation_output match ascending order paramter values
[param_values_AscendOrder,AscendOrderIdx] = sort(param_values_SimOrder);
simulation_output = simulation_output_SimOrder(:,AscendOrderIdx);

pdf_max = max(max(simulation_output));
pdf_min = min(min(simulation_output));
n_points = 100;
indx_points = linspace(pdf_min, pdf_max, n_points);
pdf_estimates = zeros(n_points, n_values);

for(j=1:n_values)

    % Apply ksdensity to estimate the pdf
    [pdf] = ksdensity(simulation_output(:,j), indx_points);
    
    % Plot pdf estimate on top of histrogram
    %hold on
    %plot(indx_points, pdf, 'lineWidth', 3)
    %hold off
    
    pdf_estimates(:,j) = pdf;
    
end

%hold off

%%
%%%%%%%%%% Step 2: Estimation of the derivative of density functions with
%%%%%%%%%% respect to the parameter
%
% Need to analyse the dependence of the (estimated) pdfs on the parameter of
% interest.
%
% Do this by estimating for each possible value of the output the derivative
% of the pdf by evaluating its slope at each of the points correponding to
% the different values of the parameter that have been tested.
%
% They used the python function matplotlib.mlab.slopes
%
% matplotlib.mlab.slopes(x, y) estimates the slope y'(x) using the slope
% obtained from a parabola through any three consecutive points
%
% Recoded into Matlab as slope.m

derivatives = zeros(n_points, n_values); %for storage

for(j=1:n_points) %x = values of parameter y=values of pdf
    derivatives(j,:) = slope(log(param_values_AscendOrder), pdf_estimates(j,:))';
end


%%
%%%%%%%% Step 3: Determination of a sensitivity coefficient
%
% 3a) Equation (5) from the paper: integrate the abs value of the drivative
% multiplied by the pdf over the output variable. We can do this numerically
% with a sum (a sum over the values of indx_points as this is the range
% of output variables).
%
% 3b) Obtain a single sensitivity coefficient for the parameter by
% integrating over the previously obtained sensitivity curve.
%
% We can do this numerically by approximating each region under the curve
% between two values of the parameter as a trapezoid and calculating the
% area and summing.

S_p = zeros(1, n_values);
for(j=1:n_values)
    S_p(j) = sum(abs(derivatives(:,j)).*pdf_estimates(:,j));
end

S = trapz(log(param_values_AscendOrder), S_p);

end
