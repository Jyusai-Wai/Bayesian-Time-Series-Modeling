% =========================================================================
% MAIN SCRIPT 1: FULL-SAMPLE ESTIMATION OF UC-SV-t MODEL
%
% DESCRIPTION: This script runs the Gibbs sampler for the Unobserved
%              Components model with Stochastic Volatility and t-errors
%              on the entire dataset. It then processes and plots the
%              posterior estimates of the key unobserved components.
% =========================================================================

%% 1. Workspace and Model Configuration
% -------------------------------------------------------------------------
clear; clc; close all;

% --- MODEL SPECIFICATION ---
spec.bound_taupi = 1;  % 1: impose bounds on inflation trend; 0: no bounds
spec.bound_tauu  = 1;  % 1: impose bounds on output trend; 0: no bounds
spec.tvp_lam     = 1;  % 1: lambda (Phillips curve slope) is time-varying; 0: constant
spec.sv          = 1;  % 1: with Stochastic Volatility; 0: constant variance
spec.use_t_errors = 1; % 1: use Student's t-errors for robustness; 0: use Normal errors

% --- MCMC SETTINGS ---
nsim   = 10000; % Number of simulations to store
burnin = 1000;  % Burn-in period

%% 2. Data Loading and Preparation
% -------------------------------------------------------------------------
disp('Loading and preparing data...');
data = load('CNdataset.csv');

y0 = data(2,6);      % Initial value for inflation
y  = data(3:end,6);  % Inflation series (e.g., CPI growth)
u0 = data(1:2,11);   % Initial values for output/unemployment
u  = data(3:end,11); % Output/unemployment series
T  = length(y);

%% 3. Run the MCMC Sampler
% This function contains the core Gibbs sampling algorithm.
[draws, acceptance_rates] = uc_sv_t_sampler(y, u, y0, u0, spec, nsim, burnin);

%% 4. Process and Plot Results
% -------------------------------------------------------------------------
fprintf('Processing posterior draws and generating plots...\n');

% Extract posterior means and 90% credible intervals (5th and 95th percentiles)
taupihat = mean(draws.taupi)';
taupilb  = quantile(draws.taupi, .05)';
taupiub  = quantile(draws.taupi, .95)';

tauuhat = mean(draws.tauu)';
tauulb  = quantile(draws.tauu, .05)';
tauuub  = quantile(draws.tauu, .95)';

rhopihat = mean(draws.rhopi)';
rhopilb  = quantile(draws.rhopi, .05)';
rhopiub  = quantile(draws.rhopi, .95)';

lamhat = mean(draws.lam)';
lamlb  = quantile(draws.lam, .05)';
lamub  = quantile(draws.lam, .95)';

hhat = mean(exp(draws.h/2))';
hlb  = quantile(exp(draws.h/2), .05)';
hub  = quantile(exp(draws.h/2), .95)';

% Calculate inflation and output gaps based on posterior means
inflation_gap = y - taupihat;
output_gap    = u - tauuhat;

% --- Plotting ---
time_axis = linspace(1992, 1992 + (T-1)/4, T)'; % Assuming quarterly data
recession_periods = {[1992, 1999], [2013, 2020]}; % Example recession periods for shading

% Plot for Trend Inflation and Inflation Gap
figure('Name', 'Inflation Dynamics', 'color', 'w', 'Position', [100, 100, 800, 300]);
subplot(1,2,1);
hold on;
for p = recession_periods; fill([p{1}(1) p{1}(1) p{1}(2) p{1}(2)], [-2 3 3 -2], [0.9 0.9 0.9], 'EdgeColor', 'none'); end
plot(time_axis, taupihat, 'LineWidth', 2, 'Color', 'black');
plot(time_axis, [taupilb taupiub], ':', 'LineWidth', 1.5, 'Color', 'red');
hold off;
xlim([1992 2024]);
title('Trend Inflation (\tau_\pi)'); box off;

subplot(1,2,2);
hold on;
for p = recession_periods; fill([p{1}(1) p{1}(1) p{1}(2) p{1}(2)], [-2 2 2 -2], [0.9 0.9 0.9], 'EdgeColor', 'none'); end
plot(time_axis, inflation_gap, 'LineWidth', 2, 'Color', 'black');
refline(0,0);
hold off;
xlim([1992 2024]);
title('Inflation Gap'); box off;

% Plot for Potential Output and Output Gap
figure('Name', 'Output Dynamics', 'color', 'w', 'Position', [100, 450, 800, 300]);
subplot(1,2,1);
hold on;
for p = recession_periods; fill([p{1}(1) p{1}(1) p{1}(2) p{1}(2)], [5 10 10 5], [0.9 0.9 0.9], 'EdgeColor', 'none'); end
plot(time_axis, tauuhat, 'LineWidth', 2, 'Color', 'black');
plot(time_axis, [tauulb tauuub], ':', 'LineWidth', 1.5, 'Color', 'red');
hold off;
xlim([1992 2024]);
title('Potential Output (\tau_y)'); box off;

subplot(1,2,2);
hold on;
for p = recession_periods; fill([p{1}(1) p{1}(1) p{1}(2) p{1}(2)], [-15 15 15 -15], [0.9 0.9 0.9], 'EdgeColor', 'none'); end
plot(time_axis, output_gap, 'LineWidth', 2, 'Color', 'black');
refline(0,0);
hold off;
xlim([1992 2024]);
title('Output Gap'); box off;

fprintf('Full-sample estimation and plotting complete.\n');
