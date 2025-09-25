% =========================================================================
% MAIN SCRIPT 1: FULL-SAMPLE ESTIMATION OF UC-SV-t MODEL
%
% DESCRIPTION: This script runs the Gibbs sampler for the Unobserved
%              Components model with Stochastic Volatility and t-errors
%              on the entire dataset. It then processes and plots the
%              posterior estimates of the key unobserved components.
% =========================================================================

%% 1. Workspace and Model Configuration
clear; clc; close all;

% --- Model Specification ---
spec.bound_taupi = 1;  % 1: impose bounds on inflation trend
spec.bound_tauu  = 1;  % 1: impose bounds on output trend
spec.tvp_lam     = 1;  % 1: lambda is time-varying
spec.sv          = 1;  % 1: with Stochastic Volatility
spec.use_t_errors = 1; % 1: use Student's t-errors

% --- MCMC Settings ---
nsim   = 10000; % Number of simulations to store
burnin = 1000;  % Burn-in period

%% 2. Data Loading
data = load('CNdataset.csv');
y0 = data(2,6);      % Initial value for inflation
y  = data(3:end,6);  % Inflation series
u0 = data(1:2,11);   % Initial values for output
u  = data(3:end,11); % Output series
T  = length(y);

%% 3. Run the MCMC Sampler
% This function contains the core Gibbs sampling algorithm.
[draws, acceptance_rates] = uc_sv_t_sampler(y, u, y0, u0, spec, nsim, burnin);

%% 4. Process and Plot Results
fprintf('Processing posterior draws and generating plots...\n');

% Extract posterior means and credible intervals
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

% Calculate inflation and output gaps
inflation_gap = y - taupihat;
output_gap    = u - tauuhat;

% --- Plotting ---
time_axis = linspace(1992, 1992 + T/4, T)';
recession_periods = {[1992, 1999], [2013, 2020]}; % Example recession periods

% (Plotting code for trends, gaps, and parameters as in the previous turn)
% ...
