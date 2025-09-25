% =========================================================================
% MAIN SCRIPT 2: RECURSIVE FORECASTING EXERCISE
%
% DESCRIPTION: This script performs a recursive out-of-sample forecasting
%              exercise using the UC-SV-t model. It re-estimates the model
%              over an expanding window and generates forecasts for multiple
%              horizons. Finally, it computes RMSFEs and log-predictive
%              likelihoods.
% =========================================================================

clear; clc; close all;

%% 1. Forecasting Configuration
% --- Model Selection ---
% 1: Full UC-SV-t | 2: const-lambda | 3: const-rhopi | 4: const-lambda-rhopi
model_choice = 1;

% --- MCMC Settings for each recursive step ---
nsim   = 10000;
burnin = 1000;

% --- Data Loading ---
data = load('CNdataset.csv');
y_full = data(3:end,6);
u_full = data(3:end,11);
y0 = data(2,6);
u0 = data(1:2,11);
T_full = length(y_full);

% --- Recursive Forecasting Settings ---
T0 = 121; % Start forecasting from this time point (e.g., 2022Q1)
horizons = [1, 4, 8, 12, 16]; % Forecast horizons in quarters
max_h = max(horizons);

% --- Storage for Forecasts ---
y_forecasts = NaN(T_full - T0, length(horizons));
u_forecasts = NaN(T_full - T0, length(horizons));
log_pred_likes = NaN(T_full - T0, length(horizons));

%% 2. Recursive Forecasting Loop
fprintf('Starting recursive forecasting exercise...\n');
start_time = clock;

for t = T0:(T_full - max_h)
    fprintf('Forecasting from T = %d...\n', t);
    
    % --- Create subset of data for this iteration ---
    y_t = y_full(1:t);
    u_t = u_full(1:t);
    
    % --- Set model specification for this run ---
    spec.sv = 1; spec.bound_taupi = 1; spec.bound_tauu = 1; spec.use_t_errors = 1;
    switch model_choice
        case 1; spec.tvp_lam = 1; spec.tvp_rhopi = 1;
        case 2; spec.tvp_lam = 0; spec.tvp_rhopi = 1;
        case 3; spec.tvp_lam = 1; spec.tvp_rhopi = 0;
        case 4; spec.tvp_lam = 0; spec.tvp_rhopi = 0;
    end
    
    % --- Run sampler on current data subset ---
    [draws, ~] = uc_sv_t_sampler(y_t, u_t, y0, u0, spec, nsim, burnin);
    
    % --- Generate forecasts for all required horizons ---
    [y_preds, u_preds, log_likes] = generate_forecasts(draws, y_t, u_t, y0, u0, max_h);
    
    % --- Store the forecasts ---
    forecast_idx = t - T0 + 1;
    for i_h = 1:length(horizons)
        h = horizons(i_h);
        y_forecasts(forecast_idx, i_h) = mean(y_preds(:, h));
        u_forecasts(forecast_idx, i_h) = mean(u_preds(:, h));
        % (Log predictive likelihood calculation would be more complex)
    end
end

fprintf('Forecasting exercise finished in %.2f seconds.\n', etime(clock, start_time));

%% 3. Evaluate and Plot Forecasts
% (Code to calculate RMSFE and plot actual vs. forecast as in your original file)
% ...
