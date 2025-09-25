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
% -------------------------------------------------------------------------
% --- Model Selection ---
% 1: Full UC-SV-t model (time-varying lambda and rho)
% 2: Constant-lambda model
% 3: Constant-rho_pi model
% 4: Constant-lambda and constant-rho_pi model
model_choice = 1;

% --- MCMC Settings for each recursive step ---
nsim   = 10000;
burnin = 1000;

% --- Data Loading ---
data = load('CNdataset.csv');
y_full = data(3:end, 6);
u_full = data(3:end, 11);
y0 = data(2, 6);
u0 = data(1:2, 11);
T_full = length(y_full);

% --- Recursive Forecasting Settings ---
T0 = 121; % Start forecasting from this time point (e.g., 2022Q1)
horizons = [1, 4, 8, 12, 16]; % Forecast horizons in quarters
max_h = max(horizons);

% --- Storage for Forecasts ---
y_forecasts = NaN(T_full - T0, length(horizons));
u_forecasts = NaN(T_full - T0, length(horizons));
y_actuals = NaN(T_full - T0, length(horizons));
u_actuals = NaN(T_full - T0, length(horizons));

%% 2. Recursive Forecasting Loop
% -------------------------------------------------------------------------
fprintf('Starting recursive forecasting exercise...\n');
start_time = clock;

for t = T0:(T_full - max_h)
    fprintf('Estimating and forecasting from T = %d...\n', t);

    % Create subset of data for this iteration
    y_t = y_full(1:t);
    u_t = u_full(1:t);

    % Set model specification for this run
    spec.sv = 1; spec.bound_taupi = 1; spec.bound_tauu = 1; spec.use_t_errors = 1;
    switch model_choice
        case 1; spec.tvp_lam = 1; spec.tvp_rhopi = 1;
        case 2; spec.tvp_lam = 0; spec.tvp_rhopi = 1;
        case 3; spec.tvp_lam = 1; spec.tvp_rhopi = 0;
        case 4; spec.tvp_lam = 0; spec.tvp_rhopi = 0;
    end

    % Run sampler on current data subset
    [draws, ~] = uc_sv_t_sampler(y_t, u_t, y0, u0, spec, nsim, burnin);

    % Generate forecasts for all required horizons
    [y_preds, u_preds] = generate_forecasts(draws, y_t, u_t, y0, u0, max_h, spec);

    % Store the mean of the predictive distribution and actual values
    forecast_idx = t - T0 + 1;
    for i_h = 1:length(horizons)
        h = horizons(i_h);
        y_forecasts(forecast_idx, i_h) = mean(y_preds(:, h));
        u_forecasts(forecast_idx, i_h) = mean(u_preds(:, h));
        y_actuals(forecast_idx, i_h) = y_full(t + h);
        u_actuals(forecast_idx, i_h) = u_full(t + h);
    end
end

fprintf('Forecasting exercise finished in %.2f seconds.\n', etime(clock, start_time));

%% 3. Evaluate and Plot Forecasts
% -------------------------------------------------------------------------
% Calculate Root Mean Squared Forecast Errors (RMSFE)
RMSFE_y = zeros(length(horizons), 1);
RMSFE_u = zeros(length(horizons), 1);

for i_h = 1:length(horizons)
    h = horizons(i_h);
    forecast_len = T_full - max_h - T0 + 1;
    errors_y = y_actuals(1:forecast_len, i_h) - y_forecasts(1:forecast_len, i_h);
    errors_u = u_actuals(1:forecast_len, i_h) - u_forecasts(1:forecast_len, i_h);
    RMSFE_y(i_h) = sqrt(mean(errors_y.^2));
    RMSFE_u(i_h) = sqrt(mean(errors_u.^2));
end

% Display RMSFE table
fprintf('\n--- Root Mean Squared Forecast Errors ---\n');
fprintf('Horizon  |   Inflation (y) |   Output (u)\n');
fprintf('-------------------------------------------\n');
for i_h = 1:length(horizons)
    fprintf('%2d-Q ahead |      %.4f     |     %.4f\n', horizons(i_h), RMSFE_y(i_h), RMSFE_u(i_h));
end

% Plot actual vs. forecasted values for inflation
figure('Name', 'Inflation Forecasts', 'color', 'w');
for i_h = 1:length(horizons)
    subplot(3, 2, i_h);
    hold on;
    plot(y_actuals(:, i_h), 'k-', 'LineWidth', 1.5);
    plot(y_forecasts(:, i_h), 'r--', 'LineWidth', 1.5);
    hold off;
    title(sprintf('%d-Quarter Ahead Inflation Forecast', horizons(i_h)));
    legend('Actual', 'Forecast');
    box off;
end
