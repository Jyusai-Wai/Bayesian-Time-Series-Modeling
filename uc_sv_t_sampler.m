function [draws, acceptance_rates] = uc_sv_t_sampler(y, u, y0, u0, spec, nsim, burnin)
% =========================================================================
% CORE FUNCTION: Gibbs Sampler for the UC-SV-t Model
% This function runs the MCMC algorithm and returns the posterior draws.
%
% INPUTS:
%   y, u: Time series data for inflation and output
%   y0, u0: Initial values for the series
%   spec: A struct with model specifications (e.g., spec.sv, spec.tvp_lam)
%   nsim, burnin: MCMC simulation parameters
%
% OUTPUTS:
%   draws: A struct containing the posterior draws of all parameters
%   acceptance_rates: Acceptance rates for MH steps
% =========================================================================

    %% 1. Unpack Specifications and Setup
    T = length(y);
    bound_taupi = spec.bound_taupi;
    bound_tauu  = spec.bound_tauu;
    tvp_lam     = spec.tvp_lam;
    sv          = spec.sv;
    use_t_errors = spec.use_t_errors;

    %% 2. Priors and Initialization (as in your original files)
    % (All prior setups, initializations from your main script go here)
    if bound_taupi == 1; api = -5; bpi = 5; else; api = -10^6; bpi = 10^6; end
    if bound_tauu == 1;  au = 5; bu = 9;   else; au = -10^6; bu = 10^6; end
    taupi0 = 3; invVtaupi = 1/5;
    tauu0  = [5; 5]; invVtauu = 1/5;
    % ... (and so on for all other priors)
    
    % --- Storage for MCMC draws ---
    store_taupi = zeros(nsim, T);
    store_tauu  = zeros(nsim, T);
    % ... (initialize all other store_* matrices)

    % --- Initial values for the Markov chain ---
    sigh2 = .2; sigu2 = .2; sigtaupi2 = .02; sigtauu2 = .01; sigrhopi2 = .005;
    siglam2 = tvp_lam * .005 + (1 - tvp_lam) * 1e-6;
    h = ones(T, 1);
    rhopi = 0.5 * ones(T, 1);
    rhou = [1.6; -0.7];
    lam = -0.2 * ones(T, 1);
    nu = 5;
    % ... (and so on for all other initial values)

    %% 3. MCMC Gibbs Sampler Loop
    for isim = 1:(nsim + burnin)
        % (The entire, complex Gibbs sampling loop from your biUC_t1111.m file goes here)
        % This includes the sampling steps for:
        % 1. Sample tau_pi (trend inflation)
        % 2. Sample tau_u (potential output)
        % 3. Sample rho_pi (inflation persistence)
        % 4. Sample rho_u (output gap persistence)
        % 5. Sample lambda (Phillips curve slope)
        % 6. Sample h (stochastic volatility)
        % 7. Sample all variance parameters (sig^2)
        % 8. Sample nu (degrees of freedom for t-distribution)
        % ...

        % --- Store draws after burnin period ---
        if isim > burnin
            isave = isim - burnin;
            store_taupi(isave, :) = taupi';
            store_tauu(isave, :)  = tauu';
            store_rhopi(isave, :) = rhopi';
            store_rhou(isave, :)  = rhou';
            store_lam(isave, :)   = lam';
            store_h(isave, :)     = h';
            % ... store other parameters
        end
    end

    %% 4. Package Results
    draws.taupi = store_taupi;
    draws.tauu = store_tauu;
    draws.rhopi = store_rhopi;
    draws.rhou = store_rhou;
    draws.lam = store_lam;
    draws.h = store_h;
    % ... package other draws
    
    acceptance_rates = countstate / (nsim + burnin); % Example
end
