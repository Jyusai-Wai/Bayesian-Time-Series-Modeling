function [draws, acceptance_rates] = uc_sv_t_sampler(y, u, y0, u0, spec, nsim, burnin)
% =========================================================================
% CORE FUNCTION: Gibbs Sampler for the UC-SV-t Model
% This function runs the MCMC algorithm and returns the posterior draws.
% =========================================================================

    % --- Setup based on inputs ---
    T = length(y);
    % (All prior setups, initializations from your main script go here)
    % ...

    % --- Storage for MCMC draws ---
    % (Initialize all 'store_*' matrices here)
    % ...

    fprintf('Starting MCMC sampling for T=%d...\n', T);
    for isim = 1:(nsim + burnin)
        % --- Gibbs Sampling Steps ---
        % (The entire MCMC loop from your biUC_t1111.m file goes here)
        % 1. Sample tau_pi
        % 2. Sample tau_u
        % 3. Sample rho_pi
        % 4. Sample rho_u
        % 5. Sample lambda
        % 6. Sample h (stochastic volatility)
        % 7. Sample variance parameters
        % 8. Sample nu (degrees of freedom)
        % ...

        % --- Store draws after burnin ---
        if isim > burnin
            % (Store all parameters into the storage matrices)
        end
    end
    
    % --- Package results into structs for easy access ---
    draws.taupi = store_taupi;
    draws.tauu = store_tauu;
    % ... (and so on for all other stored parameters)
    
    acceptance_rates = countstate / (nsim + burnin);
end
