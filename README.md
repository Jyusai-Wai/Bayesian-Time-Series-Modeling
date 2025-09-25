# Bivariate Unobserved Components Model with Stochastic Volatility (UC-SV-t) in MATLAB

This repository contains a MATLAB implementation of a Bivariate Unobserved Components (UC) model with Stochastic Volatility (SV) and Student's t-distributed errors. The model is designed for advanced macroeconomic analysis, particularly for estimating time-varying trend inflation, potential output (or NAIRU), and the Phillips curve.

The estimation is performed using Bayesian methods, specifically a Markov Chain Monte Carlo (MCMC) Gibbs sampler, based on the methodology proposed by Chan, Koop, and Potter (2016). The project includes scripts for both full-sample estimation and a comprehensive out-of-sample recursive forecasting exercise.

## Key Features

-   **Time-Varying Trends**: Both the inflation trend (`tau_pi`) and the potential output trend (`tau_y`) are modeled as random walks.
-   **Time-Varying Phillips Curve**: The slope of the Phillips curve (`lambda`) can be specified as time-varying.
-   **Stochastic Volatility (SV)**: The variance of the inflation shock is allowed to change over time, capturing periods of high and low volatility.
-   **Student's t-Errors**: The model uses t-distributed errors for robustness against outliers and fat-tailed shocks in the data.
-   **Full-Sample Estimation**: A dedicated script to estimate the model over the entire dataset and visualize the posterior distributions of the unobserved components (trends, gaps, etc.).
-   **Recursive Forecasting Exercise**: A script to systematically evaluate the model's out-of-sample forecasting performance by recursively re-estimating the model and generating forecasts for multiple horizons (1, 4, 8, 12, and 16 quarters ahead).
-   **Forecast Evaluation**: Automatically calculates Root Mean Squared Forecast Errors (RMSFE) and log-predictive likelihoods to assess forecast accuracy.

## Project Structure

-   **`run_full_sample_estimation.m`**: The main script to run a single, in-depth estimation on the full data sample and plot the results.
-   **`run_forecasting_exercise.m`**: The main script to perform the recursive out-of-sample forecasting evaluation.
-   **`uc_sv_t_sampler.m`**: A core function containing the Gibbs sampler logic, called by both main scripts.
-   **`plot_empirical_facts.m`**: A script for initial exploratory data analysis and visualization.
-   **Helper Functions**:
    -   `tnormrnd.m`: Samples from a truncated normal distribution.
    -   `sample_nu_MH.m`: Samples the degrees-of-freedom parameter for the t-distribution.
-   **Data**:
    -   `CNdataset.csv`: The dataset used for the analysis.

## How to Use

### 1. Full-Sample Estimation
    -   Open `run_full_sample_estimation.m` in MATLAB.
    -   Configure the model specifications and MCMC settings at the top of the script.
    -   Ensure `CNdataset.csv` is in the same directory.
    -   Run the script. It will produce plots of the estimated trends, gaps, and time-varying parameters.

### 2. Forecasting Exercise
    -   Open `run_forecasting_exercise.m` in MATLAB.
    -   Choose the model variant you wish to test (e.g., `model = 1` for the full model).
    -   Run the script. This is computationally intensive and may take a significant amount of time.
    -   The script will output RMSFE tables and log-predictive likelihoods to the console and generate forecast plots.
