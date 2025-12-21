*! _pretest_simulation.mata v0.1.0
*! Monte Carlo Simulation Framework for Coverage Validation
*!
*! Description:
*!   Implements Monte Carlo simulation framework for validating coverage
*!   rates of conditionally valid confidence intervals in the conditional
*!   extrapolation pre-test framework for difference-in-differences.
*!
*! Contents:
*!   - Data generating process (DGP) for panel data simulation
*!   - Single simulation run with coverage check
*!   - Full coverage simulation study
*!   - Violation pattern generators for simulation studies
*!
*! Reference:
*!   Mikhaeil, J.M. and C. Harshaw. 2025. In Defense of the Pre-Test: Valid
*!   Inference when Testing Violations of Parallel Trends for Difference-in-
*!   Differences. arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470
*!   Section 6 (Numerical Simulations)

version 17.0

mata:

// ============================================================================
// DGP DATA GENERATION
// ============================================================================

/**
 * @function _pretest_dgp_simulation
 * @brief Generate simulated panel data for coverage validation
 *
 * Implements the Data Generating Process from Section 6 (Numerical Simulations):
 *
 *   Y_it = alpha_i + beta_t + tau * D_i * post_t + nu_t * D_i * 1{t < t0} + eps_it
 *
 * Model components:
 *   - alpha_i ~ N(0, sigma_alpha^2): Unit fixed effects
 *   - beta_t = t: Linear time trend
 *   - tau: True average treatment effect on the treated (ATT)
 *   - nu_t: Pre-treatment parallel trend violations
 *   - eps_it ~ N(0, sigma_eps^2): Idiosyncratic errors
 *
 * @param n           Number of units
 * @param T           Number of time periods
 * @param t0          Treatment time (first post-treatment period)
 * @param tau         True ATT effect
 * @param p_treat     Treatment probability (proportion treated)
 * @param nu_vec      Pre-treatment violations (T_pre-1 x 1), empty for parallel trends
 * @param sigma_alpha Std dev of unit fixed effects (default: 1)
 * @param sigma_eps   Std dev of idiosyncratic errors (default: 1)
 * @param seed        Random seed for reproducibility
 *
 * @return Panel data matrix (n*T x 5):
 *         - Column 1: Y (outcome)
 *         - Column 2: D (treatment indicator)
 *         - Column 3: time (period)
 *         - Column 4: unit (unit ID)
 *         - Column 5: post (post-treatment indicator)
 *
 * @note Generates balanced panel with block-adoption design
 * @note Staggered adoption not supported
 *
 * @see Mikhaeil & Harshaw (2025), Section 6
 */
real matrix _pretest_dgp_simulation(real scalar n,
                                     real scalar T,
                                     real scalar t0,
                                     real scalar tau,
                                     real scalar p_treat,
                                     | real colvector nu_vec,
                                     real scalar sigma_alpha,
                                     real scalar sigma_eps,
                                     real scalar seed)
{
    real scalar i, t, idx, T_pre, n_treat
    real colvector alpha, D, eps_i
    real matrix data
    real scalar Y_it, post_t, nu_t
    
    // Default parameters
    if (args() < 7 | missing(sigma_alpha)) {
        sigma_alpha = 1
    }
    if (args() < 8 | missing(sigma_eps)) {
        sigma_eps = 1
    }
    if (args() < 9 | missing(seed)) {
        seed = 12345
    }
    
    // Set random seed
    rseed(seed)
    
    // Pre-treatment periods
    T_pre = t0 - 1
    
    // Default: no violations (parallel trends)
    if (args() < 6 | rows(nu_vec) == 0) {
        nu_vec = J(T_pre - 1, 1, 0)
    }
    
    // Validate nu_vec length
    if (rows(nu_vec) != T_pre - 1) {
        _error("_pretest_dgp_simulation: nu_vec length must be T_pre-1 = " + 
               strofreal(T_pre - 1))
        return(J(0, 5, .))
    }
    
    // Generate unit fixed effects: α_i ~ N(0, σ_α²)
    alpha = rnormal(n, 1, 0, sigma_alpha)
    
    // Generate treatment assignment
    D = J(n, 1, 0)
    n_treat = round(n * p_treat)
    
    // Randomly select treated units
    for (i = 1; i <= n_treat; i++) {
        D[i] = 1
    }
    // Shuffle treatment assignment
    D = D[jumble(range(1, n, 1)),]
    
    // Initialize data matrix: [Y, D, time, unit, post]
    data = J(n * T, 5, .)
    
    // Generate panel data
    idx = 0
    for (i = 1; i <= n; i++) {
        // Generate idiosyncratic errors for unit i
        eps_i = rnormal(T, 1, 0, sigma_eps)
        
        for (t = 1; t <= T; t++) {
            idx = idx + 1
            
            // Post-treatment indicator
            post_t = (t >= t0 ? 1 : 0)
            
            // Violation component for treated units in pre-treatment periods
            // ν_t applies to periods t = 2, ..., t0-1
            if (D[i] == 1 & t >= 2 & t <= T_pre) {
                nu_t = nu_vec[t - 1]  // nu_vec[1] = ν_2, etc.
            }
            else {
                nu_t = 0
            }
            
            // Generate outcome
            // Y_it = α_i + β_t + τ × D_i × post_t + ν_t × D_i + ε_it
            Y_it = alpha[i] + t + tau * D[i] * post_t + nu_t * D[i] + eps_i[t]
            
            // Store in data matrix
            data[idx, 1] = Y_it       // Y
            data[idx, 2] = D[i]       // D (treatment)
            data[idx, 3] = t          // time
            data[idx, 4] = i          // unit_id
            data[idx, 5] = post_t     // post
        }
    }
    
    return(data)
}


// ============================================================================
// SINGLE SIMULATION RUN
// ============================================================================

/**
 * @function _pretest_sim_single
 * @brief Execute single simulation replication and evaluate coverage
 *
 * Generates one dataset, runs the pre-test procedure, and checks whether
 * the resulting confidence interval covers the true ATT.
 *
 * @param tau      True ATT effect
 * @param n        Number of units
 * @param T        Number of time periods
 * @param t0       Treatment time (first post-treatment period)
 * @param p_treat  Treatment probability
 * @param nu_vec   Pre-treatment violations (T_pre-1 x 1)
 * @param M        Pre-test threshold
 * @param p_norm   Severity norm index (p >= 1)
 * @param alpha    Significance level
 * @param mode     "iterative" or "overall"
 * @param S_sims   Monte Carlo simulations for critical value
 * @param seed     Random seed for this replication
 *
 * @return Row vector (1 x 6):
 *         - [1] covers: 1 if CI contains true tau, 0 otherwise, . if pretest fails
 *         - [2] phi: pre-test result (0 = pass, 1 = fail)
 *         - [3] ci_lower: CI lower bound (. if pretest fails)
 *         - [4] ci_upper: CI upper bound (. if pretest fails)
 *         - [5] delta_bar: ATT estimate
 *         - [6] S_pre: severity estimate
 */
real rowvector _pretest_sim_single(real scalar tau,
                                   real scalar n,
                                   real scalar T,
                                   real scalar t0,
                                   real scalar p_treat,
                                   real colvector nu_vec,
                                   real scalar M,
                                   real scalar p_norm,
                                   real scalar alpha,
                                   string scalar mode,
                                   real scalar S_sims,
                                   real scalar seed)
{
    real matrix data, Sigma_hat
    real colvector Y, D, time_vec, delta_vec, nu_est, nubar_est
    real scalar delta_bar, S_pre, kappa, f_alpha, phi
    real scalar ci_lower, ci_upper, covers
    real scalar T_pre, T_post, overall_flag
    real rowvector ci_result
    
    // Generate data
    data = _pretest_dgp_simulation(n, T, t0, tau, p_treat, nu_vec, 1, 1, seed)
    
    Y = data[., 1]
    D = data[., 2]
    time_vec = data[., 3]
    
    T_pre = t0 - 1
    T_post = T - t0 + 1
    
    // Compute DID estimates
    delta_vec = _pretest_did_vector(Y, D, time_vec, t0, T)
    delta_bar = _pretest_did_avg(delta_vec)
    
    // Compute violations
    // Use _pretest_nu_vector for iterative violations
    // Use _pretest_nu_to_nubar for overall (cumulative) violations
    nu_est = _pretest_nu_vector(Y, D, time_vec, t0)
    nubar_est = _pretest_nu_to_nubar(nu_est)
    
    // Compute covariance matrix
    Sigma_hat = _pretest_covariance(Y, D, time_vec, t0, T, 1)
    
    // Compute severity
    overall_flag = (mode == "overall" ? 1 : 0)
    
    if (overall_flag == 1) {
        S_pre = _pretest_severity(nubar_est, p_norm)
    }
    else {
        S_pre = _pretest_severity(nu_est, p_norm)
    }
    
    // Compute kappa: _pretest_kappa(T_post, p, overall_mode)
    kappa = _pretest_kappa(T_post, p_norm, overall_flag)
    
    // Compute critical value
    if (overall_flag == 1) {
        real matrix Sigma_overall
        Sigma_overall = _pretest_transform_sigma(Sigma_hat, T_pre, T_post)
        f_alpha = _pretest_critical_value(Sigma_overall, alpha, S_sims,
                                          T_pre, T_post, p_norm, kappa, 1, seed + 1000)
    }
    else {
        f_alpha = _pretest_critical_value(Sigma_hat, alpha, S_sims,
                                          T_pre, T_post, p_norm, kappa, 0, seed + 1000)
    }
    
    // Pretest check
    phi = _pretest_test(S_pre, M)
    
    // Compute CI if pretest passes
    if (phi == 0) {
        if (overall_flag == 1) {
            ci_result = _pretest_ci_overall(delta_bar, S_pre, f_alpha, rows(Y))
        }
        else {
            ci_result = _pretest_ci_iterative(delta_bar, S_pre, kappa, f_alpha, rows(Y))
        }
        ci_lower = ci_result[1]
        ci_upper = ci_result[2]
        
        // Check if CI covers true tau
        covers = (ci_lower <= tau & tau <= ci_upper ? 1 : 0)
    }
    else {
        ci_lower = .
        ci_upper = .
        covers = .
    }
    
    return((covers, phi, ci_lower, ci_upper, delta_bar, S_pre))
}


// ============================================================================
// COVERAGE SIMULATION MAIN FUNCTION
// ============================================================================

/**
 * @function _pretest_coverage_simulation
 * @brief Execute full Monte Carlo coverage simulation study
 *
 * Runs multiple simulation replications to evaluate the conditional
 * coverage properties of the pre-test confidence intervals.
 *
 * @param n_sims    Number of simulation replications
 * @param tau       True ATT effect
 * @param n         Number of units per simulation
 * @param T         Number of time periods
 * @param t0        Treatment time (first post-treatment period)
 * @param p_treat   Treatment probability
 * @param nu_vec    Pre-treatment violations (T_pre-1 x 1)
 * @param M         Pre-test threshold
 * @param p_norm    Severity norm index (p >= 1)
 * @param alpha     Significance level
 * @param mode      "iterative" or "overall"
 * @param S_sims    Monte Carlo simulations for critical value (per replication)
 * @param base_seed Base random seed
 *
 * @return Row vector (1 x 6):
 *         - [1] coverage_rate: Conditional coverage P(tau in CI | pretest passed)
 *         - [2] pass_rate: Proportion of replications where pretest passed
 *         - [3] effective_coverage: pass_rate * coverage_rate
 *         - [4] n_passed: Number of replications where pretest passed
 *         - [5] n_sims: Total replications run
 *         - [6] mean_ci_width: Average CI width among passed replications
 *
 * @note Target conditional coverage: >= 1 - alpha (accounting for MC error)
 *
 * @see Mikhaeil & Harshaw (2025), Section 6
 */
real rowvector _pretest_coverage_simulation(real scalar n_sims,
                                            real scalar tau,
                                            real scalar n,
                                            real scalar T,
                                            real scalar t0,
                                            real scalar p_treat,
                                            real colvector nu_vec,
                                            real scalar M,
                                            real scalar p_norm,
                                            real scalar alpha,
                                            string scalar mode,
                                            real scalar S_sims,
                                            real scalar base_seed)
{
    real scalar i, n_passed, n_covered, total_width
    real rowvector sim_result
    real scalar coverage_rate, pass_rate, effective_coverage, mean_width
    
    // Initialize counters
    n_passed = 0
    n_covered = 0
    total_width = 0
    
    // Run simulations
    for (i = 1; i <= n_sims; i++) {
        // Run single simulation with unique seed
        sim_result = _pretest_sim_single(tau, n, T, t0, p_treat, nu_vec,
                                         M, p_norm, alpha, mode, S_sims,
                                         base_seed + i)
        
        // Check if pretest passed (phi = 0)
        if (sim_result[2] == 0) {
            n_passed = n_passed + 1
            
            // Check if CI covers true tau
            if (sim_result[1] == 1) {
                n_covered = n_covered + 1
            }
            
            // Accumulate CI width
            total_width = total_width + (sim_result[4] - sim_result[3])
        }
        
        // Progress indicator (every 100 sims)
        if (mod(i, 100) == 0) {
            printf("  Simulation %g/%g completed\n", i, n_sims)
        }
    }
    
    // Calculate rates
    if (n_passed > 0) {
        coverage_rate = n_covered / n_passed
        mean_width = total_width / n_passed
    }
    else {
        coverage_rate = .
        mean_width = .
    }
    
    pass_rate = n_passed / n_sims
    effective_coverage = pass_rate * coverage_rate
    
    return((coverage_rate, pass_rate, effective_coverage, n_passed, n_sims, mean_width))
}


// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

/**
 * @function _pretest_violation_dgp
 * @brief Generate structured violation patterns for DGP
 *
 * Creates predefined violation patterns for simulation studies.
 *
 * @param T_pre   Number of pre-treatment periods
 * @param pattern Violation pattern:
 *                - "zero": no violations (parallel trends)
 *                - "constant": constant violation at each period
 *                - "linear": linearly increasing violations
 *                - "quadratic": quadratically increasing violations
 * @param scale   Scale factor for violation magnitudes
 *
 * @return Violation vector of dimension (T_pre-1) x 1
 */
real colvector _pretest_violation_dgp(real scalar T_pre,
                                       string scalar pattern,
                                       real scalar scale)
{
    real colvector nu_vec
    real scalar n_nu, t
    
    n_nu = T_pre - 1
    
    if (n_nu < 1) {
        return(J(0, 1, .))
    }
    
    if (pattern == "zero") {
        nu_vec = J(n_nu, 1, 0)
    }
    else if (pattern == "constant") {
        nu_vec = J(n_nu, 1, scale)
    }
    else if (pattern == "linear") {
        nu_vec = J(n_nu, 1, .)
        for (t = 1; t <= n_nu; t++) {
            nu_vec[t] = scale * t / n_nu
        }
    }
    else if (pattern == "quadratic") {
        nu_vec = J(n_nu, 1, .)
        for (t = 1; t <= n_nu; t++) {
            nu_vec[t] = scale * (t / n_nu)^2
        }
    }
    else {
        _error("Unknown pattern: " + pattern)
        return(J(n_nu, 1, 0))
    }
    
    return(nu_vec)
}


end
