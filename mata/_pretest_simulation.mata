*! _pretest_simulation.mata
*! Monte Carlo Simulation Framework for Conditional Extrapolation Pre-Test
*!
*! Description:
*!   Implements Monte Carlo simulation for evaluating the finite-sample properties
*!   of conditionally valid confidence intervals under the conditional extrapolation
*!   assumption (Assumption 3). This module follows the numerical simulation design
*!   in Section 6 of Mikhaeil & Harshaw (2025).
*!
*! Simulation Design (Section 6):
*!   The Data Generating Process (DGP) constructs populations satisfying the
*!   conditional extrapolation assumption. For each Monte Carlo replication:
*!     1. Generate panel data with specified violation pattern
*!     2. Execute preliminary test (Theorem 1): phi = 1{S_hat_pre > M}
*!     3. If phi = 0 (pass), construct CI per Theorem 2 and check coverage
*!     4. Aggregate results to compute conditional coverage rates
*!
*! Key Performance Metrics:
*!   - Conditional coverage: P(tau_bar in CI | phi = 0) should be >= 1 - alpha
*!   - Probability of valid reporting: P(tau_bar in CI AND phi = 0)
*!   - Pre-test pass rate: P(phi = 0)
*!   - Mean CI width conditional on passing: E[width(CI) | phi = 0]
*!
*! Contents:
*!   _pretest_dgp_simulation()      - Data Generating Process for panel data
*!   _pretest_sim_single()          - Single Monte Carlo replication
*!   _pretest_coverage_simulation() - Full Monte Carlo simulation study
*!   _pretest_violation_dgp()       - Violation pattern generator
*!
*! Reference:
*!   Mikhaeil, J. M. and C. Harshaw. 2025. In Defense of the Pre-Test: Valid
*!   Inference when Testing Violations of Parallel Trends for Difference-in-
*!   Differences. arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470

version 17.0

mata:
mata set matastrict on

// ============================================================================
// DGP DATA GENERATION
// ============================================================================

/**
 * @function _pretest_dgp_simulation
 * @brief Generate simulated panel data following Section 6 DGP specification
 *
 * Implements the Data Generating Process from Section 6 (Numerical Simulations).
 * The DGP is constructed to satisfy the conditional extrapolation assumption
 * (Assumption 3) when S_pre <= M.
 *
 * Model Specification (Section 6):
 *
 *   Control group mean potential outcomes:
 *     E[Y^{(0)}_t | D=0] = rho * E[Y^{(0)}_{t-1} | D=0] + alpha*t + log(T)*(cos(t) + sin(t/2))
 *
 *   Treatment group mean potential outcomes under control:
 *     E[Y^{(0)}_t | D=1] = E[Y^{(0)}_t | D=0] + sum_{s=2}^{t} nu_s
 *
 *   Treatment group mean potential outcomes under treatment:
 *     E[Y^{(1)}_t | D=1] = E[Y^{(0)}_t | D=1] + tau * 1{t >= t0}
 *
 *   Parameters from Section 6: alpha = 0.3, rho = 0.7, tau = treatment effect
 *
 * Sampling (Section 6):
 *   Simulates independent draws with group-specific variances:
 *     Y_i | D=1, t >= t0 ~ N(E[Y^{(1)}_t | D=1], sigma_1^2)
 *     Y_i | D=1, t < t0  ~ N(E[Y^{(0)}_t | D=1], sigma_1^2)
 *     Y_i | D=0          ~ N(E[Y^{(0)}_t | D=0], sigma_0^2)
 *
 *   Default parameters: sigma_1 = 2.1 (treatment), sigma_0 = 1.5 (control)
 *
 * Post-Treatment Violations:
 *   This implementation sets post-treatment iterative violations to zero,
 *   carrying forward only the cumulative pre-treatment violation. This is a
 *   conservative design satisfying S_post = 0 <= S_pre, which is sufficient
 *   for validating coverage properties under Assumption 3.
 *
 * @param n           Number of units (total sample size)
 * @param T           Total number of time periods
 * @param t0          Treatment time (first post-treatment period, t0 >= 2)
 * @param tau         True average treatment effect on the treated (ATT)
 * @param p_treat     Treatment assignment probability P(D=1)
 * @param nu_vec      Iterative violations (nu_2, ..., nu_{t0-1})' of dim (T_pre-1) x 1
 * @param sigma_0     Standard deviation for control group (default: 1.5)
 * @param sigma_1     Standard deviation for treatment group (default: 2.1)
 * @param seed        Random seed for reproducibility
 *
 * @return Panel data matrix (n*T x 5):
 *         - Column 1: Y (observed outcome)
 *         - Column 2: D (treatment indicator, D in {0,1})
 *         - Column 3: time (period index, t in {1,...,T})
 *         - Column 4: unit (unit identifier)
 *         - Column 5: post (post-treatment indicator, 1{t >= t0})
 *
 * @note Block-adoption design only; staggered adoption not supported (Section 1.2)
 *
 * @see Mikhaeil & Harshaw (2025), Section 6
 */
real matrix _pretest_dgp_simulation(real scalar n,
                                     real scalar T,
                                     real scalar t0,
                                     real scalar tau,
                                     real scalar p_treat,
                                     | real colvector nu_vec,
                                     real scalar sigma_0,
                                     real scalar sigma_1,
                                     real scalar seed)
{
    real scalar i, t, idx, T_pre, n_treat
    real colvector D, nubar_vec
    real matrix data
    real scalar Y_it, post_t, nubar_t
    
    // ========================================
    // Model parameters (Section 6)
    // ========================================
    real scalar rho, alpha_coef
    rho = 0.7           // Autoregressive coefficient
    alpha_coef = 0.3    // Linear time trend coefficient
    
    // Default sampling parameters (Section 6)
    if (args() < 7 | missing(sigma_0)) {
        sigma_0 = 1.5   // Control group standard deviation
    }
    if (args() < 8 | missing(sigma_1)) {
        sigma_1 = 2.1   // Treatment group standard deviation
    }
    if (args() < 9 | missing(seed)) {
        seed = 12345
    }
    
    // Set random seed
    rseed(seed)
    
    // Pre-treatment periods
    T_pre = t0 - 1
    
    // Default: no violations (parallel trends hold)
    if (args() < 6 | rows(nu_vec) == 0) {
        nu_vec = J(T_pre - 1, 1, 0)
    }
    
    // Validate nu_vec length
    if (rows(nu_vec) != T_pre - 1) {
        _error("_pretest_dgp_simulation: nu_vec length must be T_pre-1 = " + 
               strofreal(T_pre - 1))
        return(J(0, 5, .))
    }
    
    // ========================================
    // Pre-compute control group means (Section 6)
    // ========================================
    // E[Y^{(0)}_t | D=0] = rho * E[Y^{(0)}_{t-1} | D=0] + alpha*t + log(T)*(cos(t) + sin(t/2))
    //
    // This is an AR(1) process with time-varying deterministic trend.
    // Recursive computation starting from E[Y^{(0)}_0 | D=0] = 0.
    real colvector mu_control
    mu_control = J(T, 1, 0)
    
    for (t = 1; t <= T; t++) {
        real scalar prev_mu, nonlinear_trend
        
        // Previous period mean
        if (t == 1) {
            prev_mu = 0
        }
        else {
            prev_mu = mu_control[t - 1]
        }
        
        // Nonlinear time trend component
        nonlinear_trend = ln(T) * (cos(t) + sin(t / 2))
        
        // AR(1) process with deterministic trend
        mu_control[t] = rho * prev_mu + alpha_coef * t + nonlinear_trend
    }
    
    // ========================================
    // Pre-compute cumulative violations (Section 2.1)
    // ========================================
    // nubar_t = sum_{s=2}^{t} nu_s represents the overall (cumulative) violation
    // at period t. The treatment group mean under control is:
    //   E[Y^{(0)}_t | D=1] = E[Y^{(0)}_t | D=0] + nubar_t
    nubar_vec = runningsum(nu_vec)
    
    // ========================================
    // Generate treatment assignment
    // ========================================
    D = J(n, 1, 0)
    n_treat = round(n * p_treat)
    
    // Assign first n_treat units to treatment
    for (i = 1; i <= n_treat; i++) {
        D[i] = 1
    }
    // Shuffle treatment assignment randomly
    D = D[jumble(range(1, n, 1)),]
    
    // ========================================
    // Initialize data matrix
    // ========================================
    // Columns: [Y, D, time, unit, post]
    data = J(n * T, 5, .)
    
    // ========================================
    // Generate panel data (Section 6)
    // ========================================
    idx = 0
    for (i = 1; i <= n; i++) {
        for (t = 1; t <= T; t++) {
            idx = idx + 1
            
            // Post-treatment indicator: 1{t >= t0}
            post_t = (t >= t0 ? 1 : 0)
            
            // ----------------------------------------
            // Step 1: Get control group mean for period t
            // ----------------------------------------
            real scalar mu_0t
            mu_0t = mu_control[t]
            
            // ----------------------------------------
            // Step 2: Compute cumulative violation for treated group
            // ----------------------------------------
            // nubar_t = sum_{s=2}^{t} nu_s (overall violation at period t)
            //
            // Index mapping:
            //   t = 1:                nubar_1 = 0 (no violation defined)
            //   t in {2,...,T_pre}:   nubar_t = nubar_vec[t-1]
            //   t > T_pre:            nubar_t = nubar_vec[T_pre-1] (carry forward)
            //
            // Note: Post-treatment periods carry forward the final pre-treatment
            // cumulative violation, implying zero post-treatment iterative violations.
            if (t == 1) {
                nubar_t = 0
            }
            else if (t <= T_pre) {
                nubar_t = nubar_vec[t - 1]
            }
            else {
                // Carry forward: S_post = 0 <= S_pre (conservative for coverage)
                nubar_t = nubar_vec[T_pre - 1]
            }
            
            // ----------------------------------------
            // Step 3: Generate outcome based on treatment status
            // ----------------------------------------
            real scalar mu_it, sigma_it
            
            if (D[i] == 0) {
                // Control group: Y_i ~ N(E[Y^{(0)}_t | D=0], sigma_0^2)
                mu_it = mu_0t
                sigma_it = sigma_0
            }
            else {
                // Treatment group
                if (t >= t0) {
                    // Post-treatment: Y_i ~ N(E[Y^{(1)}_t | D=1], sigma_1^2)
                    // where E[Y^{(1)}_t | D=1] = E[Y^{(0)}_t | D=0] + nubar_t + tau
                    mu_it = mu_0t + nubar_t + tau
                }
                else {
                    // Pre-treatment: Y_i ~ N(E[Y^{(0)}_t | D=1], sigma_1^2)
                    // where E[Y^{(0)}_t | D=1] = E[Y^{(0)}_t | D=0] + nubar_t
                    mu_it = mu_0t + nubar_t
                }
                sigma_it = sigma_1
            }
            
            // Draw from normal distribution
            Y_it = rnormal(1, 1, mu_it, sigma_it)
            
            // ----------------------------------------
            // Step 4: Store in data matrix
            // ----------------------------------------
            data[idx, 1] = Y_it       // Y (observed outcome)
            data[idx, 2] = D[i]       // D (treatment indicator)
            data[idx, 3] = t          // time (period index)
            data[idx, 4] = i          // unit (unit identifier)
            data[idx, 5] = post_t     // post (post-treatment indicator)
        }
    }
    
    return(data)
}


// ============================================================================
// SINGLE SIMULATION RUN
// ============================================================================

/**
 * @function _pretest_sim_single
 * @brief Execute single Monte Carlo replication with coverage evaluation
 *
 * Generates one dataset from the DGP, executes the preliminary test (Theorem 1),
 * constructs the conditionally valid confidence interval (Theorem 2), and
 * evaluates coverage. This function implements one iteration of the Monte Carlo
 * simulation study described in Section 6.
 *
 * Algorithm:
 *   1. Generate panel data from DGP with specified violations
 *   2. Estimate theta = (nu_2,...,nu_{t0-1}, delta_{t0},...,delta_T)'
 *   3. Compute S_hat_pre = severity(nu_hat, p) (Section 4.2)
 *   4. Execute preliminary test: phi = 1{S_hat_pre > M} (Theorem 1)
 *   5. If phi = 0 (pass):
 *      a. Compute kappa (Proposition 1) and f(alpha, Sigma) (Section 5.1)
 *      b. Construct CI per Theorem 2
 *      c. Check if tau_bar in CI
 *   6. Return coverage indicator and diagnostics
 *
 * @param tau       True ATT (tau_bar)
 * @param n         Number of units
 * @param T         Total number of time periods
 * @param t0        Treatment time (first post-treatment period)
 * @param p_treat   Treatment assignment probability P(D=1)
 * @param nu_vec    True iterative violations (T_pre-1 x 1)
 * @param M         Acceptable threshold for pre-test (Section 3.1)
 * @param p_norm    Lp-norm index for severity measure (p >= 1)
 * @param alpha     Significance level (e.g., 0.05 for 95% CI)
 * @param mode      Violation mode: "iterative" or "overall" (Appendix C)
 * @param S_sims    Monte Carlo simulations for critical value f(alpha, Sigma)
 * @param seed      Random seed for reproducibility
 *
 * @return Row vector (1 x 6):
 *         [1] covers:    1 if tau_bar in CI, 0 otherwise, missing if phi=1
 *         [2] phi:       Pre-test result (0 = pass, 1 = fail)
 *         [3] ci_lower:  CI lower bound (missing if phi=1)
 *         [4] ci_upper:  CI upper bound (missing if phi=1)
 *         [5] delta_bar: Average post-treatment DID estimate
 *         [6] S_pre:     Estimated pre-treatment severity
 *
 * @see Mikhaeil & Harshaw (2025), Theorems 1-2, Section 6
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
    
    // Generate data using default parameters from Section 6
    data = _pretest_dgp_simulation(n, T, t0, tau, p_treat, nu_vec, 1.5, 2.1, seed)
    
    Y = data[., 1]
    D = data[., 2]
    time_vec = data[., 3]
    
    T_pre = t0 - 1
    T_post = T - t0 + 1
    
    // Compute DID estimates
    delta_vec = _pretest_did_vector(Y, D, time_vec, t0, T)
    delta_bar = _pretest_did_avg(delta_vec)
    
    // Estimate iterative violations nu_hat = (nu_hat_2, ..., nu_hat_{t0-1})'
    // and overall violations nubar_hat (Appendix C)
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
    
    // Compute kappa constant (Proposition 1, Section 3.2)
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
    
    // Execute preliminary test (Theorem 1): phi = 1{S_hat_pre > M}
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
 * @brief Execute full Monte Carlo study for coverage validation
 *
 * Runs multiple independent simulation replications to evaluate the finite-sample
 * properties of the conditionally valid confidence intervals under Assumption 3:
 *
 *   1. Conditional coverage (Theorem 2): P(tau_bar in CI | phi = 0)
 *   2. Probability of valid reporting (Theorem 3): P(tau_bar in CI AND phi = 0)
 *   3. Pre-test pass rate: P(phi = 0)
 *
 * Theoretical Properties (Theorems 2-3):
 *   Under Assumption 3 with well-separated null (|S_pre - M| = omega(n^{-1/2})):
 *     - Conditional coverage >= 1 - alpha asymptotically
 *     - Probability of valid reporting >= 1 - alpha asymptotically
 *
 * @param n_sims    Number of Monte Carlo replications (5000 used in Section 6)
 * @param tau       True ATT (tau_bar)
 * @param n         Number of units per replication
 * @param T         Total number of time periods
 * @param t0        Treatment time (first post-treatment period)
 * @param p_treat   Treatment assignment probability P(D=1)
 * @param nu_vec    True iterative violations (T_pre-1 x 1)
 * @param M         Acceptable threshold for pre-test (Section 3.1)
 * @param p_norm    Lp-norm index for severity measure (p >= 1)
 * @param alpha     Significance level (e.g., 0.05)
 * @param mode      Violation mode: "iterative" or "overall" (Appendix C)
 * @param S_sims    MC simulations for f(alpha, Sigma) per replication
 * @param base_seed Base random seed (replication i uses seed = base_seed + i)
 *
 * @return Row vector (1 x 6):
 *         [1] coverage_rate:      P(tau in CI | phi=0) (conditional coverage)
 *         [2] pass_rate:          P(phi=0) (pre-test pass probability)
 *         [3] effective_coverage: P(tau in CI AND phi=0) (valid reporting)
 *         [4] n_passed:           Number of replications with phi=0
 *         [5] n_sims:             Total replications executed
 *         [6] mean_ci_width:      E[width(CI) | phi=0]
 *
 * Interpretation:
 *   - coverage_rate should be >= 1-alpha (e.g., >= 0.95 for alpha=0.05)
 *   - pass_rate depends on |S_pre - M| and sample size n
 *   - MC standard error: sqrt(coverage_rate*(1-coverage_rate)/n_passed)
 *
 * @see Mikhaeil & Harshaw (2025), Section 6, Figure 3
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
    // P(tau in CI AND phi=0) = P(tau in CI | phi=0) * P(phi=0)
    effective_coverage = pass_rate * coverage_rate
    
    return((coverage_rate, pass_rate, effective_coverage, n_passed, n_sims, mean_width))
}


// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

/**
 * @function _pretest_violation_dgp
 * @brief Generate structured iterative violation patterns for simulation
 *
 * Creates predefined violation patterns (nu_2, ..., nu_{t0-1}) for Monte Carlo
 * studies. These patterns facilitate investigation of CI properties under
 * various departures from parallel trends.
 *
 * Violation Semantics (Section 2.1):
 *   The returned vector represents iterative violations:
 *     nu_t = E[Y^{(0)}_t - Y^{(0)}_{t-1} | D=1] - E[Y^{(0)}_t - Y^{(0)}_{t-1} | D=0]
 *   Cumulative (overall) violations are: nubar_t = sum_{s=2}^{t} nu_s
 *
 * @param T_pre   Number of pre-treatment periods (T_pre = t0 - 1)
 * @param pattern Violation pattern type:
 *                - "zero":      nu_t = 0 for all t (exact parallel trends)
 *                - "constant":  nu_t = scale for all t (uniform deviation)
 *                - "linear":    nu_t = scale * t / (T_pre-1) (increasing)
 *                - "quadratic": nu_t = scale * (t / (T_pre-1))^2 (accelerating)
 * @param scale   Magnitude scaling factor
 *
 * @return Iterative violation vector (nu_2, ..., nu_{t0-1})' of dim (T_pre-1) x 1
 *
 * @note For nonlinear patterns as in Section 6, construct nu_vec directly
 *       using the scaling approach described therein.
 *
 * @see Mikhaeil & Harshaw (2025), Section 6
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
