*! _pretest_main.mata
*! Main Computation Engine for Conditional Extrapolation Pre-Test
*!
*! Description:
*!   Core computational module implementing the conditional extrapolation
*!   pre-test framework for difference-in-differences (DID) research designs.
*!   This framework provides valid inference conditional on passing a
*!   preliminary test for violations of parallel trends in pre-treatment periods.
*!
*! Theoretical Framework:
*!   The conditional extrapolation assumption (Assumption 3) posits that if
*!   pre-treatment violations are below an acceptable threshold M, then
*!   post-treatment violations are bounded by pre-treatment violations:
*!
*!       S_pre <= M  ==>  S_post <= S_pre
*!
*!   Under this assumption, confidence intervals are conditionally valid:
*!       P(tau_bar in CI | phi = 0) >= 1 - alpha  (asymptotically)
*!
*! Algorithm Pipeline:
*!   1. Estimate DID parameters: delta_t for t in {t0, ..., T}
*!   2. Estimate violations: nu_t for t in {2, ..., t0-1}
*!   3. Compute asymptotic covariance matrix Sigma (Assumption 1)
*!   4. Calculate severity measure S_pre = ||nu||_p (Section 3.1)
*!   5. Compute bias bound constant kappa (Proposition 1)
*!   6. Monte Carlo critical value f(alpha, Sigma) (Appendix F)
*!   7. Execute pre-test: phi = 1{S_pre > M} (Theorem 1)
*!   8. Construct conditionally valid CI (Theorem 2)
*!
*! Module Dependencies:
*!   _pretest_estimators.mata  - DID estimators (Section 2.1)
*!   _pretest_violations.mata  - Violation estimators (Section 2.1)
*!   _pretest_covariance.mata  - Covariance estimation (Section 2.2-2.3)
*!   _pretest_psi.mata         - Critical value computation (Appendix F)
*!   _pretest_ci.mata          - Confidence intervals (Section 5.1)
*!   _pretest_kappa.mata       - Kappa constant (Section 3.2)
*!
*! Reference:
*!   Mikhaeil, J.M. and C. Harshaw. 2025. In Defense of the Pre-Test: Valid
*!   Inference when Testing Violations of Parallel Trends for Difference-in-
*!   Differences. arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470

version 17.0

mata:
mata set matastrict on

// ============================================================================
// MAIN COMPUTATION FUNCTION
// ============================================================================

/**
 * @function _pretest_main_compute
 * @brief Primary computation interface for the pretest command
 *
 * Orchestrates the complete conditional extrapolation pre-test pipeline.
 * This function implements the methodology from Mikhaeil & Harshaw (2025)
 * for valid inference in difference-in-differences designs.
 *
 * Computation Steps:
 *   Step 1:  Extract parameters from Stata environment
 *   Step 2:  Load data and map time indices to consecutive integers
 *   Step 3:  Compute DID estimates delta_t for t in {t0, ..., T}
 *   Step 4:  Compute violation estimates nu_t for t in {2, ..., t0-1}
 *   Step 5:  Estimate asymptotic covariance matrix Sigma
 *   Step 6:  Calculate severity S_pre = ((1/(T_pre-1)) sum |nu_t|^p)^{1/p}
 *   Step 7:  Compute bias bound constant kappa (Proposition 1)
 *   Step 8:  Monte Carlo critical value f(alpha, Sigma) via simulation
 *   Step 9:  Execute pre-test phi = 1{S_pre > M} (Theorem 1)
 *   Step 10: Construct conditionally valid CI (Theorem 2)
 *   Step 11: Assemble parameter vector theta
 *   Step 12: Compute severity standard error via Delta Method
 *   Step 13: Compute conventional CI for comparison
 *
 * @param outcome_var   Outcome variable name (Y) in Stata dataset
 * @param treatment_var Treatment indicator variable (D in {0, 1})
 * @param time_var      Time period variable (t = 1, ..., T)
 * @param cluster_var   Cluster variable for robust SEs (empty string if none)
 *
 * Input Scalars (prefixed with __pretest_):
 *   T, t0, T_pre, T_post - Time structure: total periods, treatment time
 *   p                    - Norm index for severity (p >= 1, inf for L-infinity)
 *   alpha                - Significance level (typically 0.05)
 *   threshold            - Acceptable threshold M for pre-test
 *   sims                 - Number of Monte Carlo simulations (default: 5000)
 *   seed                 - Random seed for reproducibility
 *   is_panel             - Data type indicator (1=panel, 0=repeated cross-section)
 *   mode_overall         - Violation mode (0=iterative, 1=overall)
 *
 * Output Scalars (prefixed with __pretest_):
 *   delta_bar            - Average post-treatment DID estimate
 *   S_pre                - Pre-treatment severity measure
 *   kappa                - Bias bound constant from Proposition 1
 *   f_alpha              - Monte Carlo critical value f(alpha, Sigma)
 *   phi                  - Pre-test result: 0=PASS, 1=FAIL, .=invalid
 *   ci_lower, ci_upper   - Conditionally valid CI bounds (missing if phi != 0)
 *   se_delta_bar         - Standard error of delta_bar
 *   ci_conv_lower/upper  - Conventional CI bounds (for comparison)
 *
 * @note Results stored in Stata scalars/matrices for retrieval by ado file
 * @note Conventional CI assumes exact parallel trends (for comparison only)
 *
 * @see Mikhaeil & Harshaw (2025), Sections 2-5, Theorems 1-2, Proposition 1
 */
void _pretest_main_compute(string scalar outcome_var,
                           string scalar treatment_var,
                           string scalar time_var,
                           string scalar cluster_var)
{
    // ========================================
    // Step 1: Extract parameters from Stata
    // ========================================
    real scalar T_mata, t0_mata, T_pre_m, T_post_m
    real scalar p_m, alpha_m, threshold_m
    real scalar S_sims, seed_val, is_panel_m, overall_flag
    
    T_mata = st_numscalar("__pretest_T")
    t0_mata = st_numscalar("__pretest_t0")
    T_pre_m = st_numscalar("__pretest_T_pre")
    T_post_m = st_numscalar("__pretest_T_post")
    p_m = st_numscalar("__pretest_p")
    alpha_m = st_numscalar("__pretest_alpha")
    threshold_m = st_numscalar("__pretest_threshold")
    S_sims = st_numscalar("__pretest_sims")
    seed_val = st_numscalar("__pretest_seed")
    is_panel_m = st_numscalar("__pretest_is_panel")
    overall_flag = st_numscalar("__pretest_mode_overall")
    
    // ========================================
    // Step 2: Load data and map time indices
    // ========================================
    real colvector Y, D, time_orig, time_mapped
    real matrix time_vals_m
    real scalar i, j
    
    Y = st_data(., outcome_var)
    D = st_data(., treatment_var)
    time_orig = st_data(., time_var)
    
    // Get time mapping matrix
    time_vals_m = st_matrix("__pretest_time_vals")
    
    // Map original time values to consecutive indices
    time_mapped = J(rows(time_orig), 1, .)
    for (i = 1; i <= rows(time_orig); i++) {
        for (j = 1; j <= rows(time_vals_m); j++) {
            if (time_orig[i] == time_vals_m[j, 1]) {
                time_mapped[i] = j
                break
            }
        }
    }
    
    // ========================================
    // Step 3: DID Estimation (Section 2.1)
    // ========================================
    // The identified DID estimand for period t >= t0 (Equation 171):
    //
    //   delta_t = E[Y_t - Y_{t0} | D=1] - E[Y_t - Y_{t0} | D=0]
    //
    // Sample analog using group-time means:
    //
    //   delta_hat_t = (Y_bar_{t,1} - Y_bar_{t0,1}) - (Y_bar_{t,0} - Y_bar_{t0,0})
    //
    // The average post-treatment DID (Section 2.1):
    //
    //   delta_bar = (1/T_post) * sum_{t=t0}^T delta_t
    //
    // Under parallel trends: delta_bar = tau_bar (the average ATT)
    // Under Assumption 3:    |tau_bar - delta_bar| <= kappa * S_pre
    
    real colvector delta_vec
    real scalar delta_bar
    
    delta_vec = _pretest_did_vector(Y, D, time_mapped, t0_mata, T_mata)
    delta_bar = _pretest_did_avg(delta_vec)
    
    // Store results
    st_matrix("__pretest_delta", delta_vec)
    st_numscalar("__pretest_delta_bar", delta_bar)
    
    // ========================================
    // Step 4: Violation Estimation (Section 2.1)
    // ========================================
    // Iterative violations measure period-to-period deviations from parallel trends:
    //
    //   nu_t = E[Y^{(0)}_t - Y^{(0)}_{t-1} | D=1] - E[Y^{(0)}_t - Y^{(0)}_{t-1} | D=0]
    //
    // Sample analog (Equation 184):
    //
    //   nu_hat_t = (Y_bar_{t,1} - Y_bar_{t-1,1}) - (Y_bar_{t,0} - Y_bar_{t-1,0})
    //
    // Overall violations measure cumulative deviation (Appendix C):
    //
    //   nubar_t = sum_{s=2}^t nu_s
    //
    // Equivalence relation (Equation 191): nubar_t = sum nu_s, nu_t = nubar_t - nubar_{t-1}
    real colvector nu_vec, nubar_vec
    
    nu_vec = _pretest_nu_vector(Y, D, time_mapped, t0_mata)
    nubar_vec = _pretest_nu_to_nubar(nu_vec)
    
    st_matrix("__pretest_nu", nu_vec)
    
    // ========================================
    // Step 5: Asymptotic Covariance Estimation (Section 2.2-2.3)
    // ========================================
    // Assumption 1 (Asymptotic Normality):
    //   sqrt(n) * (theta_hat - theta) -->d N(0, Sigma)
    //
    // Assumption 2 (Consistent Variance Estimation):
    //   Sigma_hat -->p Sigma
    //
    // Parameter vector structure (Equation 205):
    //   theta = (nu_2, ..., nu_{t0-1}, delta_{t0}, ..., delta_T)'
    //
    // Covariance matrix structure:
    //   Sigma = | Sigma_nu      Sigma_nu,delta |
    //           | Sigma_delta,nu  Sigma_delta  |
    real matrix Sigma_hat, Sigma_nu, Sigma_hat_overall
    real colvector cluster_vec
    real scalar sigma_valid
    
    if (cluster_var != "") {
        cluster_vec = st_data(., cluster_var)
        Sigma_hat = _pretest_covariance(Y, D, time_mapped, t0_mata, T_mata, is_panel_m, cluster_vec)
    }
    else {
        Sigma_hat = _pretest_covariance(Y, D, time_mapped, t0_mata, T_mata, is_panel_m)
    }
    
    // Check if covariance matrix is valid (no missing values)
    sigma_valid = !hasmissing(Sigma_hat)
    
    st_matrix("__pretest_Sigma", Sigma_hat)
    
    // Extract Sigma_nu for violation component
    Sigma_nu = _pretest_extract_sigma_nu(Sigma_hat, T_pre_m)
    
    // Transform covariance matrix for overall mode if needed
    if (overall_flag == 1) {
        Sigma_hat_overall = _pretest_transform_sigma(Sigma_hat, T_pre_m, T_post_m)
    }
    
    // ========================================
    // Step 6: Severity Measurement (Section 3.1)
    // ========================================
    // The severity of pre-treatment parallel trends violations is defined as:
    //
    //   S_pre = ((1/(T_pre-1)) * sum_{t=2}^{t0-1} |nu_t|^p)^{1/p}
    //
    // This quantity aggregates violations using the Lp norm framework:
    //   p = 1:   Mean absolute violation (simple average)
    //   p = 2:   Root mean squared violation (default, sensitive to large violations)
    //   p = inf: Maximum absolute violation (worst-case single period)
    //
    // The choice of p reflects how violations are aggregated across pre-treatment periods.
    // See Section 3.1 for guidance on choosing p based on substantive considerations.
    real scalar S_pre
    
    if (overall_flag == 0) {
        S_pre = _pretest_severity(nu_vec, p_m)
    }
    else {
        S_pre = _pretest_severity(nubar_vec, p_m)
    }
    
    st_numscalar("__pretest_S_pre", S_pre)
    
    // ========================================
    // Step 7: Bias Bound Constant Kappa
    // ========================================
    // Reference: Section 3.2, Proposition 1
    //
    // Under conditional extrapolation (Assumption 3), if S_pre <= M:
    //   |tau_bar - delta_bar| <= kappa * S_pre
    //
    // Kappa formula (via Holder's inequality):
    //   kappa = ((1/T_post) * sum_{t=1}^{T_post} t^q)^{1/q}
    //   where q is the Holder conjugate: 1/p + 1/q = 1
    //
    // Note: In overall mode, kappa = 1 (violations directly bound bias)
    real scalar kappa
    kappa = _pretest_kappa(T_post_m, p_m, overall_flag)
    
    st_numscalar("__pretest_kappa", kappa)
    
    // ========================================
    // Step 8: Critical Value via Monte Carlo
    // ========================================
    // Reference: Section 5.1, Theorem 2; Appendix D
    //
    // The critical value f(alpha, Sigma) is defined as:
    //   f(alpha, Sigma) = inf{c : P(psi(Z) >= c) <= alpha}
    //   where Z ~ N(0, Sigma)
    //
    // Monte Carlo algorithm:
    //   1. Draw Z^(s) ~ N(0, Sigma) for s = 1,...,S
    //   2. Compute psi^(s) = |mean(Z_post)| + kappa * severity(Z_pre)
    //   3. Return (1-alpha) quantile of {psi^(s)}
    real scalar f_alpha
    
    // Only compute critical value if covariance matrix is valid
    if (sigma_valid) {
        if (overall_flag == 1) {
            f_alpha = _pretest_critical_value(Sigma_hat_overall, alpha_m, S_sims, T_pre_m, T_post_m, p_m, kappa, 1, seed_val)
        }
        else {
            f_alpha = _pretest_critical_value(Sigma_hat, alpha_m, S_sims, T_pre_m, T_post_m, p_m, kappa, 0, seed_val)
        }
    }
    else {
        // Covariance matrix invalid (some periods have no observations), set f_alpha to missing
        f_alpha = .
    }
    
    st_numscalar("__pretest_f_alpha", f_alpha)
    
    // ========================================
    // Step 9: Pre-Test Evaluation
    // ========================================
    // Reference: Section 4.2, Theorem 1 (Asymptotic Consistency)
    //
    // Pre-test function:
    //   phi = 1{S_pre > M}
    //
    // Interpretation:
    //   phi = 0 (PASS): S_pre <= M, extrapolation condition satisfied,
    //                   proceed with conditionally valid inference
    //   phi = 1 (FAIL): S_pre > M, extrapolation not justified,
    //                   DID analysis may be invalid
    real scalar M, phi, pretest_passed, n
    
    n = rows(Y)
    
    // Default threshold: auto-compute (rarely used)
    if (threshold_m == -999) {
        M = S_pre + 0.01
    }
    else {
        M = threshold_m
    }
    
    phi = _pretest_test(S_pre, M)
    
    // Only pass pretest if all computations are valid
    // Conditions: sigma_valid, f_alpha not missing, S_pre not missing, phi=0
    if (phi == 0 & sigma_valid & !missing(f_alpha) & !missing(S_pre)) {
        pretest_passed = 1
    }
    else {
        pretest_passed = 0
        // If failure due to data issues, set phi to missing
        if (!sigma_valid | missing(f_alpha) | missing(S_pre)) {
            phi = .  // Missing: indicates data issues preventing valid inference
        }
    }
    
    st_numscalar("__pretest_phi", phi)
    st_numscalar("__pretest_n", n)
    
    // ========================================
    // Step 10: Conditionally Valid Confidence Interval
    // ========================================
    // Reference: Section 5.1, Theorem 2
    //
    // The conditionally valid CI for the average post-treatment ATT:
    //
    // Iterative mode:
    //   CI = delta_bar +/- {kappa * S_pre + f(alpha, Sigma) / sqrt(n)}
    //
    // Overall mode:
    //   CI = delta_bar +/- {S_pre + f(alpha, Sigma) / sqrt(n)}
    //
    // Property: P(tau_bar in CI | phi = 0) >= 1 - alpha asymptotically
    real scalar ci_lower, ci_upper
    real rowvector ci_result
    
    if (pretest_passed == 1) {
        if (overall_flag == 1) {
            ci_result = _pretest_ci_overall(delta_bar, S_pre, f_alpha, n)
        }
        else {
            ci_result = _pretest_ci_iterative(delta_bar, S_pre, kappa, f_alpha, n)
        }
        ci_lower = ci_result[1]
        ci_upper = ci_result[2]
    }
    else {
        ci_lower = .
        ci_upper = .
    }
    
    st_numscalar("__pretest_ci_lower", ci_lower)
    st_numscalar("__pretest_ci_upper", ci_upper)
    
    // ========================================
    // Step 11: Parameter Vector Theta
    // ========================================
    // Reference: Section 2.2, Equation for theta
    //
    // The full parameter vector for asymptotic analysis:
    //   theta = (nu_2, ..., nu_{t0-1}, delta_{t0}, ..., delta_T)'
    //
    // Dimension: dim(theta) = (T_pre - 1) + T_post = T - 1
    real colvector theta_vec
    
    theta_vec = nu_vec \ delta_vec
    st_matrix("__pretest_theta", theta_vec)
    
    // ========================================
    // Step 12: Severity Standard Error (Delta Method)
    // ========================================
    // Reference: Section 5.1, variance estimation
    //
    // The Delta Method approximation for SE(S_pre):
    //   Var(S_pre) approx grad_g' * Sigma_nu * grad_g
    //   SE(S_pre) = sqrt(Var(S_pre))
    //
    // where grad_g = partial S_pre / partial nu is the gradient vector
    //
    // Note: Not applicable for p = infinity (L-inf norm is non-differentiable)
    real scalar S_pre_se
    
    // Delta Method not applicable when p = infinity, return missing
    if (p_m >= 1e10) {
        S_pre_se = .
    }
    else {
        S_pre_se = _pretest_severity_se(nu_vec, p_m, Sigma_nu)
    }
    
    st_numscalar("__pretest_S_pre_se", S_pre_se)
    
    // ========================================
    // Step 13: Conventional DID Confidence Interval
    // ========================================
    // Standard CI assuming exact parallel trends (for comparison):
    //
    //   CI_conv = delta_bar +/- z_{1-alpha/2} * SE(delta_bar)
    //
    // where:
    //   SE(delta_bar) = sqrt(Var(delta_bar))
    //   Var(delta_bar) = (1/T_post^2) * 1' * Sigma_delta * 1 / n
    //
    // Note: This CI is only valid under exact parallel trends;
    //       use the conditional CI from Step 10 under Assumption 3
    real scalar se_delta_bar, conv_ci_lower, conv_ci_upper
    real rowvector conv_ci_result
    real matrix Sigma_delta
    real scalar var_delta_bar, i_d, j_d
    
    // Extract delta portion from Sigma_hat (lower-right T_post x T_post submatrix)
    // Sigma_hat structure: | Sigma_nu      Sigma_nu,delta |
    //                      | Sigma_delta,nu  Sigma_delta  |
    // nu portion occupies T_pre-1 rows/cols, delta portion occupies T_post rows/cols
    real scalar n_nu_dim, delta_start_idx, delta_end_idx
    
    n_nu_dim = T_pre_m - 1  // Dimension of iterative violations
    delta_start_idx = n_nu_dim + 1
    delta_end_idx = rows(Sigma_hat)
    
    // If covariance matrix is invalid, set se_delta_bar to missing
    if (!sigma_valid) {
        se_delta_bar = .
    }
    else if (delta_end_idx >= delta_start_idx & T_post_m > 0) {
        Sigma_delta = Sigma_hat[delta_start_idx..delta_end_idx, delta_start_idx..delta_end_idx]
        
        // Var(delta_bar) = (1/T_post^2) * 1'*Sigma_delta*1 / n
        // i.e., sum of all elements in Sigma_delta divided by (T_post^2 * n)
        var_delta_bar = 0
        real scalar valid_count
        valid_count = 0
        for (i_d = 1; i_d <= rows(Sigma_delta); i_d++) {
            for (j_d = 1; j_d <= cols(Sigma_delta); j_d++) {
                if (!missing(Sigma_delta[i_d, j_d])) {
                    var_delta_bar = var_delta_bar + Sigma_delta[i_d, j_d]
                    valid_count++
                }
            }
        }
        
        // If no valid elements, set to missing
        if (valid_count == 0) {
            se_delta_bar = .
        }
        else {
            var_delta_bar = var_delta_bar / (T_post_m * T_post_m * n)
            
            // Handle numerical issues
            if (var_delta_bar < 0) {
                var_delta_bar = 0
            }
            
            se_delta_bar = sqrt(var_delta_bar)
        }
    }
    else {
        se_delta_bar = .
    }
    
    // Compute conventional CI
    conv_ci_result = _pretest_ci_conventional(delta_bar, se_delta_bar, alpha_m)
    conv_ci_lower = conv_ci_result[1]
    conv_ci_upper = conv_ci_result[2]
    
    st_numscalar("__pretest_se_delta_bar", se_delta_bar)
    st_numscalar("__pretest_ci_conv_lower", conv_ci_lower)
    st_numscalar("__pretest_ci_conv_upper", conv_ci_upper)
}

end
