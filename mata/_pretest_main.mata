*! _pretest_main.mata
*! Main Computation Engine for Conditional Extrapolation Pre-Test
*!
*! Purpose:
*!   Core module implementing the conditional extrapolation pre-test framework
*!   for difference-in-differences designs. Provides valid inference conditional
*!   on passing the preliminary test for pre-treatment parallel trends violations.
*!
*! Theoretical Foundation:
*!   Conditional Extrapolation Assumption (Assumption 3, Section 3.1):
*!     If S_pre <= M, then S_post <= S_pre
*!   
*!   This implies conditional validity (Theorem 2, Section 5.1):
*!     lim inf P(tau_bar in CI | phi = 0) >= 1 - alpha
*!
*! Algorithm Pipeline:
*!   Steps 1-3: Data preparation and DID estimation (Section 2.1)
*!   Steps 4-5: Violation estimation and covariance (Section 2.2)
*!   Step 6:    Severity measure S_pre (Section 3.1)
*!   Step 7:    Bias bound constant kappa (Proposition 1, Section 3.2)
*!   Step 8:    Critical value f(alpha, Sigma) (Appendix F)
*!   Step 9:    Preliminary test phi = 1{S_hat_pre > M} (Theorem 1, Section 4.2)
*!   Step 10:   Conditionally valid confidence interval (Theorem 2, Section 5.1)
*!
*! Dependencies:
*!   _pretest_estimators, _pretest_violations, _pretest_covariance,
*!   _pretest_psi, _pretest_ci, _pretest_kappa
*!
*! Reference:
*!   Mikhaeil, J. M. and C. Harshaw (2025). In Defense of the Pre-Test:
*!   Valid Inference when Testing Violations of Parallel Trends for
*!   Difference-in-Differences. arXiv:2510.26470.

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
    // Step 2.5: Filter observations with missing values
    // ========================================
    // Exclude observations where Y, D, or time_mapped are missing to ensure
    // consistent sample size reporting and valid group-time cell estimation.
    // Note: Must filter AFTER time mapping because mapping can introduce
    // additional missing values if time_orig values are not in time_vals_m.
    real colvector valid_rows
    real scalar n_total, n_valid
    
    n_total = rows(Y)
    valid_rows = !rowmissing(Y) :& !rowmissing(D) :& !rowmissing(time_mapped)
    n_valid = sum(valid_rows)
    
    // Filter to valid observations only
    if (n_valid < n_total) {
        Y = select(Y, valid_rows)
        D = select(D, valid_rows)
        time_mapped = select(time_mapped, valid_rows)
    }
    
    // ========================================
    // Step 3: DID Estimation (Section 2.1)
    // ========================================
    // The identified difference-in-differences estimand for period t >= t0:
    //
    //   delta_t = E[Y_t - Y_{t0} | D=1] - E[Y_t - Y_{t0} | D=0]
    //
    // Sample analog using group-time cell means:
    //
    //   delta_hat_t = (Y_bar_{t,1} - Y_bar_{t0,1}) - (Y_bar_{t,0} - Y_bar_{t0,0})
    //
    // The average post-treatment DID estimand (Section 2.1):
    //
    //   delta_bar = (1/T_post) * sum_{t=t0}^{T} delta_t
    //
    // Under exact parallel trends: delta_bar = tau_bar (the average ATT)
    // Under Assumption 3 (Conditional Extrapolation):
    //   |tau_bar - delta_bar| <= kappa * S_pre  (Proposition 1)
    
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
    // Sample analog (Section 2.2):
    //
    //   nu_hat_t = (Y_bar_{t,1} - Y_bar_{t-1,1}) - (Y_bar_{t,0} - Y_bar_{t-1,0})
    //
    // Overall violations measure cumulative deviation from reference period t0:
    //
    //   nubar_t = sum_{s=2}^{t} nu_s
    //
    // Equivalence relation (Section 2.1): nubar_t = sum nu_s, nu_t = nubar_t - nubar_{t-1}
    real colvector nu_vec, nubar_vec
    
    nu_vec = _pretest_nu_vector(Y, D, time_mapped, t0_mata)
    nubar_vec = _pretest_nu_to_nubar(nu_vec)
    
    st_matrix("__pretest_nu", nu_vec)
    
    // ========================================
    // Step 5: Asymptotic Covariance Estimation (Sections 2.2-2.3)
    // ========================================
    // Assumption 1 (Asymptotic Normality, Section 2.3):
    //   sqrt(n) * (theta_hat - theta) -->d N(0, Sigma)
    //
    // Assumption 2 (Consistent Variance Estimation, Section 2.3):
    //   Sigma_hat -->p Sigma
    //
    // Parameter vector structure (Section 2.2):
    //   theta = (nu_2, ..., nu_{t0-1}, delta_{t0}, ..., delta_T)'
    //
    // Covariance matrix partitioned structure:
    //   Sigma = | Sigma_nu        Sigma_nu,delta |
    //           | Sigma_delta,nu  Sigma_delta    |
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
    // Step 7: Bias Bound Constant Kappa (Section 3.2, Proposition 1)
    // ========================================
    // Under conditional extrapolation (Assumption 3), if S_pre <= M, then
    // the bias of the DID estimator is bounded:
    //
    //   |tau_bar - delta_bar| <= kappa * S_pre
    //
    // Kappa formula (derived via Hoelder's inequality):
    //
    //   kappa = ((1/T_post) * sum_{t=1}^{T_post} t^q)^{1/q}
    //
    // where q is the Hoelder conjugate exponent: 1/p + 1/q = 1
    //
    // Special cases:
    //   p = 1:   kappa = T_post (maximum time weight)
    //   p = 2:   kappa = sqrt((T_post+1)(2*T_post+1)/6) (closed-form)
    //   p = inf: kappa = (T_post+1)/2 (average time weight)
    //
    // Note: In overall mode (nubar), kappa = 1 regardless of p, because
    // cumulative violations directly bound bias without time-weight amplification.
    real scalar kappa
    kappa = _pretest_kappa(T_post_m, p_m, overall_flag)
    
    st_numscalar("__pretest_kappa", kappa)
    
    // ========================================
    // Step 8: Critical Value via Monte Carlo (Section 5.1, Appendix F)
    // ========================================
    // The critical value f(alpha, Sigma) is defined as (Section 5.1):
    //
    //   f(alpha, Sigma) = inf{c : P(psi(Z) >= c) <= alpha}  where Z ~ N(0, Sigma)
    //
    // The psi function captures joint statistical uncertainty (Section 5.1):
    //
    //   psi(x) = |mean(x_post)| + kappa * ((1/(T_pre-1)) sum |x_pre|^p)^{1/p}
    //
    // Monte Carlo simulation algorithm (Appendix F):
    //   1. Draw Z^{(s)} ~ N(0, Sigma) via Cholesky decomposition: Z = L * W, W ~ N(0, I)
    //   2. Compute psi^{(s)} = psi(Z^{(s)}) for each simulation s = 1, ..., S
    //   3. Return the ceiling((1-alpha)*S)-th order statistic as f(alpha, Sigma)
    //
    // Note: Regularization Sigma + epsilon*I ensures positive definiteness for Cholesky.
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
    // Step 9: Pre-Test Evaluation (Section 4.2, Theorem 1)
    // ========================================
    // The preliminary test evaluates whether extrapolation is justified:
    //
    //   phi = 1{S_pre > M}
    //
    // Theorem 1 (Asymptotic Consistency): Under Assumption 1, this test is
    // asymptotically consistent for separation s_n = omega(n^{-1/2}), meaning
    // both Type I and Type II errors vanish as n -> infinity when the true
    // severity is bounded away from the threshold M.
    //
    // Interpretation:
    //   phi = 0 (PASS): S_pre <= M, extrapolation condition satisfied,
    //                   proceed with conditionally valid CI (Theorem 2)
    //   phi = 1 (FAIL): S_pre > M, extrapolation not justified,
    //                   DID analysis requires alternative assumptions
    //   phi = .       : Computation failed due to data issues
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
    // Step 10: Conditionally Valid Confidence Interval (Theorem 2, Section 5.1)
    // ========================================
    // Theorem 2 (Conditional Validity): Under Assumptions 1-3 and the
    // well-separated null hypothesis (s_n = omega(n^{-1/2})), the following
    // confidence interval achieves asymptotic conditional coverage:
    //
    //   lim inf P(tau_bar in CI | phi = 0) >= 1 - alpha
    //
    // Iterative mode (Section 5.1):
    //   CI = delta_bar +/- {kappa * S_pre + f(alpha, Sigma) / sqrt(n)}
    //        |______________|   |__________________________|
    //         worst-case bias    statistical uncertainty
    //
    // Overall mode (Appendix C):
    //   CI = delta_bar +/- {S_pre + f(alpha, Sigma) / sqrt(n)}
    //        Note: kappa = 1 for overall violations (cumulative form)
    //
    // Width scaling behavior (Corollary 1): O(max{S_pre, n^{-1/2}})
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
    // Step 11: Parameter Vector Theta (Section 2.2)
    // ========================================
    // The full parameter vector for asymptotic analysis:
    //
    //   theta = (nu_2, ..., nu_{t0-1}, delta_{t0}, ..., delta_T)'
    //
    // Dimension: dim(theta) = (T_pre - 1) + T_post = T - 1
    //
    // In overall mode (Appendix C), replace nu with nubar:
    //   theta^{Delta} = (nubar_2, ..., nubar_{t0-1}, delta_{t0}, ..., delta_T)'
    real colvector theta_vec
    
    theta_vec = nu_vec \ delta_vec
    st_matrix("__pretest_theta", theta_vec)
    
    // ========================================
    // Step 12: Severity Standard Error (Delta Method)
    // ========================================
    // The Delta Method provides an approximation for the variance of S_pre:
    //
    //   Var(S_pre) approx grad_g' * Sigma_nu * grad_g
    //   SE(S_pre) = sqrt(Var(S_pre))
    //
    // where grad_g = partial S_pre / partial nu is the gradient vector:
    //
    //   dg/dnu_t = |nu_t|^{p-1} * sign(nu_t) / [n_pre * S_pre^{p-1}]
    //
    // Note: The Delta Method requires differentiability of S_pre with respect to nu.
    // For p = infinity (L-inf norm), the max function is non-differentiable,
    // so SE(S_pre) returns missing. Bootstrap methods may be used instead.
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
    // Step 13: Conventional DID Confidence Interval (for comparison)
    // ========================================
    // Standard CI assuming exact parallel trends (classical DID inference):
    //
    //   CI_conv = delta_bar +/- z_{1-alpha/2} * SE(delta_bar)
    //
    // Variance estimation:
    //   Var(delta_bar) = (1/T_post^2) * 1' * Sigma_delta * 1 / n
    //   SE(delta_bar) = sqrt(Var(delta_bar))
    //
    // Important: This conventional CI is provided for comparison only.
    // It is only valid under exact parallel trends (S_pre = S_post = 0).
    // Under violations of parallel trends, use the conditionally valid CI
    // from Step 10, which accounts for worst-case bias via kappa * S_pre.
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
