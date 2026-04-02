*! _pretest_ci.mata
*! Conditionally Valid Confidence Intervals for the Pre-Test DID Framework
*!
*! Description:
*!   Implements confidence intervals for the Average Treatment Effect on the
*!   Treated (ATT) that are asymptotically valid conditional on passing the
*!   pre-test for parallel trends violations. This approach resolves the 
*!   conditional coverage distortion identified by Roth (2022).
*!
*! Theoretical Foundation:
*!   Under the Conditional Extrapolation Assumption (Assumption 3), if the
*!   pre-treatment severity S_pre <= M, then post-treatment violations are
*!   bounded by S_pre. The CI accounts for both bias and sampling uncertainty.
*!
*! Conditional Coverage Property (Theorem 2):
*!   Under Assumptions 1-3 with well-separated null (s_n = omega(n^{-1/2})):
*!     lim inf_{n -> inf} P(tau_bar in CI | phi = 0) >= 1 - alpha
*!
*! CI Construction:
*!   Iterative mode: CI = delta_bar +/- {kappa * S_pre + f(alpha, Sigma)/sqrt(n)}
*!   Overall mode:   CI = delta_bar +/- {S_pre + f(alpha, Sigma)/sqrt(n)}
*!
*! Width Properties (Corollary 1):
*!   - If S_pre = O(n^{-1/2}): width = O(n^{-1/2}), same rate as conventional CI
*!   - If S_pre = O(1): width = O(1), interval does not shrink with sample size
*!
*! Reference:
*!   Mikhaeil, J. M. and C. Harshaw (2025). "In Defense of the Pre-Test: Valid
*!   Inference when Testing Violations of Parallel Trends for Difference-in-
*!   Differences." arXiv preprint arXiv:2510.26470.
*!   See Section 5.1, Theorem 2; Appendix C for overall mode.

version 17.0

mata:
mata set matastrict on

// ============================================================================
// ITERATIVE MODE CONFIDENCE INTERVAL
// ============================================================================

/**
 * @brief Construct conditionally valid CI for iterative violation mode
 *
 * Implements the confidence interval from Theorem 2:
 *
 *   CI = [delta_bar_hat - h, delta_bar_hat + h]
 *
 * where the half-width h = kappa * S_pre_hat + f(alpha, Sigma_hat) / sqrt(n)
 *
 * The half-width consists of two components:
 *   1. Bias adjustment: kappa * S_pre_hat bounds the worst-case bias
 *      |tau_bar - delta_bar| <= kappa * S_pre (Proposition 1)
 *   2. Sampling uncertainty: f(alpha, Sigma)/sqrt(n) accounts for the
 *      statistical fluctuations in estimating delta_bar and S_pre
 *
 * The constant kappa = (T_post^{-1} * sum_{t=1}^{T_post} t^q)^{1/q} where
 * q satisfies 1/p + 1/q = 1 (Holder conjugate).
 *
 * @param delta_bar_hat  Average post-treatment DID estimate
 * @param S_pre_hat      Estimated pre-treatment severity
 * @param kappa          Bias bound constant (depends on T_post and p)
 * @param f_alpha        Critical value f(alpha, Sigma_hat) from Monte Carlo
 * @param n              Total sample size
 *
 * @return 1x2 row vector (ci_lower, ci_upper), or (., .) if inputs invalid
 *
 * @note Valid only conditional on phi = 0 (pre-test passed)
 * @see Mikhaeil & Harshaw (2025), Section 5.1, Theorem 2
 */
real rowvector _pretest_ci_iterative(real scalar delta_bar_hat,
                                      real scalar S_pre_hat,
                                      real scalar kappa,
                                      real scalar f_alpha,
                                      real scalar n)
{
    real scalar half_width, ci_lower, ci_upper
    
    // Validate inputs
    if (missing(delta_bar_hat) | missing(S_pre_hat) | missing(kappa) | 
        missing(f_alpha) | missing(n)) {
        return((., .))
    }
    
    if (n <= 0) {
        errprintf("Error: n must be positive\n")
        return((., .))
    }
    
    // Half-width: h = kappa * S_pre_hat + f(alpha, Sigma_hat) / sqrt(n)
    half_width = kappa * S_pre_hat + f_alpha / sqrt(n)
    
    ci_lower = delta_bar_hat - half_width
    ci_upper = delta_bar_hat + half_width
    
    return((ci_lower, ci_upper))
}


// ============================================================================
// OVERALL MODE CONFIDENCE INTERVAL
// ============================================================================

/**
 * @brief Construct conditionally valid CI for overall violation mode
 *
 * Implements the confidence interval from Appendix C.3:
 *
 *   CI = [delta_bar_hat - h, delta_bar_hat + h]
 *
 * where the half-width h = S_pre_hat^Delta + f^Delta(alpha, Sigma^Delta) / sqrt(n)
 *
 * Key difference from iterative mode: The overall mode uses the severity
 * measure based on cumulative (overall) violations nu_bar_t rather than
 * iterative violations nu_t. The bias bound simplifies to |tau_bar - delta_bar|
 * <= S_pre^Delta without the kappa multiplier (see Appendix C.1, Proposition).
 *
 * @param delta_bar_hat      Average post-treatment DID estimate
 * @param S_pre_hat_overall  Overall mode severity estimate S_pre^Delta
 * @param f_alpha_overall    Critical value f^Delta(alpha, Sigma^Delta)
 * @param n                  Total sample size
 *
 * @return 1x2 row vector (ci_lower, ci_upper), or (., .) if inputs invalid
 *
 * @note No kappa multiplier in overall mode (fundamental distinction)
 * @see Mikhaeil & Harshaw (2025), Appendix C.3
 */
real rowvector _pretest_ci_overall(real scalar delta_bar_hat,
                                    real scalar S_pre_hat_overall,
                                    real scalar f_alpha_overall,
                                    real scalar n)
{
    real scalar half_width, ci_lower, ci_upper
    
    // Validate inputs
    if (missing(delta_bar_hat) | missing(S_pre_hat_overall) | 
        missing(f_alpha_overall) | missing(n)) {
        return((., .))
    }
    
    if (n <= 0) {
        errprintf("Error: n must be positive\n")
        return((., .))
    }
    
    // Half-width: h = S_pre^Delta + f^Delta(alpha, Sigma^Delta) / sqrt(n)
    // Note: No kappa multiplier in overall mode (see Appendix C.1)
    half_width = S_pre_hat_overall + f_alpha_overall / sqrt(n)
    
    ci_lower = delta_bar_hat - half_width
    ci_upper = delta_bar_hat + half_width
    
    return((ci_lower, ci_upper))
}


// ============================================================================
// CONVENTIONAL DID CONFIDENCE INTERVAL
// ============================================================================

/**
 * @brief Construct conventional DID confidence interval
 *
 * Standard CI under the assumption of exact parallel trends:
 *
 *   CI = [delta_bar - z_{1-alpha/2} * SE, delta_bar + z_{1-alpha/2} * SE]
 *
 * This interval assumes parallel trends hold exactly (nu_t = 0 for all t),
 * so that tau_bar = delta_bar. Provided for comparison with the conditionally
 * valid intervals which account for potential parallel trends violations.
 *
 * Warning: This interval may suffer from coverage distortion when used after
 * preliminary testing for parallel trends (Roth, 2022).
 *
 * @param delta_bar  DID estimate (average post-treatment ATT under PT)
 * @param se_delta   Standard error of the DID estimate
 * @param alpha      Significance level (default: 0.05)
 *
 * @return 1x2 row vector (ci_lower, ci_upper), or (., .) if inputs invalid
 *
 * @note Uses standard normal critical value z_{1-alpha/2}
 */
real rowvector _pretest_ci_conventional(real scalar delta_bar,
                                         real scalar se_delta,
                                         | real scalar alpha)
{
    real scalar z, ci_lower, ci_upper
    
    // Default alpha
    if (args() < 3 | missing(alpha)) {
        alpha = 0.05
    }
    
    // Validate inputs
    if (missing(delta_bar) | missing(se_delta)) {
        return((., .))
    }
    
    if (se_delta < 0) {
        errprintf("Error: se_delta must be non-negative\n")
        return((., .))
    }
    
    if (alpha <= 0 | alpha >= 1) {
        errprintf("Error: alpha must be in (0, 1)\n")
        return((., .))
    }
    
    // Standard normal critical value: z_{1-alpha/2}
    z = invnormal(1 - alpha/2)
    
    ci_lower = delta_bar - z * se_delta
    ci_upper = delta_bar + z * se_delta
    
    return((ci_lower, ci_upper))
}


// ============================================================================
// PRETEST CHECK FUNCTION
// ============================================================================

/**
 * @brief Evaluate the pre-test for extrapolation condition
 *
 * Implements the preliminary test from Section 4.2:
 *
 *   phi = 1{S_pre_hat > M}
 *
 * The test determines whether extrapolation of parallel trends violations
 * from pre- to post-treatment is justified. Under the Conditional
 * Extrapolation Assumption, if S_pre <= M, then S_post <= S_pre.
 *
 * Theorem 1 establishes that this test is asymptotically consistent for
 * separation s_n = omega(n^{-1/2}), meaning both Type I and Type II errors
 * vanish as sample size grows for well-separated hypotheses.
 *
 * @param S_pre_hat  Estimated pre-treatment severity
 * @param M          Acceptable threshold for extrapolation
 *
 * @return 1 if pre-test passes (phi = 0, S_pre_hat <= M, extrapolation valid)
 *         0 if pre-test fails (phi = 1, S_pre_hat > M, extrapolation invalid)
 *         . if inputs are missing
 *
 * @note Boundary case: S_pre_hat = M is classified as PASS
 * @see Mikhaeil & Harshaw (2025), Section 4.2, Theorem 1
 */
real scalar _pretest_check(real scalar S_pre_hat, real scalar M)
{
    // Test statistic: phi = 1{S_pre_hat > M}
    // Return convention: 1 = pass (phi=0), 0 = fail (phi=1)
    
    if (missing(S_pre_hat) | missing(M)) {
        return(.)
    }
    
    // S_pre_hat <= M: extrapolation condition holds, return pass (1)
    // S_pre_hat > M: extrapolation condition fails, return fail (0)
    if (S_pre_hat <= M) {
        return(1)
    } else {
        return(0)
    }
}


// ============================================================================
// MAIN CONFIDENCE INTERVAL FUNCTION (WITH CONDITIONAL CHECK)
// ============================================================================

/**
 * @brief Main interface for conditionally valid CI with pre-test check
 *
 * This is the primary function for computing confidence intervals under the
 * conditional extrapolation framework. It first evaluates the pre-test to
 * determine if extrapolation is justified, then computes the appropriate CI.
 *
 * Workflow:
 *   1. Evaluate pre-test: phi = 1{S_pre_hat > M}
 *   2. If phi = 1 (fail): return missing CI bounds
 *   3. If phi = 0 (pass): compute CI based on specified mode
 *
 * @param delta_bar_hat  Average post-treatment DID estimate
 * @param S_pre_hat      Estimated pre-treatment severity
 * @param M              Acceptable threshold for extrapolation
 * @param kappa          Bias bound constant (used only in iterative mode)
 * @param f_alpha        Critical value from Monte Carlo simulation
 * @param n              Total sample size
 * @param mode           Violation mode: "iterative" or "overall"
 *
 * @return 1x3 row vector: (ci_lower, ci_upper, pretest_passed)
 *         - pretest_passed: 1 if passed, 0 if failed, . if error
 *         - CI bounds are missing (.) if pre-test fails or on error
 *
 * @see Mikhaeil & Harshaw (2025), Section 5.1, Theorem 2
 */
real rowvector _pretest_ci(real scalar delta_bar_hat,
                            real scalar S_pre_hat,
                            real scalar M,
                            real scalar kappa,
                            real scalar f_alpha,
                            real scalar n,
                            string scalar mode)
{
    real scalar pretest_passed, ci_lower, ci_upper
    real rowvector ci_result
    
    // Step 1: Evaluate pre-test phi = 1{S_pre_hat > M}
    pretest_passed = _pretest_check(S_pre_hat, M)
    
    // Step 2: If pre-test fails (phi = 1), extrapolation not justified
    if (pretest_passed == 0) {
        return((., ., 0))
    }
    
    // Handle missing pre-test result
    if (missing(pretest_passed)) {
        return((., ., .))
    }
    
    // Step 3: Pre-test passed (phi = 0), compute CI based on mode
    if (mode == "iterative") {
        ci_result = _pretest_ci_iterative(delta_bar_hat, S_pre_hat, 
                                          kappa, f_alpha, n)
        ci_lower = ci_result[1]
        ci_upper = ci_result[2]
    }
    else if (mode == "overall") {
        // Overall mode: no kappa multiplier (Appendix C)
        ci_result = _pretest_ci_overall(delta_bar_hat, S_pre_hat, 
                                        f_alpha, n)
        ci_lower = ci_result[1]
        ci_upper = ci_result[2]
    }
    else {
        errprintf("Error: Invalid mode '%s'. Use 'iterative' or 'overall'.\n", mode)
        return((., ., -1))
    }
    
    return((ci_lower, ci_upper, 1))
}


// ============================================================================
// HALF-WIDTH CALCULATION FUNCTIONS
// ============================================================================

/**
 * @brief Calculate CI half-width for iterative mode
 *
 * Computes the half-width of the conditionally valid CI:
 *   h = kappa * S_pre_hat + f(alpha, Sigma_hat) / sqrt(n)
 *
 * @param S_pre_hat  Estimated pre-treatment severity
 * @param kappa      Bias bound constant from Proposition 1
 * @param f_alpha    Critical value from Monte Carlo simulation
 * @param n          Total sample size
 *
 * @return Half-width h, or missing (.) if n <= 0
 */
real scalar _pretest_ci_halfwidth_iterative(real scalar S_pre_hat,
                                             real scalar kappa,
                                             real scalar f_alpha,
                                             real scalar n)
{
    if (n <= 0 | missing(n)) {
        return(.)
    }
    
    return(kappa * S_pre_hat + f_alpha / sqrt(n))
}


/**
 * @brief Calculate CI half-width for overall mode
 *
 * Computes the half-width of the conditionally valid CI:
 *   h = S_pre_hat^Delta + f^Delta(alpha, Sigma^Delta) / sqrt(n)
 *
 * Note: No kappa multiplier in overall mode (Appendix C.1).
 *
 * @param S_pre_hat_overall  Overall mode severity estimate S_pre^Delta
 * @param f_alpha_overall    Critical value f^Delta(alpha, Sigma^Delta)
 * @param n                  Total sample size
 *
 * @return Half-width h, or missing (.) if n <= 0
 */
real scalar _pretest_ci_halfwidth_overall(real scalar S_pre_hat_overall,
                                           real scalar f_alpha_overall,
                                           real scalar n)
{
    if (n <= 0 | missing(n)) {
        return(.)
    }
    
    // No kappa multiplier in overall mode
    return(S_pre_hat_overall + f_alpha_overall / sqrt(n))
}


// ============================================================================
// COMPREHENSIVE CI FUNCTION WITH ALL RESULTS
// ============================================================================

/**
 * @brief Comprehensive CI calculation returning all results
 *
 * Computes both conditionally valid and conventional CIs for comparison.
 * This function facilitates the comparison illustrated in Figure 1 of the
 * paper, where both interval types are shown side-by-side.
 *
 * @param delta_bar_hat  Average post-treatment DID estimate
 * @param se_delta       Standard error of the DID estimate
 * @param S_pre_hat      Estimated pre-treatment severity
 * @param M              Acceptable threshold for extrapolation
 * @param kappa          Bias bound constant (iterative mode only)
 * @param f_alpha        Critical value from Monte Carlo simulation
 * @param n              Total sample size
 * @param alpha          Significance level (default: 0.05)
 * @param mode           Violation mode: "iterative" or "overall"
 *
 * @return 1x8 row vector:
 *         [1-2] Conditionally valid CI bounds (lower, upper)
 *         [3]   Conditionally valid CI half-width
 *         [4]   Pre-test result (1=passed, 0=failed)
 *         [5-6] Conventional CI bounds (lower, upper)
 *         [7]   Conventional CI half-width
 *         [8]   Status code (0=success, 1=pretest failed, -1=error)
 */
real rowvector _pretest_ci_results(real scalar delta_bar_hat,
                                    real scalar se_delta,
                                    real scalar S_pre_hat,
                                    real scalar M,
                                    real scalar kappa,
                                    real scalar f_alpha,
                                    real scalar n,
                                    real scalar alpha,
                                    string scalar mode)
{
    real scalar pretest_passed, ci_lower, ci_upper, half_width
    real scalar conv_ci_lower, conv_ci_upper, conv_half_width, status
    real rowvector ci_result, conv_result
    
    // Default alpha
    if (missing(alpha)) {
        alpha = 0.05
    }
    
    // Compute conventional CI (always computed for comparison)
    conv_result = _pretest_ci_conventional(delta_bar_hat, se_delta, alpha)
    conv_ci_lower = conv_result[1]
    conv_ci_upper = conv_result[2]
    conv_half_width = (conv_ci_upper - conv_ci_lower) / 2
    
    // Check pretest
    pretest_passed = _pretest_check(S_pre_hat, M)
    
    // If pretest failed
    if (pretest_passed == 0) {
        return((., ., ., 0, conv_ci_lower, conv_ci_upper, conv_half_width, 1))
    }
    
    // Compute conditionally valid CI
    if (mode == "iterative") {
        ci_result = _pretest_ci_iterative(delta_bar_hat, S_pre_hat, 
                                          kappa, f_alpha, n)
        half_width = _pretest_ci_halfwidth_iterative(S_pre_hat, kappa, f_alpha, n)
    }
    else if (mode == "overall") {
        ci_result = _pretest_ci_overall(delta_bar_hat, S_pre_hat, f_alpha, n)
        half_width = _pretest_ci_halfwidth_overall(S_pre_hat, f_alpha, n)
    }
    else {
        errprintf("Error: Invalid mode '%s'\n", mode)
        return((., ., ., ., conv_ci_lower, conv_ci_upper, conv_half_width, -1))
    }
    
    ci_lower = ci_result[1]
    ci_upper = ci_result[2]
    status = 0
    
    return((ci_lower, ci_upper, half_width, 1, 
            conv_ci_lower, conv_ci_upper, conv_half_width, status))
}


// ============================================================================
// WRAPPER FOR ADO FILE
// ============================================================================

/**
 * @brief Wrapper function for Stata ado file interface
 *
 * Designed to be called from Stata ado files via mata: statements.
 * Parses string arguments to real scalars and stores results in
 * Stata scalars for retrieval by the calling ado file.
 *
 * Stored Stata scalars:
 *   - pretest_passed: 1 if passed, 0 if failed
 *   - ci_lower: Lower bound of CI (missing if pre-test failed)
 *   - ci_upper: Upper bound of CI (missing if pre-test failed)
 *
 * @param delta_bar_str  DID estimate as string
 * @param S_pre_str      Severity estimate as string
 * @param M_str          Threshold as string
 * @param kappa_str      Bias bound constant as string
 * @param f_alpha_str    Critical value as string
 * @param n_str          Sample size as string
 * @param alpha_str      Significance level as string
 * @param mode           Violation mode: "iterative" or "overall"
 */
void _pretest_ci_wrapper(string scalar delta_bar_str,
                          string scalar S_pre_str,
                          string scalar M_str,
                          string scalar kappa_str,
                          string scalar f_alpha_str,
                          string scalar n_str,
                          string scalar alpha_str,
                          string scalar mode)
{
    real scalar delta_bar_hat, S_pre_hat, M, kappa, f_alpha, n, alpha
    real scalar pretest_passed, ci_lower, ci_upper
    real rowvector result
    
    // Parse string arguments to real scalars
    delta_bar_hat = strtoreal(delta_bar_str)
    S_pre_hat = strtoreal(S_pre_str)
    M = strtoreal(M_str)
    kappa = strtoreal(kappa_str)
    f_alpha = strtoreal(f_alpha_str)
    n = strtoreal(n_str)
    alpha = strtoreal(alpha_str)
    
    // Compute CI
    result = _pretest_ci(delta_bar_hat, S_pre_hat, M, kappa, f_alpha, n, mode)
    
    ci_lower = result[1]
    ci_upper = result[2]
    pretest_passed = result[3]
    
    // Store results in Stata
    st_numscalar("pretest_passed", pretest_passed)
    st_numscalar("ci_lower", ci_lower)
    st_numscalar("ci_upper", ci_upper)
}


end
