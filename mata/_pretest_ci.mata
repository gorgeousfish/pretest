*! _pretest_ci.mata v0.1.0
*! Conditionally Valid Confidence Intervals for Pre-Test Framework
*!
*! Description:
*!   Implements confidence interval construction for the conditional
*!   extrapolation pre-test framework in difference-in-differences designs.
*!   These CIs are valid conditional on passing the pre-test (phi = 0).
*!
*! Main Result (Section 5.1, Theorem 2):
*!   Under Assumptions 1-3, the CI:
*!     delta_bar +/- {kappa * S_pre + f(alpha, Sigma) / sqrt(n)}
*!   achieves conditional coverage:
*!     P(tau_bar in CI | phi = 0) >= 1 - alpha  asymptotically
*!
*! Contents:
*!   Iterative mode: CI = delta_bar +/- {kappa*S_pre + f(alpha,Sigma)/sqrt(n)}
*!   Overall mode:   CI = delta_bar +/- {S_pre + f(alpha,Sigma)/sqrt(n)}
*!   Conventional:   CI = delta_bar +/- z_{1-alpha/2} * SE(delta_bar)
*!
*! Reference:
*!   Mikhaeil, J.M. and C. Harshaw. 2025. "In Defense of the Pre-Test: Valid
*!   Inference when Testing Violations of Parallel Trends for Difference-in-
*!   Differences." arXiv preprint arXiv:2510.26470.
*!   URL: https://arxiv.org/abs/2510.26470

version 17.0

mata:

// ============================================================================
// ITERATIVE MODE CONFIDENCE INTERVAL
// ============================================================================

/**
 * @function _pretest_ci_iterative
 * @brief Construct conditionally valid CI for iterative violation mode
 *
 * Implements the confidence interval from Section 5.1, Theorem 2:
 *
 *   CI = [delta_bar - h, delta_bar + h]
 *
 * where the half-width h = kappa * S_pre + f(alpha, Sigma) / sqrt(n)
 *
 * The two terms in h correspond to:
 *   - kappa * S_pre: Worst-case bias bound (Proposition 1)
 *   - f(alpha, Sigma) / sqrt(n): Statistical uncertainty from estimation
 *
 * @param delta_bar_hat Average post-treatment DID estimate delta_bar_hat
 * @param S_pre_hat     Estimated pre-treatment severity S_pre_hat
 * @param kappa         Bias bound constant from Proposition 1
 * @param f_alpha       Critical value f(alpha, Sigma) from Monte Carlo
 * @param n             Total sample size
 *
 * @return Row vector (ci_lower, ci_upper), or (., .) if inputs invalid
 *
 * @note This CI is only valid conditional on phi = 0 (pre-test passed)
 *
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
    
    // Formula: half_width = κ·Ŝ_pre + f(α,Σ̂)/√n
    half_width = kappa * S_pre_hat + f_alpha / sqrt(n)
    
    ci_lower = delta_bar_hat - half_width
    ci_upper = delta_bar_hat + half_width
    
    return((ci_lower, ci_upper))
}


// ============================================================================
// OVERALL MODE CONFIDENCE INTERVAL
// ============================================================================

/**
 * @function _pretest_ci_overall
 * @brief Calculate conditionally valid CI in overall mode
 *
 * Implements the confidence interval from Appendix C:
 *
 *   CI = delta_bar +/- {S_pre + f(alpha, Sigma) / sqrt(n)}
 *
 * IMPORTANT: Overall mode does NOT include kappa multiplier.
 *
 * @param delta_bar_hat     Average post-treatment DID estimate
 * @param S_pre_hat_overall Overall mode severity estimate
 * @param f_alpha_overall   Overall mode critical value
 * @param n                 Sample size
 *
 * @return Row vector (1 x 2) containing (ci_lower, ci_upper)
 *
 * @note No kappa multiplier in overall mode (key difference from iterative)
 *
 * @see Mikhaeil & Harshaw (2025), Appendix C
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
    
    // Formula: half_width = Ŝ_pre^Δ + f^Δ(α,Σ̂^Δ)/√n
    // NOTE: NO κ multiplier in overall mode!
    half_width = S_pre_hat_overall + f_alpha_overall / sqrt(n)
    
    ci_lower = delta_bar_hat - half_width
    ci_upper = delta_bar_hat + half_width
    
    return((ci_lower, ci_upper))
}


// ============================================================================
// CONVENTIONAL DID CONFIDENCE INTERVAL
// ============================================================================

/**
 * @function _pretest_ci_conventional
 * @brief Calculate conventional DID confidence interval
 *
 * Standard CI assuming parallel trends:
 *
 *   CI = delta_bar +/- z_{1-alpha/2} * SE(delta_bar)
 *
 * Provided for comparison with conditionally valid intervals.
 *
 * @param delta_bar ATT estimate
 * @param se_delta  Standard error of ATT
 * @param alpha     Significance level (default: 0.05)
 *
 * @return Row vector (1 x 2) containing (ci_lower, ci_upper)
 *
 * @note Uses standard normal critical value (z_{0.975} = 1.96 for alpha = 0.05)
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
    
    // Standard normal critical value: z_{1-α/2}
    z = invnormal(1 - alpha/2)
    
    ci_lower = delta_bar - z * se_delta
    ci_upper = delta_bar + z * se_delta
    
    return((ci_lower, ci_upper))
}


// ============================================================================
// PRETEST CHECK FUNCTION
// ============================================================================

/**
 * @function _pretest_check
 * @brief Check if pre-test passes (extrapolation justified)
 *
 * Evaluates the pre-test condition:
 *
 *   phi = 1{S_pre > M}
 *
 * @param S_pre_hat Estimated pre-treatment severity
 * @param M         Acceptable threshold
 *
 * @return 1 if pre-test passes (phi = 0, S_pre <= M),
 *         0 if pre-test fails (phi = 1, S_pre > M)
 *
 * @note S_pre = M is considered PASS (uses <= comparison)
 *
 * @see Mikhaeil & Harshaw (2025), Section 4.2
 */
real scalar _pretest_check(real scalar S_pre_hat, real scalar M)
{
    // φ = 1{Ŝ_pre > M}
    // Return: 1 if pass (φ=0), 0 if fail (φ=1)
    
    if (missing(S_pre_hat) | missing(M)) {
        return(.)
    }
    
    // S_pre ≤ M: pass (return 1)
    // S_pre > M: fail (return 0)
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
 * @function _pretest_ci
 * @brief Main CI function with pre-test condition check
 *
 * Primary interface for computing conditionally valid confidence intervals.
 * Checks pre-test condition first, then computes appropriate CI if passed.
 *
 * @param delta_bar_hat ATT estimate
 * @param S_pre_hat     Severity estimate
 * @param M             Pre-test threshold
 * @param kappa         Time weight constant (iterative mode)
 * @param f_alpha       Adjusted critical value
 * @param n             Sample size
 * @param mode          "iterative" or "overall"
 *
 * @return Row vector (1 x 3): (ci_lower, ci_upper, pretest_passed)
 *         - pretest_passed: 1 if passed, 0 if failed
 *         - CI bounds are missing (.) if pre-test fails
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
    
    // Check pretest: φ = 1{Ŝ_pre > M}
    pretest_passed = _pretest_check(S_pre_hat, M)
    
    // If pretest failed (φ = 1), return missing values
    if (pretest_passed == 0) {
        return((., ., 0))
    }
    
    // Handle missing pretest result
    if (missing(pretest_passed)) {
        return((., ., .))
    }
    
    // Pretest passed (φ = 0), compute CI based on mode
    if (mode == "iterative") {
        ci_result = _pretest_ci_iterative(delta_bar_hat, S_pre_hat, 
                                          kappa, f_alpha, n)
        ci_lower = ci_result[1]
        ci_upper = ci_result[2]
    }
    else if (mode == "overall") {
        // Overall mode: no κ multiplier
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
 * @function _pretest_ci_halfwidth_iterative
 * @brief Calculate CI half-width for iterative mode
 *
 * @param S_pre_hat Severity estimate
 * @param kappa     Time weight constant
 * @param f_alpha   Critical value
 * @param n         Sample size
 *
 * @return Half-width: kappa * S_pre + f_alpha / sqrt(n)
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
 * @function _pretest_ci_halfwidth_overall
 * @brief Calculate CI half-width for overall mode
 *
 * @param S_pre_hat_overall Overall severity estimate
 * @param f_alpha_overall   Overall mode critical value
 * @param n                 Sample size
 *
 * @return Half-width: S_pre + f_alpha / sqrt(n) (no kappa)
 */
real scalar _pretest_ci_halfwidth_overall(real scalar S_pre_hat_overall,
                                           real scalar f_alpha_overall,
                                           real scalar n)
{
    if (n <= 0 | missing(n)) {
        return(.)
    }
    
    // NO κ multiplier in overall mode
    return(S_pre_hat_overall + f_alpha_overall / sqrt(n))
}


// ============================================================================
// COMPREHENSIVE CI FUNCTION WITH ALL RESULTS
// ============================================================================

/**
 * @function _pretest_ci_results
 * @brief Comprehensive CI calculation returning all results
 *
 * Computes both conditionally valid and conventional CIs for comparison.
 *
 * @param delta_bar_hat ATT estimate
 * @param se_delta      Standard error of ATT
 * @param S_pre_hat     Severity estimate
 * @param M             Pre-test threshold
 * @param kappa         Time weight constant
 * @param f_alpha       Adjusted critical value
 * @param n             Sample size
 * @param alpha         Significance level (default: 0.05)
 * @param mode          "iterative" or "overall"
 *
 * @return Row vector (1 x 8):
 *         [1-2] Conditional CI bounds
 *         [3]   Conditional CI half-width
 *         [4]   Pre-test passed (1/0)
 *         [5-6] Conventional CI bounds
 *         [7]   Conventional CI half-width
 *         [8]   Status (0=success, 1=pretest failed, -1=error)
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
 * @function _pretest_ci_wrapper
 * @brief Wrapper function for ado file interface
 *
 * Designed to be called from Stata ado files. Parses string arguments
 * and stores results in Stata scalars.
 *
 * @param delta_bar_str ATT estimate (string)
 * @param S_pre_str     Severity estimate (string)
 * @param M_str         Pre-test threshold (string)
 * @param kappa_str     Time weight constant (string)
 * @param f_alpha_str   Critical value (string)
 * @param n_str         Sample size (string)
 * @param alpha_str     Significance level (string)
 * @param mode          "iterative" or "overall"
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
