*! _pretest_violations.mata v0.1.0
*! Parallel Trends Violations Estimation Module
*!
*! Description:
*!   Implements estimators for violations of parallel trends in the conditional
*!   extrapolation pre-test framework for difference-in-differences designs.
*!   Provides both iterative (period-to-period) and overall (cumulative)
*!   violation measures as defined in the theoretical framework.
*!
*! Mathematical Framework (Section 2.1):
*!
*!   Iterative Violation (successive period differences):
*!     nu_t = E[Y^{(0)}_t - Y^{(0)}_{t-1} | D=1] - E[Y^{(0)}_t - Y^{(0)}_{t-1} | D=0]
*!
*!   Overall Violation (cumulative deviation from parallel trends):
*!     nubar_t = sum_{s=2}^t nu_s
*!
*!   Relationship: nubar_t = sum nu_s  and  nu_t = nubar_t - nubar_{t-1}
*!
*! Severity Measure (Section 3.1):
*!   S_pre = ((1/(T_pre-1)) * sum_{t=2}^{t0-1} |nu_t|^p)^{1/p}
*!
*!   This measures the magnitude of pre-treatment parallel trends violations.
*!   Under the conditional extrapolation assumption (Assumption 3):
*!     If S_pre <= M, then S_post <= S_pre
*!
*! Pre-Test Function (Section 4.2, Theorem 1):
*!   phi = 1{S_pre > M}
*!
*!   phi = 0 (PASS): Extrapolation condition satisfied, proceed with inference
*!   phi = 1 (FAIL): Violations too severe, extrapolation not justified
*!
*! Reference:
*!   Mikhaeil, J. M. and C. Harshaw. 2025. In Defense of the Pre-Test: Valid
*!   Inference when Testing Violations of Parallel Trends for Difference-in-
*!   Differences. arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470

version 17.0

mata:

// ============================================================================
// ITERATIVE VIOLATIONS
// ============================================================================

/**
 * @function _pretest_nu_t
 * @brief Compute single-period iterative violation of parallel trends
 *
 * Implements the iterative violation estimator from Section 2.1:
 *
 *   nu_hat_t = [Y_bar_{1,t} - Y_bar_{1,t-1}] - [Y_bar_{0,t} - Y_bar_{0,t-1}]
 *
 * where Y_bar_{d,t} denotes the sample mean of outcomes for group D=d at time t.
 *
 * @param Y    Outcome vector (n x 1)
 * @param D    Treatment indicator vector (n x 1), values in {0, 1}
 * @param time Time period vector (n x 1)
 * @param t    Target time period for violation calculation
 * @param t0   Treatment time (first post-treatment period)
 *
 * @return Estimated iterative violation nu_hat_t, or missing (.) if
 *         insufficient data in any required group-time cell
 *
 * @note Valid for t in {2, ..., t0-1} (pre-treatment periods only)
 *
 * @see Mikhaeil & Harshaw (2025), Section 2.1, Equation (184)
 */
real scalar _pretest_nu_t(real colvector Y, real colvector D, 
                          real colvector time, real scalar t, real scalar t0)
{
    real scalar Y1t, Y1t_1, Y0t, Y0t_1
    
    // Validate t is in valid range for iterative violations
    if (t < 2 | t >= t0) {
        return(.)
    }
    
    // Get group means at time t and t-1
    Y1t   = _pretest_group_mean(Y, D, time, t, 1)
    Y1t_1 = _pretest_group_mean(Y, D, time, t - 1, 1)
    Y0t   = _pretest_group_mean(Y, D, time, t, 0)
    Y0t_1 = _pretest_group_mean(Y, D, time, t - 1, 0)
    
    // Check for missing values
    if (missing(Y1t) | missing(Y1t_1) | missing(Y0t) | missing(Y0t_1)) {
        return(.)
    }
    
    // nu_hat_t = (Y_bar_{1,t} - Y_bar_{1,t-1}) - (Y_bar_{0,t} - Y_bar_{0,t-1})
    return((Y1t - Y1t_1) - (Y0t - Y0t_1))
}


/**
 * @function _pretest_nu_vector
 * @brief Compute full vector of iterative violations
 *
 * Constructs the vector of pre-treatment iterative violations:
 *
 *   nu = (nu_hat_2, nu_hat_3, ..., nu_hat_{t0-1})'
 *
 * This vector forms part of the parameter vector theta in the asymptotic
 * framework (Section 2.2).
 *
 * @param Y    Outcome vector (n x 1)
 * @param D    Treatment indicator vector (n x 1)
 * @param time Time period vector (n x 1)
 * @param t0   Treatment time (first post-treatment period)
 *
 * @return Column vector of dimension (T_pre - 1) x 1, where T_pre = t0 - 1
 *
 * @see Mikhaeil & Harshaw (2025), Section 2.1
 */
real colvector _pretest_nu_vector(real colvector Y, real colvector D,
                                  real colvector time, real scalar t0)
{
    real scalar T_pre, t, idx
    real colvector nu
    
    // T_pre = t0 - 2 (number of iterative violations)
    T_pre = t0 - 2
    
    if (T_pre < 1) {
        // Edge case: no pre-treatment violations possible
        return(J(0, 1, .))
    }
    
    nu = J(T_pre, 1, .)
    
    // Compute nu_hat_t for t = 2, 3, ..., t0-1
    idx = 1
    for (t = 2; t <= t0 - 1; t++) {
        nu[idx] = _pretest_nu_t(Y, D, time, t, t0)
        idx++
    }
    
    return(nu)
}


// ============================================================================
// OVERALL VIOLATIONS
// ============================================================================

/**
 * @function _pretest_nubar_t
 * @brief Compute single-period overall violation via cumulative summation
 *
 * Implements the overall violation from Section 2.1:
 *
 *   nubar_t = sum_{s=2}^t nu_s
 *
 * Overall violations are computed as cumulative sums of iterative violations,
 * ensuring theoretical consistency with the framework where nubar_t represents
 * cumulative deviation from parallel trends relative to the baseline period.
 *
 * @param nu_vec Iterative violation vector (nu_hat_2, ..., nu_hat_{t0-1})'
 * @param k      Index (1-based): nu_vec[k] corresponds to nu_hat_{k+1}
 *
 * @return Overall violation nubar_hat_{k+1} = sum(nu_vec[1..k])
 *
 * @example
 *   nu_vec = (0.5, 1.0, -0.3)', k = 2
 *   Returns nubar_hat_3 = 0.5 + 1.0 = 1.5
 *
 * @see Mikhaeil & Harshaw (2025), Section 2.1, Equation (191)
 */
real scalar _pretest_nubar_t(real colvector nu_vec, real scalar k)
{
    if (k < 1 | k > rows(nu_vec)) {
        _error("Index k out of bounds in _pretest_nubar_t")
        return(.)
    }
    
    return(sum(nu_vec[1..k]))
}


/**
 * @function _pretest_nubar_vector
 * @brief Compute full vector of overall violations
 *
 * Constructs the vector of pre-treatment overall violations:
 *
 *   nubar = (nubar_hat_2, nubar_hat_3, ..., nubar_hat_{t0-1})'
 *
 * where nubar_t = sum_{s=2}^t nu_s represents cumulative deviation from
 * parallel trends. Computed as cumulative sums of iterative violations
 * for numerical consistency.
 *
 * @param Y    Outcome vector (n x 1)
 * @param D    Treatment indicator vector (n x 1)
 * @param time Time period vector (n x 1)
 * @param t0   Treatment time (first post-treatment period)
 *
 * @return Column vector of dimension (T_pre - 1) x 1
 *
 * @see Mikhaeil & Harshaw (2025), Section 2.1, Appendix D
 */
real colvector _pretest_nubar_vector(real colvector Y, real colvector D,
                                     real colvector time, real scalar t0)
{
    real colvector nu_vec
    
    if (t0 < 3) {
        return(J(0, 1, .))
    }
    
    // Step 1: Compute iterative violations (directly estimable)
    nu_vec = _pretest_nu_vector(Y, D, time, t0)
    
    // Step 2: Convert to overall violations via cumulative sum
    // runningsum((nu_2, nu_3, ...)) = (nu_bar_2, nu_bar_3, ...)
    return(runningsum(nu_vec))
}


// ============================================================================
// CONVERSION FUNCTIONS
// ============================================================================

/**
 * @function _pretest_nu_to_nubar
 * @brief Convert iterative violations to overall violations
 *
 * Applies cumulative summation to transform iterative violations into
 * overall violations:
 *
 *   nubar_t = sum_{s=2}^t nu_s
 *
 * @param nu Iterative violation vector (nu_hat_2, ..., nu_hat_{t0-1})'
 *
 * @return Overall violations (nubar_hat_2, ..., nubar_hat_{t0-1})'
 *
 * @example
 *   Input:  (0.5, 1.0, -0.3)'
 *   Output: (0.5, 1.5, 1.2)' (cumulative sums)
 *
 * @see Mikhaeil & Harshaw (2025), Section 2.1, Appendix D
 */
real colvector _pretest_nu_to_nubar(real colvector nu)
{
    if (rows(nu) < 1) {
        return(J(0, 1, .))
    }
    
    // Simple: runningsum((a,b,c)) = (a, a+b, a+b+c)
    return(runningsum(nu))
}


/**
 * @function _pretest_nubar_to_nu
 * @brief Convert overall violations to iterative violations
 *
 * Applies first-differencing to transform overall violations back to
 * iterative violations:
 *
 *   nu_t = nubar_t - nubar_{t-1}
 *
 * This is the inverse operation of _pretest_nu_to_nubar().
 *
 * @param nubar Overall violation vector (nubar_1, ..., nubar_{t0-1})'
 *              of dimension (t0-1) x 1
 *
 * @return Iterative violation vector (nu_2, ..., nu_{t0-1})'
 *         of dimension (T_pre) x 1, where T_pre = t0 - 2
 *
 * @note Uses forward differencing: nu[t-1] = nubar[t] - nubar[t-1]
 *
 * @see Mikhaeil & Harshaw (2025), Appendix D, Equation (191)
 */
real colvector _pretest_nubar_to_nu(real colvector nubar)
{
    real scalar t0_minus_1, T_pre, t
    real colvector nu
    
    t0_minus_1 = rows(nubar)
    
    if (t0_minus_1 < 2) {
        return(J(0, 1, .))
    }
    
    T_pre = t0_minus_1 - 1
    nu = J(T_pre, 1, .)
    
    // nu_t = nu_bar_t - nu_bar_{t-1} for t = 2, ..., t0-1
    // nu[1] = nu_2 = nu_bar_2 - nu_bar_1 = nubar[2] - nubar[1]
    // nu[t-1] = nu_t = nu_bar_t - nu_bar_{t-1} = nubar[t] - nubar[t-1]
    for (t = 2; t <= t0_minus_1; t++) {
        nu[t - 1] = nubar[t] - nubar[t - 1]
    }
    
    return(nu)
}


// ============================================================================
// SEVERITY MEASURE
// ============================================================================

/**
 * @function _pretest_severity
 * @brief Calculate the severity measure S_pre for pre-treatment violations
 *
 * Implements the severity measure from Section 3.1:
 *
 *   S_pre = ((1/(T_pre-1)) * sum_{t=2}^{t0-1} |nu_t|^p)^{1/p}
 *
 * This quantity measures the magnitude of parallel trends violations in the
 * pre-treatment period. Under the conditional extrapolation assumption
 * (Assumption 3), if S_pre <= M, then S_post <= S_pre.
 *
 * Special cases:
 *   p = 1:   Mean absolute violation (L1 norm / n)
 *   p = 2:   Root mean squared violation (L2 norm / sqrt(n))
 *   p = inf: Maximum absolute violation (L-infinity norm)
 *
 * @param nu_vec Violation vector (iterative nu_t or overall nubar_t)
 * @param p      Norm index p >= 1, or missing/1e10 for L-infinity
 *
 * @return Estimated severity S_pre_hat
 *
 * @note Default p = 2 (root mean square) if not specified
 * @note For p = infinity, uses max(|nu|) directly
 * @note Numerical stability: Uses scaling for large p to avoid underflow
 *
 * @see Mikhaeil & Harshaw (2025), Section 3.1
 */
real scalar _pretest_severity(real colvector nu_vec, | real scalar p)
{
    real scalar n, i, sum_abs_p
    
    // Default p = 2
    if (args() < 2 | missing(p)) {
        p = 2
    }
    
    n = rows(nu_vec)
    
    if (n == 0) {
        return(0)
    }
    
    // Check for missing values
    if (missing(nu_vec)) {
        return(.)
    }
    
    // p = infinity case (including p >= 1e10)
    if (p >= 1e10) {
        return(max(abs(nu_vec)))
    }
    
    // p = 1 case: MEAN of absolute values
    // S_pre = mean(|nu_hat_t|)
    if (p == 1) {
        return(mean(abs(nu_vec)))
    }
    
    // p = 2 case: Root MEAN square
    // S_pre = sqrt(mean(nu_hat_t^2))
    if (p == 2) {
        return(sqrt(mean(nu_vec:^2)))
    }
    
    // General p-norm
    if (p <= 0) {
        return(.)  // Invalid p
    }
    
    // General p-norm: S_pre = ((1/n) * sum(|nu_hat_t|^p))^(1/p)
    // 
    // Numerical stability note: For large p (e.g., p >= 100), direct computation
    // causes underflow when |nu| < 1 (e.g., 0.2^1000 ≈ 0).
    // 
    // Solution: Use LOG-SPACE computation to avoid underflow.
    // 
    // Mathematical derivation:
    //   S_pre = (mean(|nu|^p))^(1/p)
    //   log(S_pre) = (1/p) * log(mean(|nu|^p))
    //              = (1/p) * log(mean(exp(p * log(|nu|))))
    // 
    // Using log-sum-exp trick with max_log = max(p * log(|nu|)) = p * log(max(|nu|)):
    //   log(S_pre) = (1/p) * [max_log + log(mean(exp(p*log(|nu|) - max_log)))]
    //              = log(max|nu|) + (1/p) * log(mean(exp(p * log(|nu|/max|nu|))))
    // 
    // Since |nu|/max|nu| <= 1, log(|nu|/max|nu|) <= 0, and p*log(...) <= 0
    // So exp(p*log(...)) is in (0, 1], which won't overflow or underflow to 0.
    
    real colvector abs_nu, log_abs_nu, log_ratios, exp_terms
    real scalar max_abs, log_max_abs, mean_exp, log_S_pre
    
    abs_nu = abs(nu_vec)
    max_abs = max(abs_nu)
    
    // Handle zero vector case
    if (max_abs == 0 | max_abs < 1e-300) {
        return(0)
    }
    
    // Handle zeros in abs_nu (replace with very small positive to avoid log(0))
    abs_nu = abs_nu :+ (abs_nu :== 0) * 1e-300
    
    // For p >= 1e4, directly return L∞ (max norm) as the limit
    if (p >= 1e4) {
        return(max_abs)
    }
    
    // Log-space computation for numerical stability (works for ALL p >= 3)
    // This avoids both underflow (small^large_p) and overflow issues
    
    log_max_abs = log(max_abs)
    log_abs_nu = log(abs_nu)
    
    // log_ratios = log(|nu| / max|nu|) = log(|nu|) - log(max|nu|)
    // All values <= 0, with maximum = 0
    log_ratios = log_abs_nu :- log_max_abs
    
    // exp_terms = exp(p * log_ratios) = (|nu|/max|nu|)^p
    // All values in (0, 1] since log_ratios <= 0 and p > 0
    exp_terms = exp(p :* log_ratios)
    
    // mean_exp = mean((|nu|/max|nu|)^p)
    mean_exp = mean(exp_terms)
    
    // S_pre = max|nu| * mean_exp^(1/p)
    //       = exp(log_max_abs + (1/p) * log(mean_exp))
    if (mean_exp <= 0 | missing(mean_exp)) {
        // Fallback to L∞ if numerical issues occur
        return(max_abs)
    }
    
    log_S_pre = log_max_abs + (1/p) * log(mean_exp)
    
    return(exp(log_S_pre))
}


// ============================================================================
// PRETEST FUNCTION
// ============================================================================

/**
 * @function _pretest_test
 * @brief Execute the preliminary test for conditional extrapolation
 *
 * Implements the pre-test function from Section 4.2, Theorem 1:
 *
 *   phi = 1{S_pre > M}
 *
 * where M is the acceptable threshold for violation severity chosen by
 * the analyst based on substantive considerations (see Section 3.1).
 *
 * @param S_pre Estimated severity measure S_pre_hat
 * @param M     Acceptable threshold M > 0 (analyst-specified)
 *
 * @return Pre-test result phi in {0, 1, .}:
 *         phi = 0: PASS - S_pre <= M, extrapolation condition satisfied
 *         phi = 1: FAIL - S_pre > M, extrapolation not justified
 *         phi = .: Missing input, cannot evaluate
 *
 * @note Strict inequality: S_pre = M exactly yields phi = 0 (PASS)
 * @note Theorem 1 establishes asymptotic consistency of this test
 *
 * @see Mikhaeil & Harshaw (2025), Section 4.2, Theorem 1
 */
real scalar _pretest_test(real scalar S_pre, real scalar M)
{
    // Validate M
    if (M <= 0 | missing(M)) {
        _error("Threshold M must be positive")
        return(.)
    }
    
    if (missing(S_pre)) {
        return(.)
    }
    
    // phi = 1{S_pre > M}
    // Note: S_pre = M returns 0 (pass) due to strict inequality
    return(S_pre > M)
}


// ============================================================================
// SEVERITY GRADIENT (for Delta Method)
// ============================================================================

/**
 * @function _pretest_severity_gradient
 * @brief Compute gradient of severity function (for Delta Method)
 *
 * Computes the gradient vector of the severity function g(nu) = S_pre
 * with respect to the violation vector nu:
 *
 *   dg/dnu_t = |nu_t|^{p-1} * sign(nu_t) / [n * S_pre^{p-1}]
 *
 * Used for standard error estimation via the Delta Method.
 *
 * @param nu Violation vector
 * @param p  Norm index (default: 2)
 *
 * @return Gradient vector with same dimension as nu
 *
 * @note Returns missing for p = infinity (max not differentiable)
 * @note Returns zero vector if S_pre approx 0 (numerical stability)
 *
 * @see Mikhaeil & Harshaw (2025), Section 5.1
 */
real colvector _pretest_severity_gradient(real colvector nu, | real scalar p)
{
    real scalar n_pre, S_pre
    real colvector grad
    
    // Default p = 2
    if (args() < 2 | missing(p)) {
        p = 2
    }
    
    n_pre = rows(nu)
    
    if (n_pre == 0) {
        return(J(0, 1, 0))
    }
    
    // p = ∞ case: not differentiable
    if (p >= 1e10) {
        return(J(n_pre, 1, .))  // Return missing for p=∞
    }
    
    // Compute severity
    S_pre = _pretest_severity(nu, p)
    
    // Handle near-zero severity (numerical stability)
    if (S_pre < 1e-10) {
        return(J(n_pre, 1, 0))
    }
    
    // Gradient formula:
    // ∂g/∂ν_t = |ν_t|^{p-1} × sgn(ν_t) / [n_pre × S_pre^{p-1}]
    //
    // Note: n_pre = T_pre - 1 (number of violations)
    // The denominator n_pre comes from the mean() in severity definition
    
    grad = (abs(nu):^(p-1) :* sign(nu)) :/ (n_pre * S_pre^(p-1))
    
    return(grad)
}


// ============================================================================
// SEVERITY STANDARD ERROR (Delta Method)
// ============================================================================

/**
 * @function _pretest_severity_se
 * @brief Estimate standard error of severity via Delta Method
 *
 * Computes the standard error of the severity estimator using
 * the Delta Method:
 *
 *   SE(S_pre_hat) = sqrt(grad_g' * Sigma_nu * grad_g)
 *
 * where grad_g is the gradient of the severity function.
 *
 * @param nu       Violation vector
 * @param p        Norm index
 * @param Sigma_nu Covariance matrix of violations (n_pre x n_pre)
 *
 * @return Standard error SE(S_pre_hat), or:
 *         - Missing (.) if p = infinity (use bootstrap instead)
 *         - Zero if S_pre approximately 0
 *
 * @note Delta Method requires differentiability; not valid for p = inf
 *
 * @see Mikhaeil & Harshaw (2025), Section 5.1
 */
real scalar _pretest_severity_se(real colvector nu, real scalar p, 
                                  real matrix Sigma_nu)
{
    real scalar n_pre, S_pre, variance
    real colvector grad
    
    n_pre = rows(nu)
    
    // Input validation
    if (n_pre == 0) {
        return(0)
    }
    
    // Check covariance matrix dimensions
    if (rows(Sigma_nu) != n_pre | cols(Sigma_nu) != n_pre) {
        _error("Sigma_nu dimension mismatch: expected " + strofreal(n_pre) + 
               "x" + strofreal(n_pre))
        return(.)
    }
    
    // p = ∞ case: Delta Method not applicable
    if (p >= 1e10) {
        return(.)  // Return missing, user should use Bootstrap
    }
    
    // Compute severity
    S_pre = _pretest_severity(nu, p)
    
    // Handle near-zero severity
    if (S_pre < 1e-10) {
        return(0)
    }
    
    // Compute gradient
    grad = _pretest_severity_gradient(nu, p)
    
    // Check gradient computation succeeded
    if (missing(grad)) {
        return(.)
    }
    
    // Delta Method: Var(g(θ̂)) ≈ ∇g' Σ ∇g
    // SE = sqrt(∇g' × Σ × ∇g)
    variance = grad' * Sigma_nu * grad
    
    // Handle potential numerical issues
    if (variance < 0) {
        // Numerical precision issue, treat as zero
        return(0)
    }
    
    return(sqrt(variance))
}


// ============================================================================
// COVARIANCE MATRIX EXTRACTION
// ============================================================================

/**
 * @function _pretest_extract_sigma_nu
 * @brief Extract violation covariance submatrix from full Sigma
 *
 * Extracts the upper-left block Sigma_nu from the full covariance matrix:
 *
 *   Sigma = | Sigma_nu       Sigma_nu,delta |
 *           | Sigma_delta,nu Sigma_delta    |
 *
 * where Sigma_nu is the (T_pre-1) x (T_pre-1) covariance of violations.
 *
 * @param Sigma_hat Full covariance matrix from _pretest_covariance()
 * @param T_pre     Number of pre-treatment periods (t0 - 1)
 *
 * @return Sigma_nu of dimension (T_pre-1) x (T_pre-1)
 *
 * @note For overall mode, first transform Sigma using _pretest_transform_sigma()
 *
 * @see Mikhaeil & Harshaw (2025), Section 2.2
 */
real matrix _pretest_extract_sigma_nu(real matrix Sigma_hat, real scalar T_pre)
{
    real scalar n_nu
    
    // Number of violations = T_pre - 1
    n_nu = T_pre - 1
    
    if (n_nu < 1) {
        return(J(0, 0, .))
    }
    
    // Validate Sigma_hat dimensions
    if (rows(Sigma_hat) < n_nu | cols(Sigma_hat) < n_nu) {
        _error("Sigma_hat too small: need at least " + strofreal(n_nu) + 
               " rows/cols")
        return(J(0, 0, .))
    }
    
    // Extract upper-left block: Σ_nu = Σ̂[1..n_nu, 1..n_nu]
    return(Sigma_hat[1..n_nu, 1..n_nu])
}


// ============================================================================
// THETA VECTOR CONSTRUCTION
// ============================================================================

/**
 * @function _pretest_theta_vector
 * @brief Construct full parameter vector theta_hat
 *
 * Assembles the parameter vector for asymptotic analysis (Section 2.2):
 *
 *   theta = (nu_2, ..., nu_{t0-1}, delta_{t0}, ..., delta_T)'
 *
 * In overall mode, nu is replaced by nubar (cumulative violations).
 *
 * @param Y            Outcome vector (n x 1)
 * @param D            Treatment indicator vector (n x 1)
 * @param time         Time period vector (n x 1)
 * @param t0           Treatment time (first post-treatment period)
 * @param T            Total number of time periods
 * @param overall_mode 0 for iterative (default), 1 for overall
 *
 * @return Parameter vector theta_hat of dimension (T-1) x 1
 *
 * @see Mikhaeil & Harshaw (2025), Section 2.2, Equation (205)
 */
real colvector _pretest_theta_vector(real colvector Y, real colvector D,
                                     real colvector time, real scalar t0,
                                     real scalar T, | real scalar overall_mode)
{
    real colvector violations, did_vec, theta
    real scalar T_post
    
    // Default to iterative mode
    if (args() < 6 | missing(overall_mode)) {
        overall_mode = 0
    }
    
    T_post = T - t0
    
    // Get violations
    if (overall_mode == 1) {
        violations = _pretest_nubar_vector(Y, D, time, t0)
    } else {
        violations = _pretest_nu_vector(Y, D, time, t0)
    }
    
    // Get DID estimators (delta_{t0}, delta_{t0+1}, ..., delta_T)
    // Paper Section 2.2 Eq. (205): theta = (..., delta_{t0}, ..., delta_T)
    // Includes delta_{t0} even though it equals 0 by definition
    did_vec = _pretest_did_vector(Y, D, time, t0, T)
    
    // Concatenate: θ̂ = (violations', did_vec')'
    theta = violations \ did_vec
    
    return(theta)
}


// ============================================================================
// CONVENIENCE FUNCTIONS
// ============================================================================

/**
 * @function _pretest_T_pre
 * @brief Calculate number of pre-treatment violation parameters
 *
 * @param t0           Treatment time
 * @param overall_mode 0 for iterative (default), 1 for overall
 *
 * @return Number of violation parameters:
 *         - Iterative: T_pre - 1 = t0 - 2
 *         - Overall: t0 - 1
 */
real scalar _pretest_T_pre(real scalar t0, | real scalar overall_mode)
{
    if (args() < 2 | missing(overall_mode)) {
        overall_mode = 0
    }
    
    if (overall_mode == 1) {
        return(t0 - 1)
    } else {
        return(t0 - 2)
    }
}


/**
 * @function _pretest_T_post
 * @brief Calculate number of post-treatment periods
 *
 * Computes T_post following the paper's definition (Section 2.1):
 *   T_post = T - T_pre = T - (t0 - 1) = T - t0 + 1
 *
 * This includes the treatment period t0 itself as a post-treatment period.
 *
 * @param t0 Treatment time (first post-treatment period)
 * @param T  Total number of time periods
 *
 * @return Number of post-treatment periods: T_post = T - t0 + 1
 *
 * @see Mikhaeil & Harshaw (2025), Section 2.1
 */
real scalar _pretest_T_post(real scalar t0, real scalar T)
{
    return(T - t0 + 1)
}


end
