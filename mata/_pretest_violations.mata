*! _pretest_violations.mata
*! Parallel Trends Violations Estimation Module
*!
*! Purpose:
*!   Implements estimators for parallel trends violations within the conditional
*!   extrapolation pre-test framework. Provides both iterative (nu_t) and overall
*!   (nubar_t) violation measures, severity computation, and pre-test evaluation.
*!
*! Theoretical Framework:
*!
*!   Iterative Violation (Section 2.1):
*!     The period-to-period deviation from parallel trends is defined as:
*!       nu_t = E[Y^{(0)}_t - Y^{(0)}_{t-1} | D=1] - E[Y^{(0)}_t - Y^{(0)}_{t-1} | D=0]
*!     This measures the differential change in untreated potential outcomes
*!     between treatment and control groups at each time period.
*!
*!   Overall Violation (Section 2.1):
*!     The cumulative deviation relative to the baseline period:
*!       nubar_t = sum_{s=2}^t nu_s
*!     Equivalence relations: nubar_t = sum nu_s and nu_t = nubar_t - nubar_{t-1}
*!
*!   Severity Measure (Section 3.1):
*!     The severity of pre-treatment parallel trends violations is:
*!       S_pre = ( (1/(T_pre-1)) * sum_{t=2}^{t0-1} |nu_t|^p )^{1/p}
*!     This normalized Lp-norm quantifies violation magnitude. Under the
*!     conditional extrapolation assumption (Assumption 3):
*!       If S_pre <= M, then S_post <= S_pre
*!
*!   Pre-Test (Section 4.2, Theorem 1):
*!     phi = 1{S_pre_hat > M}
*!     phi = 0: PASS - extrapolation condition satisfied, proceed with inference
*!     phi = 1: FAIL - violations exceed threshold, extrapolation not justified
*!     Asymptotic consistency holds for separation s_n = omega(n^{-1/2}).
*!
*! Reference:
*!   Mikhaeil, J. M. and C. Harshaw. 2025. "In Defense of the Pre-Test: Valid
*!   Inference when Testing Violations of Parallel Trends for Difference-in-
*!   Differences." arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470

version 17.0

mata:
mata set matastrict on

// ============================================================================
// ITERATIVE VIOLATIONS
// ============================================================================

/**
 * @brief Compute single-period iterative violation of parallel trends
 *
 * Implements the iterative violation estimator (Section 2.1):
 *
 *   nu_hat_t = [Y_bar_{1,t} - Y_bar_{1,t-1}] - [Y_bar_{0,t} - Y_bar_{0,t-1}]
 *
 * where Y_bar_{d,t} denotes the sample mean of outcomes for group D=d at time t.
 * This difference-in-differences of adjacent periods captures the deviation
 * from parallel trends at each time point.
 *
 * @param Y    Outcome vector (n x 1)
 * @param D    Treatment indicator (n x 1), binary {0, 1}
 * @param time Time period vector (n x 1)
 * @param t    Target period for violation calculation, t in {2, ..., t0-1}
 * @param t0   Treatment onset time (first post-treatment period)
 *
 * @return Scalar nu_hat_t, or missing (.) if data insufficient
 *
 * @see Mikhaeil & Harshaw (2025), Section 2.1
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
 * @brief Compute full vector of iterative violations
 *
 * Constructs the pre-treatment iterative violation vector:
 *
 *   nu = (nu_hat_2, nu_hat_3, ..., nu_hat_{t0-1})'
 *
 * This vector constitutes the first (t0-2) elements of the parameter vector
 * theta defined in Section 2.2 for asymptotic analysis.
 *
 * @param Y    Outcome vector (n x 1)
 * @param D    Treatment indicator (n x 1)
 * @param time Time period vector (n x 1)
 * @param t0   Treatment onset time
 *
 * @return Column vector of dimension (t0-2) x 1
 *
 * @see Mikhaeil & Harshaw (2025), Section 2.2
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
 * @brief Compute single-period overall violation via cumulative summation
 *
 * Implements the overall violation (Section 2.1):
 *
 *   nubar_t = sum_{s=2}^t nu_s
 *
 * Overall violations represent cumulative deviation from parallel trends
 * relative to the baseline period t=1.
 *
 * @param nu_vec Iterative violation vector (nu_hat_2, ..., nu_hat_{t0-1})'
 * @param k      Index (1-based): nu_vec[k] corresponds to nu_hat_{k+1}
 *
 * @return Overall violation nubar_hat_{k+1} = sum(nu_vec[1..k])
 *
 * @see Mikhaeil & Harshaw (2025), Section 2.1
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
 * @brief Compute full vector of overall violations
 *
 * Constructs the pre-treatment overall violation vector:
 *
 *   nubar = (nubar_hat_2, nubar_hat_3, ..., nubar_hat_{t0-1})'
 *
 * where nubar_t = sum_{s=2}^t nu_s represents cumulative deviation from
 * parallel trends. Computed via cumulative summation of iterative violations.
 *
 * @param Y    Outcome vector (n x 1)
 * @param D    Treatment indicator (n x 1)
 * @param time Time period vector (n x 1)
 * @param t0   Treatment onset time
 *
 * @return Column vector of dimension (t0-2) x 1
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
 * @brief Convert iterative violations to overall violations
 *
 * Applies cumulative summation per the equivalence relation (Section 2.1):
 *
 *   nubar_t = sum_{s=2}^t nu_s
 *
 * @param nu Iterative violation vector (nu_hat_2, ..., nu_hat_{t0-1})'
 *
 * @return Overall violations (nubar_hat_2, ..., nubar_hat_{t0-1})'
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
 * @brief Convert overall violations to iterative violations
 *
 * Applies first-differencing per the equivalence relation (Section 2.1):
 *
 *   nu_t = nubar_t - nubar_{t-1}
 *
 * This is the inverse operation of _pretest_nu_to_nubar().
 *
 * @param nubar Overall violation vector (nubar_1, nubar_2, ..., nubar_{t0-1})'
 *              of dimension (t0-1) x 1, where nubar_1 = 0 by convention
 *
 * @return Iterative violation vector (nu_2, nu_3, ..., nu_{t0-1})'
 *         of dimension (t0-2) x 1
 *
 * @note nubar_1 = 0 serves as the baseline for differencing
 *
 * @see Mikhaeil & Harshaw (2025), Section 2.1, Appendix D
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
 * @brief Calculate the severity measure S_pre for pre-treatment violations
 *
 * Implements the severity measure (Section 3.1):
 *
 *   S_pre = ( (1/(T_pre-1)) * sum_{t=2}^{t0-1} |nu_t|^p )^{1/p}
 *
 * This normalized Lp-norm quantifies violation magnitude and plays a central
 * role in the conditional extrapolation assumption (Assumption 3):
 *   If S_pre <= M, then S_post <= S_pre
 *
 * Interpretation of p (Section 3.1):
 *   - p = 1:   Mean absolute violation, robust to outliers
 *   - p = 2:   Root mean squared violation (default)
 *   - p = inf: Maximum absolute violation, worst-case bound
 *
 * @param nu_vec Violation vector (iterative or overall)
 * @param p      Norm index p >= 1; use p >= 1e10 for L-infinity
 *
 * @return Severity S_pre_hat, 0 if empty, missing (.) if invalid
 *
 * @note Uses log-space computation for p >= 3 to ensure numerical stability
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
    // Numerical stability: For large p (e.g., p >= 100), direct computation
    // causes underflow when |nu| < 1 (e.g., 0.2^1000 -> 0).
    //
    // Solution: Log-space computation using the identity:
    //   S_pre = max|nu| * ( mean((|nu|/max|nu|)^p) )^{1/p}
    //
    // Since |nu|/max|nu| <= 1, the ratio^p stays in (0,1], avoiding underflow.
    
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
 * @brief Execute the preliminary test for conditional extrapolation
 *
 * Implements the pre-test indicator function (Section 4.2, Theorem 1):
 *
 *   phi = 1{S_pre_hat > M}
 *
 * where M is the analyst-specified acceptable threshold for violation severity.
 *
 * @param S_pre Estimated severity S_pre_hat
 * @param M     Acceptable threshold M > 0
 *
 * @return phi in {0, 1, .}:
 *         0: PASS - S_pre <= M, extrapolation justified
 *         1: FAIL - S_pre > M, extrapolation not justified
 *         .: Invalid input
 *
 * @note Boundary case: S_pre = M exactly yields phi = 0 (PASS)
 * @note Theorem 1: Asymptotically consistent for s_n = omega(n^{-1/2})
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
 * @brief Compute gradient of severity function for Delta Method inference
 *
 * Computes the gradient of g(nu) = S_pre with respect to nu:
 *
 *   dg/dnu_t = |nu_t|^{p-1} * sign(nu_t) / [n * S_pre^{p-1}]
 *
 * Required for standard error estimation via the Delta Method.
 *
 * @param nu Violation vector
 * @param p  Norm index (default: 2)
 *
 * @return Gradient vector (same dimension as nu)
 *
 * @note Returns missing for p = infinity (non-differentiable at maximum)
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
 * @brief Estimate standard error of severity via Delta Method
 *
 * Applies the Delta Method for asymptotic variance:
 *
 *   SE(S_pre_hat) = sqrt(grad_g' * Sigma_nu * grad_g)
 *
 * @param nu       Violation vector
 * @param p        Norm index
 * @param Sigma_nu Covariance matrix of violations
 *
 * @return Standard error SE(S_pre_hat), or missing for p = infinity
 *
 * @note For p = infinity, use bootstrap methods instead
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
 * @brief Extract violation covariance submatrix from full Sigma
 *
 * Extracts the upper-left block from the partitioned covariance:
 *
 *   Sigma = | Sigma_nu       Sigma_nu,delta |
 *           | Sigma_delta,nu Sigma_delta    |
 *
 * @param Sigma_hat Full covariance matrix
 * @param T_pre     Number of pre-treatment periods (t0 - 1)
 *
 * @return Sigma_nu of dimension (T_pre-1) x (T_pre-1)
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
 * @brief Construct full parameter vector theta_hat
 *
 * Assembles the parameter vector for asymptotic analysis (Section 2.2):
 *
 *   theta = (nu_2, ..., nu_{t0-1}, delta_{t0}, ..., delta_T)'
 *
 * @param Y            Outcome vector (n x 1)
 * @param D            Treatment indicator (n x 1)
 * @param time         Time period vector (n x 1)
 * @param t0           Treatment onset time
 * @param T            Total number of periods
 * @param overall_mode 0 = iterative (default), 1 = overall (nubar)
 *
 * @return theta_hat of dimension (T-1) x 1
 *
 * @see Mikhaeil & Harshaw (2025), Section 2.2
 */
real colvector _pretest_theta_vector(real colvector Y, real colvector D,
                                     real colvector time, real scalar t0,
                                     real scalar T, | real scalar overall_mode)
{
    real colvector violations, did_vec, theta
    
    // Default to iterative mode
    if (args() < 6 | missing(overall_mode)) {
        overall_mode = 0
    }
    
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
 * @brief Calculate number of pre-treatment violation parameters
 *
 * @param t0           Treatment onset time
 * @param overall_mode 0 = iterative (default), 1 = overall
 *
 * @return Number of violation parameters:
 *         - Iterative mode: t0 - 2 (for nu_2, ..., nu_{t0-1})
 *         - Overall mode: t0 - 1 (for nubar_1, ..., nubar_{t0-1})
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
 * @brief Calculate number of post-treatment periods
 *
 * Per Section 2.1: T_post = T - T_pre = T - (t0 - 1) = T - t0 + 1
 *
 * The treatment period t0 is included as a post-treatment period.
 *
 * @param t0 Treatment onset time (first post-treatment period)
 * @param T  Total number of time periods
 *
 * @return T_post = T - t0 + 1
 *
 * @see Mikhaeil & Harshaw (2025), Section 2.1
 */
real scalar _pretest_T_post(real scalar t0, real scalar T)
{
    return(T - t0 + 1)
}


end
