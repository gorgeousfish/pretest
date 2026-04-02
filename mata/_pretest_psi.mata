*! _pretest_psi.mata
*! Psi Function and Monte Carlo Critical Value Computation
*!
*! Description:
*!   Implements the psi function and Monte Carlo critical value computation
*!   for conditionally valid confidence intervals under the conditional
*!   extrapolation assumption in difference-in-differences designs.
*!
*! Psi Function (Theorem 2):
*!   Iterative mode: psi(x) = |mean(x_post)| + kappa * S(x_pre, p)
*!   Overall mode:   psi^Delta(x) = |mean(x_post)| + S(x_pre, p)
*!
*! where S(v, p) = (mean(|v|^p))^{1/p} is the severity measure.
*!
*! Critical Value (Appendix D.5):
*!   f(alpha, Sigma) = inf{c : P(psi(Z) >= c) <= alpha}, Z ~ N(0, Sigma)
*!
*! Monte Carlo Algorithm:
*!   1. Draw Z^{(s)} ~ N(0, Sigma) via Cholesky decomposition
*!   2. Compute psi^{(s)} = psi(Z^{(s)}) for s = 1, ..., S
*!   3. Return the ceiling((1-alpha)*S)-th order statistic
*!
*! Reference:
*!   Mikhaeil, J. M. and C. Harshaw. 2025. In Defense of the Pre-Test:
*!   Valid Inference when Testing Violations of Parallel Trends for
*!   Difference-in-Differences. arXiv preprint arXiv:2510.26470.

version 17.0

mata:
mata set matastrict on

// ============================================================================
// PSI FUNCTION (ITERATIVE MODE)
// ============================================================================

/**
 * _pretest_psi()
 *
 * Computes the psi function for critical value calculation in iterative mode.
 *
 * Mathematical Definition (Theorem 2):
 *
 *   psi(x) = |mean(x_post)| + kappa * S(x_pre, p)
 *
 * where S(v, p) = (mean(|v|^p))^{1/p} is the severity measure (Section 3.1).
 *
 * The psi function characterizes the distribution of the worst-case error
 * when the estimation error vector x ~ N(0, Sigma) is drawn from the
 * asymptotic distribution of sqrt(n) * (theta_hat - theta).
 *
 * Parameters:
 *   x      - Estimation error vector, dimension (T-1) x 1
 *            x[1..(T_pre-1)]: pre-treatment iterative violation errors
 *            x[T_pre..(T-1)]: post-treatment DID estimator errors
 *   T_pre  - Number of pre-treatment periods (= t0 - 1)
 *   T_post - Number of post-treatment periods (= T - t0 + 1)
 *   p      - L_p norm index for severity measure (p >= 1; p >= 1e10 for L_inf)
 *   kappa  - Bias amplification constant from Proposition 1
 *
 * Returns:
 *   Scalar psi value = |mean(x_post)| + kappa * severity(x_pre, p)
 *
 * Notes:
 *   - Dimension requirement: rows(x) = (T_pre - 1) + T_post = T - 1
 *   - Edge case: when T_pre < 2, returns only the post-treatment component
 *   - Uses log-space computation for numerical stability when p >= 3
 *
 * Reference:
 *   Mikhaeil & Harshaw (2025), Section 5.1, Theorem 2
 */
real scalar _pretest_psi(real colvector x, real scalar T_pre,
                         real scalar T_post, real scalar p,
                         real scalar kappa)
{
    real scalar term1, term2, n_pre
    real colvector x_pre, x_post
    
    // The input vector x has dimension (T-1), partitioned as:
    //   x[1..(T_pre-1)]: pre-treatment iterative violation errors
    //   x[T_pre..(T-1)]: post-treatment DID estimator errors
    n_pre = T_pre - 1
    
    if (rows(x) != n_pre + T_post) {
        _error("_pretest_psi: dimension mismatch - expected " + 
               strofreal(n_pre + T_post) + " but got " + strofreal(rows(x)))
        return(.)
    }
    
    // Edge case: when T_pre = 1, there are no pre-treatment violations
    // In this case, the severity term is undefined and we return only
    // the post-treatment component |mean(x_post)|
    if (n_pre < 1) {
        x_post = x[1..T_post]
        return(abs(mean(x_post)))
    }
    
    // Partition the error vector following the paper's notation (Section 2.2)
    x_pre = x[1..n_pre]
    x_post = x[(n_pre + 1)..rows(x)]
    
    // Term 1: Absolute mean of post-treatment DID errors
    // Corresponds to the estimation error in delta_bar_hat
    term1 = abs(mean(x_post))
    
    // Term 2: Bias bound = kappa * severity(x_pre, p)
    // The severity measure S(v,p) = (mean(|v|^p))^{1/p} quantifies
    // the magnitude of pre-treatment parallel trend violations
    if (p >= 1e10) {
        // L_infinity norm: S(v, inf) = max(|v|)
        term2 = kappa * max(abs(x_pre))
    } else {
        // General L_p norm with numerical stability considerations
        real colvector abs_x, log_abs_x, log_ratios, exp_terms
        real scalar max_abs, log_max_abs, mean_exp, severity
        
        abs_x = abs(x_pre)
        max_abs = max(abs_x)
        
        if (max_abs == 0 | max_abs < 1e-300) {
            // All violations are zero; severity is zero
            term2 = 0
        } else if (p >= 3) {
            // For large p, use log-space computation to prevent overflow
            // of |x|^p when |x| > 1 and underflow when |x| < 1
            abs_x = abs_x :+ (abs_x :== 0) * 1e-300  // Handle exact zeros
            
            log_max_abs = log(max_abs)
            log_abs_x = log(abs_x)
            log_ratios = log_abs_x :- log_max_abs    // Normalized: all <= 0
            exp_terms = exp(p :* log_ratios)         // Scaled: all in (0, 1]
            mean_exp = mean(exp_terms)
            
            if (mean_exp <= 0 | missing(mean_exp)) {
                severity = max_abs  // Fallback for numerical edge cases
            } else {
                // severity = max_abs * (mean_exp)^{1/p}
                severity = exp(log_max_abs + (1/p) * log(mean_exp))
            }
            term2 = kappa * severity
        } else {
            // Standard computation for moderate p
            // S(v,p) = (mean(|v|^p))^{1/p}
            term2 = kappa * (mean(abs_x:^p))^(1/p)
        }
    }
    
    return(term1 + term2)
}


// ============================================================================
// PSI DELTA FUNCTION (OVERALL MODE)
// ============================================================================

/**
 * _pretest_psi_delta()
 *
 * Computes the psi^Delta function for critical value calculation in overall mode.
 *
 * Mathematical Definition (Appendix C):
 *
 *   psi^Delta(x) = |mean(x_post)| + S(x_pre, p)
 *
 * where S(v, p) = (mean(|v|^p))^{1/p} is the severity measure.
 *
 * Key Distinction from Iterative Mode:
 *   In overall mode, cumulative violations (nubar_t) directly bound the
 *   post-treatment bias without the kappa amplification factor:
 *
 *     |tau_bar - delta_bar| <= S_pre^Delta     (Appendix C, Proposition)
 *
 *   This contrasts with iterative mode where:
 *
 *     |tau_bar - delta_bar| <= kappa * S_pre   (Proposition 1)
 *
 * Parameters:
 *   x      - Estimation error vector, dimension (T-1) x 1
 *            x[1..(T_pre-1)]: pre-treatment cumulative violation errors
 *            x[T_pre..(T-1)]: post-treatment DID estimator errors
 *   T_pre  - Number of pre-treatment periods (= t0 - 1)
 *   T_post - Number of post-treatment periods (= T - t0 + 1)
 *   p      - L_p norm index for severity measure (p >= 1; p >= 1e10 for L_inf)
 *
 * Returns:
 *   Scalar psi^Delta value = |mean(x_post)| + severity(x_pre, p)
 *
 * Notes:
 *   - No kappa multiplier: this is the key distinction from iterative mode
 *   - Overall violations represent cumulative deviation from parallel trends
 *   - The covariance matrix Sigma^Delta differs from the iterative Sigma
 *
 * Reference:
 *   Mikhaeil & Harshaw (2025), Appendix C
 */
real scalar _pretest_psi_delta(real colvector x, real scalar T_pre,
                               real scalar T_post, real scalar p)
{
    real scalar term1, term2, n_pre
    real colvector x_pre, x_post
    
    // Validate input dimension
    n_pre = T_pre - 1
    
    if (rows(x) != n_pre + T_post) {
        _error("_pretest_psi_delta: dimension mismatch - expected " + 
               strofreal(n_pre + T_post) + " but got " + strofreal(rows(x)))
        return(.)
    }
    
    // Edge case: no pre-treatment violations available
    if (n_pre < 1) {
        x_post = x[1..T_post]
        return(abs(mean(x_post)))
    }
    
    // Partition the error vector
    x_pre = x[1..n_pre]
    x_post = x[(n_pre + 1)..rows(x)]
    
    // Term 1: Absolute mean of post-treatment DID errors
    term1 = abs(mean(x_post))
    
    // Term 2: Severity measure S(x_pre, p) WITHOUT kappa multiplier
    // This is the key distinction from iterative mode (see Appendix C)
    if (p >= 1e10) {
        // L_infinity norm: S(v, inf) = max(|v|)
        term2 = max(abs(x_pre))
    } else {
        // General L_p norm with numerical stability considerations
        real colvector abs_x2, log_abs_x2, log_ratios2, exp_terms2
        real scalar max_abs2, log_max_abs2, mean_exp2
        
        abs_x2 = abs(x_pre)
        max_abs2 = max(abs_x2)
        
        if (max_abs2 == 0 | max_abs2 < 1e-300) {
            // All violations are zero
            term2 = 0
        } else if (p >= 3) {
            // Log-space computation for numerical stability
            abs_x2 = abs_x2 :+ (abs_x2 :== 0) * 1e-300
            
            log_max_abs2 = log(max_abs2)
            log_abs_x2 = log(abs_x2)
            log_ratios2 = log_abs_x2 :- log_max_abs2
            exp_terms2 = exp(p :* log_ratios2)
            mean_exp2 = mean(exp_terms2)
            
            if (mean_exp2 <= 0 | missing(mean_exp2)) {
                term2 = max_abs2
            } else {
                term2 = exp(log_max_abs2 + (1/p) * log(mean_exp2))
            }
        } else {
            // Standard computation for moderate p
            term2 = (mean(abs_x2:^p))^(1/p)
        }
    }
    
    return(term1 + term2)
}


// ============================================================================
// MONTE CARLO CRITICAL VALUE
// ============================================================================

/**
 * _pretest_critical_value()
 *
 * Computes the Monte Carlo critical value f(alpha, Sigma) for confidence
 * interval construction under the conditional extrapolation assumption.
 *
 * Mathematical Definition (Appendix F):
 *
 *   f(alpha, Sigma) = inf{c : P(psi(Z) >= c) <= alpha}
 *
 * where Z ~ N(0, Sigma) represents the asymptotic distribution of the
 * scaled estimation error sqrt(n) * (theta_hat - theta).
 *
 * The critical value determines the width of the confidence interval:
 *
 *   CI = delta_bar_hat +/- {kappa * S_pre_hat + f(alpha, Sigma_hat) / sqrt(n)}
 *
 * Algorithm (Monte Carlo Simulation):
 *   1. Regularize: Sigma_reg = Sigma + epsilon*I for numerical stability
 *   2. Cholesky decomposition: L such that Sigma_reg = LL'
 *   3. For s = 1, ..., S:
 *      a. Draw W^{(s)} ~ N(0, I_dim) (standard normal)
 *      b. Transform: Z^{(s)} = L * W^{(s)} ~ N(0, Sigma_reg)
 *      c. Evaluate: psi^{(s)} = psi(Z^{(s)})
 *   4. Return ceil((1-alpha)*S)-th order statistic of {psi^{(s)}}
 *
 * Parameters:
 *   Sigma_hat    - Estimated asymptotic covariance matrix, (T-1) x (T-1)
 *   alpha        - Significance level in (0, 1), typically 0.05
 *   S            - Number of Monte Carlo draws (recommended: >= 5000)
 *   T_pre        - Number of pre-treatment periods (= t0 - 1)
 *   T_post       - Number of post-treatment periods (= T - t0 + 1)
 *   p            - L_p norm index for severity (p >= 1; p >= 1e10 for L_inf)
 *   kappa        - Bias bound constant from _pretest_kappa()
 *   overall_mode - (optional) 0 = iterative mode (default), 1 = overall mode
 *   seed         - (optional) Random seed for reproducibility (default: 12345)
 *
 * Returns:
 *   Scalar critical value f(alpha, Sigma_hat) for CI construction
 *
 * Notes:
 *   - Cholesky decomposition: Z = LW with W ~ N(0,I) yields Z ~ N(0,Sigma)
 *   - Regularization epsilon = 1e-10 handles near-singular matrices
 *   - Ceiling index ensures conservative critical value estimation
 *   - The quantile function F^{-1}(1-alpha) is continuous (Lemma in Appendix D)
 *
 * Reference:
 *   Mikhaeil & Harshaw (2025), Appendix F (Section D.5)
 */
real scalar _pretest_critical_value(real matrix Sigma_hat,
                                    real scalar alpha,
                                    real scalar S,
                                    real scalar T_pre,
                                    real scalar T_post,
                                    real scalar p,
                                    real scalar kappa,
                                    | real scalar overall_mode,
                                    real scalar seed)
{
    real scalar i, dim, f_alpha, quantile_idx
    real matrix L, Sigma_reg
    real colvector W, Z, psi_values
    
    // Handle optional arguments
    if (args() < 8 | missing(overall_mode)) {
        overall_mode = 0
    }
    if (args() < 9 | missing(seed)) {
        seed = 12345
    }
    
    // Input validation
    dim = rows(Sigma_hat)
    if (dim != cols(Sigma_hat)) {
        _error("_pretest_critical_value: Sigma_hat must be square")
        return(.)
    }
    
    if (dim != (T_pre - 1) + T_post) {
        _error("_pretest_critical_value: Sigma_hat dimension mismatch - expected " + 
               strofreal((T_pre - 1) + T_post) + " but got " + strofreal(dim))
        return(.)
    }
    
    // Initialize random number generator for reproducibility
    rseed(seed)
    
    // Step 1: Regularize covariance matrix to ensure positive definiteness
    // This handles near-singular matrices that may arise from small samples
    Sigma_reg = Sigma_hat + 1e-10 * I(dim)
    
    // Step 2: Cholesky decomposition for sampling from N(0, Sigma)
    // Mata's cholesky() returns lower triangular L where Sigma = L * L'
    L = cholesky(Sigma_reg)
    
    // Step 3-4: Monte Carlo simulation
    psi_values = J(S, 1, .)
    
    for (i = 1; i <= S; i++) {
        // Draw from standard multivariate normal: W ~ N(0, I)
        W = rnormal(dim, 1, 0, 1)
        
        // Transform to target distribution: Z = L * W ~ N(0, Sigma_reg)
        Z = L * W
        
        // Evaluate psi function at the simulated point
        if (overall_mode == 1) {
            psi_values[i] = _pretest_psi_delta(Z, T_pre, T_post, p)
        } else {
            psi_values[i] = _pretest_psi(Z, T_pre, T_post, p, kappa)
        }
    }
    
    // Step 5: Extract (1-alpha) quantile from empirical distribution
    _sort(psi_values, 1)  // Sort in ascending order
    
    // Use ceiling for conservative critical value
    // This ensures nominal coverage is achieved asymptotically
    quantile_idx = ceil(S * (1 - alpha))
    
    // Bound index to valid range
    if (quantile_idx < 1) quantile_idx = 1
    if (quantile_idx > S) quantile_idx = S
    
    f_alpha = psi_values[quantile_idx]
    
    return(f_alpha)
}


// ============================================================================
// CONVENIENCE FUNCTIONS
// ============================================================================

/**
 * _pretest_critical_value_simple()
 *
 * Simplified interface for computing the critical value with default parameters.
 *
 * Uses the following defaults:
 *   - alpha = 0.05 (5% significance level)
 *   - S = 5000 (Monte Carlo draws)
 *   - p = 2 (L_2 norm for severity)
 *   - seed = 12345 (random seed)
 *   - mode = iterative
 *
 * Parameters:
 *   Sigma_hat - Asymptotic covariance matrix, (T-1) x (T-1)
 *   T_pre     - Number of pre-treatment periods
 *   T_post    - Number of post-treatment periods
 *   kappa     - Bias amplification constant
 *
 * Returns:
 *   Critical value f(0.05, Sigma_hat) at 5% significance level
 */
real scalar _pretest_critical_value_simple(real matrix Sigma_hat,
                                           real scalar T_pre,
                                           real scalar T_post,
                                           real scalar kappa)
{
    return(_pretest_critical_value(Sigma_hat, 0.05, 5000, T_pre, T_post, 2, kappa, 0, 12345))
}


/**
 * _pretest_verify_cholesky()
 *
 * Diagnostic function to verify Cholesky decomposition correctness.
 *
 * Checks that L * L' = Sigma within the specified tolerance, where L is
 * the lower triangular Cholesky factor. Used for unit testing and
 * validating numerical stability of the covariance matrix operations.
 *
 * Parameters:
 *   Sigma - Covariance matrix to decompose
 *   tol   - (optional) Tolerance for comparison (default: 1e-10)
 *
 * Returns:
 *   1 if max|L*L' - Sigma| < tol, 0 otherwise
 */
real scalar _pretest_verify_cholesky(real matrix Sigma, | real scalar tol)
{
    real matrix L, Sigma_reconstructed
    real scalar max_error
    
    if (args() < 2 | missing(tol)) {
        tol = 1e-10
    }
    
    // Compute Cholesky factor with regularization
    L = cholesky(Sigma + 1e-10 * I(rows(Sigma)))
    
    // Reconstruct the original matrix
    Sigma_reconstructed = L * L'
    
    // Verify reconstruction accuracy
    max_error = max(abs(Sigma_reconstructed - Sigma))
    
    return(max_error < tol)
}


end
