*! _pretest_psi.mata v0.1.0
*! Psi Function and Monte Carlo Critical Value Computation
*!
*! Implements the psi function and Monte Carlo simulation for computing
*! critical values in the conditional extrapolation pre-test framework.
*!
*! Reference:
*!   Mikhaeil & Harshaw (2025), "In Defense of the Pre-Test"
*!   arXiv:2510.26470, Section 5.1, Appendix D

version 17.0

mata:

// ============================================================================
// PSI FUNCTION (ITERATIVE MODE)
// ============================================================================

/**
 * @function _pretest_psi
 * @brief Calculate the psi function for iterative mode critical value computation
 *
 * Implements the psi function from Section 5.1, Theorem 2:
 *
 *   psi(x) = |mean(x_post)| + kappa * S(x_pre, p)
 *
 * where S(v, p) = (mean(|v|^p))^{1/p} is the severity measure.
 *
 * This function measures the "worst-case" deviation when the estimation error
 * vector x is drawn from the asymptotic distribution N(0, Sigma).
 *
 * @param x      Estimation error vector of dimension (T-1) x 1:
 *               - x[1..(T_pre-1)]: pre-treatment violation errors (nu_hat - nu)
 *               - x[T_pre..(T-1)]: post-treatment DID errors (delta_hat - delta)
 * @param T_pre  Number of pre-treatment periods (t0 - 1)
 * @param T_post Number of post-treatment periods (T - t0 + 1)
 * @param p      Norm index for severity (p >= 1, or >= 1e10 for infinity)
 * @param kappa  Bias bound constant from _pretest_kappa()
 *
 * @return Psi value: |mean(x_post)| + kappa * severity(x_pre, p)
 *
 * @note Dimension requirement: rows(x) = (T_pre - 1) + T_post = T - 1
 * @note When T_pre < 2, returns only the post-treatment component
 *
 * @see Mikhaeil & Harshaw (2025), Section 5.1, Theorem 2
 */
real scalar _pretest_psi(real colvector x, real scalar T_pre,
                         real scalar T_post, real scalar p,
                         real scalar kappa)
{
    real scalar term1, term2, n_pre
    real colvector x_pre, x_post
    
    // Dimension of x should be T-1 = (T_pre-1) + T_post
    // Validate
    n_pre = T_pre - 1  // Number of pre-treatment elements in x
    
    if (rows(x) != n_pre + T_post) {
        _error("_pretest_psi: dimension mismatch - expected " + 
               strofreal(n_pre + T_post) + " but got " + strofreal(rows(x)))
        return(.)
    }
    
    // Edge case: no pre-treatment violations (n_pre = 0)
    if (n_pre < 1) {
        // Only post-treatment component
        x_post = x[1..T_post]
        return(abs(mean(x_post)))
    }
    
    // Extract x_pre and x_post (1-indexed, consistent with paper)
    // x_pre = x[1..(T_pre-1)] - pre-treatment violation errors
    // x_post = x[T_pre..rows(x)] - post-treatment DID errors
    x_pre = x[1..n_pre]
    x_post = x[(n_pre + 1)..rows(x)]
    
    // Term 1: |mean(x_post)|
    term1 = abs(mean(x_post))
    
    // Term 2: κ × severity(x_pre, p)
    // severity = (mean(|x_pre|^p))^(1/p)
    // Note: For large p (approaching infinity), use max norm to avoid numerical underflow
    if (p >= 1e10) {
        // p = ∞: use max(|x_pre|)
        term2 = kappa * max(abs(x_pre))
    } else {
        // General p-norm with numerical stability
        real colvector abs_x
        real scalar max_abs, severity
        
        abs_x = abs(x_pre)
        max_abs = max(abs_x)
        
        if (max_abs == 0 | max_abs < 1e-300) {
            term2 = 0
        } else if (p >= 100) {
            // For large p, use scaling to avoid underflow
            real colvector scaled
            real scalar mean_scaled_p
            scaled = abs_x :/ max_abs
            mean_scaled_p = mean(scaled:^p)
            if (mean_scaled_p == 0 | mean_scaled_p < 1e-300) {
                severity = max_abs  // Fall back to L∞
            } else {
                severity = max_abs * mean_scaled_p^(1/p)
            }
            term2 = kappa * severity
        } else {
            // Standard computation for moderate p
            term2 = kappa * (mean(abs_x:^p))^(1/p)
        }
    }
    
    return(term1 + term2)
}


// ============================================================================
// PSI DELTA FUNCTION (OVERALL MODE)
// ============================================================================

/**
 * @function _pretest_psi_delta
 * @brief Calculate the psi^Delta function for overall mode critical value computation
 *
 * Implements the overall mode psi function from Appendix C:
 *
 *   psi^Delta(x) = |mean(x_post)| + S(x_pre, p)
 *
 * Key difference from iterative mode: NO kappa multiplier. In overall mode,
 * cumulative violations (nubar) directly bound post-treatment bias without
 * requiring time-weight amplification.
 *
 * @param x      Estimation error vector of dimension (T-1) x 1
 * @param T_pre  Number of pre-treatment periods (t0 - 1)
 * @param T_post Number of post-treatment periods (T - t0 + 1)
 * @param p      Norm index for severity (p >= 1, or >= 1e10 for infinity)
 *
 * @return Psi^Delta value: |mean(x_post)| + severity(x_pre, p)
 *
 * @note No kappa multiplier - this is the key distinction from iterative mode
 * @note Overall violations represent cumulative deviation from parallel trends
 *
 * @see Mikhaeil & Harshaw (2025), Appendix C
 */
real scalar _pretest_psi_delta(real colvector x, real scalar T_pre,
                               real scalar T_post, real scalar p)
{
    real scalar term1, term2, n_pre
    real colvector x_pre, x_post
    
    // Dimension of x should be T-1 = (T_pre-1) + T_post
    n_pre = T_pre - 1
    
    if (rows(x) != n_pre + T_post) {
        _error("_pretest_psi_delta: dimension mismatch - expected " + 
               strofreal(n_pre + T_post) + " but got " + strofreal(rows(x)))
        return(.)
    }
    
    // Edge case: no pre-treatment violations
    if (n_pre < 1) {
        x_post = x[1..T_post]
        return(abs(mean(x_post)))
    }
    
    // Extract x_pre and x_post
    x_pre = x[1..n_pre]
    x_post = x[(n_pre + 1)..rows(x)]
    
    // Term 1: |mean(x_post)|
    term1 = abs(mean(x_post))
    
    // Term 2: severity(x_pre, p) - NO κ multiplier
    // Note: For large p (approaching infinity), use max norm to avoid numerical underflow
    if (p >= 1e10) {
        // p = ∞: use max(|x_pre|)
        term2 = max(abs(x_pre))
    } else {
        // General p-norm with numerical stability
        real colvector abs_x2
        real scalar max_abs2
        
        abs_x2 = abs(x_pre)
        max_abs2 = max(abs_x2)
        
        if (max_abs2 == 0 | max_abs2 < 1e-300) {
            term2 = 0
        } else if (p >= 100) {
            // For large p, use scaling to avoid underflow
            real colvector scaled2
            real scalar mean_scaled_p2
            scaled2 = abs_x2 :/ max_abs2
            mean_scaled_p2 = mean(scaled2:^p)
            if (mean_scaled_p2 == 0 | mean_scaled_p2 < 1e-300) {
                term2 = max_abs2  // Fall back to L∞
            } else {
                term2 = max_abs2 * mean_scaled_p2^(1/p)
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
 * @function _pretest_critical_value
 * @brief Compute Monte Carlo critical value f(alpha, Sigma)
 *
 * Implements the critical value computation from Appendix D:
 *
 *   f(alpha, Sigma) = inf{c : Pr{psi(Z) >= c} <= alpha}, Z ~ N(0, Sigma)
 *
 * Monte Carlo approximation algorithm:
 *   1. Regularize: Sigma_reg = Sigma + epsilon*I (ensure positive definiteness)
 *   2. Cholesky: L = cholesky(Sigma_reg) where Sigma_reg = LL'
 *   3. Sample: W^(s) ~ N(0, I), Z^(s) = L * W^(s) for s = 1,...,S
 *   4. Compute: psi^(s) = psi(Z^(s))
 *   5. Quantile: f = Q_{1-alpha}({psi^(s)})
 *
 * @param Sigma_hat    Asymptotic covariance matrix of dimension (T-1) x (T-1)
 * @param alpha        Significance level (default: 0.05)
 * @param S            Number of Monte Carlo simulations (default: 5000)
 * @param T_pre        Number of pre-treatment periods (t0 - 1)
 * @param T_post       Number of post-treatment periods (T - t0 + 1)
 * @param p            Norm index for severity measure
 * @param kappa        Bias bound constant (from _pretest_kappa)
 * @param overall_mode 0 for iterative (default), 1 for overall mode
 * @param seed         Random seed for reproducibility (default: 12345)
 *
 * @return Critical value f(alpha, Sigma_hat) such that Pr{psi(Z) >= f} <= alpha
 *
 * @note Uses Cholesky decomposition: Z = L*W where W ~ N(0,I) gives Z ~ N(0,Sigma)
 * @note Regularization epsilon = 1e-10 ensures numerical stability
 * @note Quantile uses ceiling for conservative (larger) critical values
 *
 * @see Mikhaeil & Harshaw (2025), Appendix D
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
    
    // Default values for optional arguments
    if (args() < 8 | missing(overall_mode)) {
        overall_mode = 0
    }
    if (args() < 9 | missing(seed)) {
        seed = 12345
    }
    
    // Validate inputs
    dim = rows(Sigma_hat)
    if (dim != cols(Sigma_hat)) {
        _error("_pretest_critical_value: Sigma_hat must be square")
        return(.)
    }
    
    // Expected dimension: T-1 = (T_pre-1) + T_post
    if (dim != (T_pre - 1) + T_post) {
        _error("_pretest_critical_value: Sigma_hat dimension mismatch - expected " + 
               strofreal((T_pre - 1) + T_post) + " but got " + strofreal(dim))
        return(.)
    }
    
    // Set random seed for reproducibility
    rseed(seed)
    
    // Step 1: Regularize covariance matrix
    // Add small positive value to diagonal to ensure positive definiteness
    Sigma_reg = Sigma_hat + 1e-10 * I(dim)
    
    // Step 2: Cholesky decomposition
    // Mata cholesky(A) returns lower triangular L where A = LL'
    L = cholesky(Sigma_reg)
    
    // Step 3-4: Monte Carlo sampling
    psi_values = J(S, 1, .)
    
    for (i = 1; i <= S; i++) {
        // Draw standard normal vector W ~ N(0, I)
        W = rnormal(dim, 1, 0, 1)
        
        // Transform: Z = L × W ~ N(0, Σ_reg)
        Z = L * W
        
        // Compute ψ value
        if (overall_mode == 1) {
            psi_values[i] = _pretest_psi_delta(Z, T_pre, T_post, p)
        } else {
            psi_values[i] = _pretest_psi(Z, T_pre, T_post, p, kappa)
        }
    }
    
    // Step 5: Compute (1-α) quantile
    // Sort ascending
    _sort(psi_values, 1)
    
    // Quantile index: ⌈S(1-α)⌉ (ceiling for conservative estimate)
    quantile_idx = ceil(S * (1 - alpha))
    
    // Ensure index is within bounds
    if (quantile_idx < 1) quantile_idx = 1
    if (quantile_idx > S) quantile_idx = S
    
    f_alpha = psi_values[quantile_idx]
    
    return(f_alpha)
}


// ============================================================================
// CONVENIENCE FUNCTIONS
// ============================================================================

/**
 * @function _pretest_critical_value_simple
 * @brief Simplified interface for critical value with default parameters
 *
 * Convenience wrapper for _pretest_critical_value() using defaults:
 * alpha = 0.05, S = 5000, p = 2, seed = 12345, iterative mode.
 *
 * @param Sigma_hat Asymptotic covariance matrix (T-1) x (T-1)
 * @param T_pre     Number of pre-treatment periods
 * @param T_post    Number of post-treatment periods
 * @param kappa     Bias bound constant
 *
 * @return Critical value f(0.05, Sigma_hat) at 5% significance level
 */
real scalar _pretest_critical_value_simple(real matrix Sigma_hat,
                                           real scalar T_pre,
                                           real scalar T_post,
                                           real scalar kappa)
{
    return(_pretest_critical_value(Sigma_hat, 0.05, 5000, T_pre, T_post, 2, kappa, 0, 12345))
}


/**
 * @function _pretest_verify_cholesky
 * @brief Verify Cholesky decomposition correctness (for testing)
 *
 * Diagnostic function that verifies L * L' = Sigma within tolerance.
 * Used for unit testing and validating numerical stability.
 *
 * @param Sigma Covariance matrix to decompose
 * @param tol   Tolerance for comparison (default: 1e-10)
 *
 * @return 1 if max|LL' - Sigma| < tol, 0 otherwise
 */
real scalar _pretest_verify_cholesky(real matrix Sigma, | real scalar tol)
{
    real matrix L, Sigma_reconstructed
    real scalar max_error
    
    if (args() < 2 | missing(tol)) {
        tol = 1e-10
    }
    
    // Add small regularization
    L = cholesky(Sigma + 1e-10 * I(rows(Sigma)))
    
    // Reconstruct
    Sigma_reconstructed = L * L'
    
    // Check max absolute difference
    max_error = max(abs(Sigma_reconstructed - Sigma))
    
    return(max_error < tol)
}


end
