*! _pretest_kappa.mata v0.1.0
*! Kappa Constant for Bias Bound
*!
*! Implements kappa constant calculation for the conditional extrapolation
*! pre-test bias bounds.
*!
*! Reference:
*!   Mikhaeil & Harshaw (2025), "In Defense of the Pre-Test"
*!   arXiv:2510.26470, Section 3.2, Proposition 1

version 17.0

mata:

// ============================================================================
// KAPPA CONSTANT
// ============================================================================

/**
 * @function _pretest_kappa
 * @brief Calculate bias bound constant kappa for the conditional extrapolation framework
 *
 * Computes the kappa constant from Proposition 1 (Section 3.2) that bounds the
 * extrapolation of pre-treatment violations to post-treatment bias:
 *
 *   kappa = ((1/T_post) * sum_{t=1}^{T_post} t^q)^{1/q}
 *
 * where q = p/(p-1) is the Holder conjugate exponent satisfying 1/p + 1/q = 1.
 *
 * Special cases:
 *   - p = 1:   kappa = T_post (q = infinity, max norm)
 *   - p = 2:   kappa = sqrt((T_post + 1) * (2*T_post + 1) / 6) (closed-form)
 *   - p = inf: kappa = (T_post + 1) / 2 (q = 1, sum/mean)
 *
 * In overall mode, kappa = 1 regardless of p, as overall violations directly
 * bound bias without time-weight amplification.
 *
 * @param T_post       Number of post-treatment periods (T - t0 + 1)
 * @param p            Norm index for severity measure (p >= 1, or missing/1e10 for infinity)
 * @param overall_mode 0 for iterative violations (default), 1 for overall violations
 *
 * @return Bias bound constant kappa, or missing (.) if T_post < 1 or p < 1
 *
 * @note Derived from Holder's inequality: |sum(a_t * b_t)| <= ||a||_p * ||b||_q
 * @note kappa weights how iterative violations accumulate over post-treatment periods
 *
 * @see Mikhaeil & Harshaw (2025), Section 3.2, Proposition 1
 */
real scalar _pretest_kappa(real scalar T_post, | real scalar p, 
                           real scalar overall_mode)
{
    real scalar q, kappa
    real colvector t_vec
    
    // Default values
    if (args() < 2 | missing(p)) {
        p = 2
    }
    if (args() < 3 | missing(overall_mode)) {
        overall_mode = 0
    }
    
    // Overall mode: κ = 1
    if (overall_mode == 1) {
        return(1)
    }
    
    // Validate T_post
    if (T_post < 1) {
        return(.)
    }
    
    // Validate p (must be >= 1 for valid norm)
    if (p < 1) {
        return(.)
    }
    
    // p = 1 (q = infinity): kappa = max(t) = T_post  
    if (p == 1) {
        return(T_post)
    }
    
    // p = infinity (q = 1): kappa = (1/T_post) * sum(t) = (T_post+1)/2
    if (p >= 1e10) {
        return((T_post + 1) / 2)
    }
    
    // General case: compute ((1/T_post) * sum_{t=1}^{T_post} t^q)^(1/q)
    // From Holder's inequality: 1/p + 1/q = 1 => q = p/(p-1)
    q = p / (p - 1)
    t_vec = 1::T_post
    kappa = (mean(t_vec:^q))^(1/q)
    
    return(kappa)
}

end