*! _pretest_kappa.mata v0.1.0
*! Bias Bound Constant Kappa for Conditional Extrapolation Framework
*!
*! Description:
*!   Computes the kappa constant that bounds the bias between the average
*!   treatment effect on the treated (ATT) and the difference-in-differences
*!   estimand under the conditional extrapolation assumption.
*!
*! Mathematical Foundation (Proposition 1, Section 3.2):
*!   Under Assumption 3 (Conditional Extrapolation), if S_pre <= M, then:
*!
*!       |tau_bar - delta_bar| <= kappa * S_pre
*!
*!   where kappa = ((1/T_post) * sum_{t=1}^{T_post} t^q)^{1/q}
*!   and q satisfies the Holder conjugate relation: 1/p + 1/q = 1.
*!
*! Reference:
*!   Mikhaeil, J. M. and C. Harshaw. 2025. In Defense of the Pre-Test: Valid
*!   Inference when Testing Violations of Parallel Trends for Difference-in-
*!   Differences. arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470

version 17.0

mata:

// ============================================================================
// KAPPA CONSTANT
// ============================================================================

/**
 * @function _pretest_kappa
 * @brief Compute bias bound constant kappa (Proposition 1)
 *
 * Computes the kappa constant from Proposition 1 (Section 3.2) that establishes
 * a sharp upper bound on the bias between the causal estimand (tau_bar) and the
 * identified estimand (delta_bar) under the conditional extrapolation assumption.
 *
 * @formula
 *   kappa = ((1/T_post) * sum_{t=1}^{T_post} t^q)^{1/q}
 *
 *   where q = p/(p-1) is the Holder conjugate exponent satisfying 1/p + 1/q = 1.
 *
 * @derivation
 *   The bound arises from applying Holder's inequality to the bias decomposition:
 *     |tau_bar - delta_bar| = |mean(nu_{t0},...,nu_T)|
 *                           <= ||weights||_q * ||nu||_p
 *                           = kappa * S_post
 *                           <= kappa * S_pre  (by Assumption 3)
 *
 * @special_cases
 *   p = 1:   kappa = T_post (q = inf, weights bounded by max weight = T_post)
 *   p = 2:   kappa = sqrt((T_post+1)(2*T_post+1)/6) (closed-form RMS of 1:T_post)
 *   p = inf: kappa = (T_post+1)/2 (q = 1, arithmetic mean of 1:T_post)
 *
 * @overall_mode
 *   When overall_mode = 1, returns kappa = 1. In overall mode, cumulative violations
 *   (nubar_t) directly bound the post-treatment bias without requiring the time-weight
 *   amplification factor, as shown in Appendix C.
 *
 * @param T_post       Number of post-treatment periods: T_post = T - t0 + 1
 * @param p            Lp-norm index for severity measure (p >= 1; use 1e10 for L-infinity)
 * @param overall_mode 0 = iterative violations (default), 1 = overall violations
 *
 * @return kappa >= 1 (bias bound constant), or missing (.) if inputs invalid
 *
 * @note Proposition 1 establishes this bound is sharp (attainable in worst case)
 * @note The kappa constant is deterministic given (T_post, p) and does not depend on data
 *
 * @see Mikhaeil & Harshaw (2025), Section 3.2, Proposition 1
 * @see Holder's inequality: |<a,b>| <= ||a||_p * ||b||_q for 1/p + 1/q = 1
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
    
    // Special case: p = 1 implies q = infinity (L-infinity dual norm)
    // kappa = max_{t in 1:T_post}(t) = T_post
    if (p == 1) {
        return(T_post)
    }
    
    // Special case: p = infinity implies q = 1 (L1 dual norm)
    // kappa = (1/T_post) * sum_{t=1}^{T_post} t = (T_post + 1) / 2
    if (p >= 1e10) {
        return((T_post + 1) / 2)
    }
    
    // General case: compute kappa via Holder conjugate exponent
    // Holder relation: 1/p + 1/q = 1  =>  q = p / (p - 1)
    // Formula: kappa = ((1/T_post) * sum_{t=1}^{T_post} t^q)^{1/q}
    q = p / (p - 1)
    t_vec = 1::T_post
    kappa = (mean(t_vec:^q))^(1/q)
    
    return(kappa)
}

end