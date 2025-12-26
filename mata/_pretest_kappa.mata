*! _pretest_kappa.mata
*! Bias Amplification Factor for Conditional Extrapolation Inference
*!
*! Computes the kappa constant from Proposition 1 of Mikhaeil & Harshaw (2025),
*! which provides a sharp upper bound on the bias |tau_bar - delta_bar| under
*! the conditional extrapolation assumption:
*!
*!     |tau_bar - delta_bar| <= kappa * S_pre
*!
*! The kappa constant measures how iterative violations accumulate into bias:
*!
*!     kappa = ( (1/T_post) * sum_{t=1}^{T_post} t^q )^{1/q}
*!
*! where q is the Hoelder conjugate of p, satisfying 1/p + 1/q = 1.
*!
*! Reference:
*!   Mikhaeil, J. M. and C. Harshaw. 2025. In Defense of the Pre-Test: Valid
*!   Inference when Testing Violations of Parallel Trends for Difference-in-
*!   Differences. arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470

version 17.0

mata:
mata set matastrict on

// ---------------------------------------------------------------------------
// _pretest_kappa: Bias amplification factor (Proposition 1)
// ---------------------------------------------------------------------------
//
// Arguments:
//     T_post       Number of post-treatment periods (>= 1)
//     p            Lp-norm exponent for severity (>= 1; use 1e10 for infinity)
//     overall_mode 0 = iterative violations (default), 1 = overall violations
//
// Returns:
//     kappa >= 1, or missing (.) if inputs invalid
//     In overall mode, returns 1 (Appendix C)
//
// Special cases:
//     p = 1:   kappa = T_post           (worst-case accumulation)
//     p = 2:   kappa = sqrt(mean(t^2))  (Euclidean, self-dual)
//     p = inf: kappa = (T_post + 1) / 2 (tightest bound)
//
// Properties:
//     - kappa >= 1, with equality when T_post = 1
//     - Strictly increasing in T_post
//     - Strictly decreasing in p
// ---------------------------------------------------------------------------
real scalar _pretest_kappa(real scalar T_post, | real scalar p, 
                           real scalar overall_mode)
{
    real scalar q, kappa
    real colvector t_vec
    
    // Default parameter values
    if (args() < 2 | missing(p)) p = 2
    if (args() < 3 | missing(overall_mode)) overall_mode = 0
    
    // Overall mode: kappa = 1 (Appendix C)
    if (overall_mode == 1) return(1)
    
    // Input validation
    if (T_post < 1 | p < 1) return(.)
    
    // p = 1: Hoelder conjugate q = infinity, kappa = max(t) = T_post
    if (p == 1) return(T_post)
    
    // p = infinity: Hoelder conjugate q = 1, kappa = mean(t) = (T_post + 1) / 2
    if (p >= 1e10) return((T_post + 1) / 2)
    
    // General case: q = p / (p - 1), kappa = (mean(t^q))^(1/q)
    q = p / (p - 1)
    t_vec = 1::T_post
    kappa = (mean(t_vec:^q))^(1/q)
    
    return(kappa)
}

end