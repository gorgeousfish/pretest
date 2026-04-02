*! _pretest_utils.mata
*! Utility functions for the pretest package
*!
*! Part of pretest: Conditional Extrapolation Pre-Test for Difference-in-Differences
*!
*! Description:
*!   Core utility functions providing numerical stability, input validation,
*!   and time variable mapping for the conditional extrapolation framework.
*!
*! Contents:
*!   _pretest_check_dim()          - Vector dimension validation
*!   _pretest_safe_div()           - Safe division with zero-denominator handling
*!   _pretest_lp_norm()            - Normalized Lp-norm for severity measure
*!   _pretest_make_index_map()     - Time-to-index mapping for non-consecutive periods
*!   _pretest_find_time_index()    - Index lookup for original time values
*!   _pretest_is_consecutive()     - Check for gaps in time variable
*!   _pretest_validate_positive()  - Strictly positive scalar validation
*!   _pretest_conjugate_exponent() - Hoelder conjugate exponent computation
*!
*! Reference:
*!   Mikhaeil, J. M. and C. Harshaw. 2025. In Defense of the Pre-Test: Valid
*!   Inference when Testing Violations of Parallel Trends for Difference-in-
*!   Differences. arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470

version 17.0

mata:
mata set matastrict on

// ============================================================================
// DIMENSION VALIDATION
// ============================================================================

/*
 * _pretest_check_dim()
 *
 * Validate that a column vector has the expected number of rows.
 *
 * Arguments:
 *   v            - Column vector to validate
 *   expected_dim - Required number of rows
 *
 * Returns:
 *   1 if rows(v) == expected_dim, 0 otherwise (with error message)
 */
real scalar _pretest_check_dim(real colvector v, real scalar expected_dim) {
    if (rows(v) != expected_dim) {
        errprintf("Dimension mismatch: expected %f rows, got %f\n", 
                  expected_dim, rows(v))
        return(0)
    }
    return(1)
}

// ============================================================================
// SAFE ARITHMETIC
// ============================================================================

/*
 * _pretest_safe_div()
 *
 * Perform division with protection against near-zero denominators.
 *
 * Arguments:
 *   num   - Numerator
 *   denom - Denominator
 *
 * Returns:
 *   num/denom if |denom| >= 1e-15, missing value (.) otherwise
 *
 * Note:
 *   Threshold 1e-15 is chosen to approximate IEEE 754 double-precision
 *   machine epsilon, preventing floating-point overflow.
 */
real scalar _pretest_safe_div(real scalar num, real scalar denom) {
    if (abs(denom) < 1e-15) {
        return(.)
    }
    return(num / denom)
}

// ============================================================================
// LP NORM COMPUTATION
// ============================================================================

/*
 * _pretest_lp_norm()
 *
 * Compute the normalized Lp norm of a vector.
 *
 * Implements the severity measure from Mikhaeil & Harshaw (2025), Section 3.1:
 *
 *   S_pre = ( (1/(T_pre-1)) * sum_{t=2}^{t0-1} |nu_t|^p )^{1/p}
 *
 * The normalized form uses mean rather than sum, making the measure
 * comparable across different numbers of pre-treatment periods.
 *
 * Arguments:
 *   v - Column vector of iterative violations (or other values)
 *   p - Norm index, p >= 1; use missing (.) or p >= 1e10 for L-infinity
 *
 * Returns:
 *   Normalized Lp norm: ( mean(|v|^p) )^{1/p}
 *     p = 1:   Mean absolute deviation
 *     p = 2:   Root mean square
 *     p = inf: Maximum absolute value
 *
 * Returns 0 for empty vector, missing (.) if p < 1.
 *
 * Reference:
 *   Mikhaeil & Harshaw (2025), Section 3.1
 */
real scalar _pretest_lp_norm(real colvector v, real scalar p) {
    real scalar n, result
    
    n = rows(v)
    if (n == 0) return(0)
    
    // Validate p >= 1 (required for valid Lp norm; violates triangle inequality otherwise)
    if (p < 1 & p != .) {
        errprintf("Error: p must be >= 1 for valid Lp norm (got %f)\n", p)
        return(.)
    }
    
    // L-infinity norm: max|v_i|
    if (p == . | p >= 1e10) {
        return(max(abs(v)))
    }
    
    // Standard normalized Lp norm: ( (1/n) * sum(|v_i|^p) )^{1/p}
    result = (mean(abs(v):^p))^(1/p)
    
    return(result)
}

// ============================================================================
// TIME INDEX MAPPING
// ============================================================================

/*
 * _pretest_make_index_map()
 *
 * Create a mapping from original time values to consecutive indices.
 *
 * Panel data may have non-consecutive time values (e.g., 2000, 2002, 2004).
 * This function creates a lookup table mapping such values to indices 1, 2, ..., T
 * for internal computations.
 *
 * Arguments:
 *   orig_times - Sorted unique time values (T x 1 column vector)
 *
 * Returns:
 *   T x 2 matrix where column 1 = original times, column 2 = indices 1..T
 *
 * Example:
 *   Input:  (2000, 2002, 2004, 2006)'
 *   Output: (2000, 1 \ 2002, 2 \ 2004, 3 \ 2006, 4)
 */
real matrix _pretest_make_index_map(real colvector orig_times) {
    real scalar T
    real matrix index_map
    
    T = rows(orig_times)
    index_map = J(T, 2, .)
    
    index_map[., 1] = orig_times
    index_map[., 2] = (1::T)
    
    return(index_map)
}

/*
 * _pretest_find_time_index()
 *
 * Look up the consecutive index for an original time value.
 *
 * Arguments:
 *   index_map - Mapping matrix from _pretest_make_index_map()
 *   orig_time - Original time value to look up
 *
 * Returns:
 *   Consecutive index (1..T) if found, missing (.) if not found
 */
real scalar _pretest_find_time_index(real matrix index_map, real scalar orig_time) {
    real scalar i, T
    
    T = rows(index_map)
    for (i = 1; i <= T; i++) {
        if (index_map[i, 1] == orig_time) {
            return(index_map[i, 2])
        }
    }
    
    return(.)
}

/*
 * _pretest_is_consecutive()
 *
 * Check whether time values form a consecutive integer sequence.
 *
 * Used to determine if time-to-index mapping is needed. If times are
 * consecutive (e.g., 1, 2, 3, 4), they can be used directly as indices.
 *
 * Arguments:
 *   times - Sorted time values (column vector)
 *
 * Returns:
 *   1 if consecutive (no gaps), 0 if gaps exist
 */
real scalar _pretest_is_consecutive(real colvector times) {
    real scalar T, t_min, t_max
    
    T = rows(times)
    if (T <= 1) return(1)
    
    t_min = min(times)
    t_max = max(times)
    
    // Consecutive iff range equals count
    return((t_max - t_min + 1) == T)
}

// ============================================================================
// INPUT VALIDATION
// ============================================================================

/*
 * _pretest_validate_positive()
 *
 * Validate that a scalar is strictly positive.
 *
 * Arguments:
 *   x    - Scalar value to validate
 *   name - Variable name for error message
 *
 * Returns:
 *   1 if x > 0, 0 otherwise (with error message)
 */
real scalar _pretest_validate_positive(real scalar x, string scalar name) {
    if (missing(x) | x <= 0) {
        errprintf("%s must be positive (got %f)\n", name, x)
        return(0)
    }
    return(1)
}

// ============================================================================
// CONJUGATE EXPONENT
// ============================================================================

/*
 * _pretest_conjugate_exponent()
 *
 * Compute the Hoelder conjugate exponent q from p.
 *
 * The conjugate exponent satisfies the Hoelder relation: 1/p + 1/q = 1.
 * This is used in the kappa constant computation (Proposition 1):
 *
 *   kappa = ( (1/T_post) * sum_{t=1}^{T_post} t^q )^{1/q}
 *
 * where kappa bounds the bias |tau_bar - delta_bar| <= kappa * S_pre.
 *
 * Arguments:
 *   p - Norm index (p >= 1)
 *
 * Returns:
 *   q = p/(p-1) for p > 1
 *   q = 1       for p = infinity (encoded as missing or >= 1e10)
 *   q = .       for p = 1 (representing infinity)
 *
 * Reference:
 *   Mikhaeil & Harshaw (2025), Section 3.2, Proposition 1
 */
real scalar _pretest_conjugate_exponent(real scalar p) {
    // p = infinity => q = 1
    if (p == . | p >= 1e10) {
        return(1)
    }
    
    // p = 1 => q = infinity (represented as missing)
    if (abs(p - 1) < 1e-10) {
        return(.)
    }
    
    // General case: q = p/(p-1) from 1/p + 1/q = 1
    return(p / (p - 1))
}

end
