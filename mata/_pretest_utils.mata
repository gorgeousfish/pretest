*! _pretest_utils.mata v0.1.0
*! Utility Functions for Pretest Package
*!
*! Description:
*!   Provides fundamental utility functions supporting the pretest Stata package.
*!   These utilities handle numerical stability, input validation, and time
*!   variable mapping for the conditional extrapolation pre-test framework.
*!
*! Functions:
*!   _pretest_check_dim()          - Validate vector/matrix dimensions
*!   _pretest_safe_div()           - Division with zero-denominator protection
*!   _pretest_lp_norm()            - Normalized Lp-norm computation
*!   _pretest_make_index_map()     - Map non-consecutive time values to indices
*!   _pretest_find_time_index()    - Lookup consecutive index for time value
*!   _pretest_is_consecutive()     - Test if time variable has gaps
*!   _pretest_validate_positive()  - Validate strictly positive scalar input
*!   _pretest_conjugate_exponent() - Compute Holder conjugate q from p
*!
*! Naming Convention:
*!   All internal functions use _pretest_ prefix to avoid namespace conflicts
*!   with user programs and other Stata packages.
*!
*! Reference:
*!   Mikhaeil, J. M. and C. Harshaw. 2025. In Defense of the Pre-Test: Valid
*!   Inference when Testing Violations of Parallel Trends for Difference-in-
*!   Differences. arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470

version 17.0
mata:

// ============================================================================
// DIMENSION VALIDATION
// ============================================================================

/**
 * @function _pretest_check_dim
 * @brief Validate vector dimensions
 *
 * Verifies that input vector dimension matches the expected size.
 * Used for input validation before vector operations.
 *
 * @param v            Column vector to check
 * @param expected_dim Expected number of rows
 *
 * @return 1 if dimension matches, 0 if mismatch (prints error message)
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

/**
 * @function _pretest_safe_div
 * @brief Perform safe division avoiding division by zero
 *
 * Returns num/denom if |denom| >= 1e-15, otherwise returns missing (.).
 * Threshold 1e-15 approximates IEEE 754 double precision machine epsilon.
 *
 * @param num   Numerator
 * @param denom Denominator
 *
 * @return num/denom or missing (.) if denominator is near zero
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

/**
 * @function _pretest_lp_norm
 * @brief Compute normalized Lp norm of a vector
 *
 * Calculates the normalized Lp norm:
 *
 *   ||v||_p = ((1/n) * sum(|v_i|^p))^{1/p}
 *
 * This matches the severity measure definition in Section 3.1:
 *   S_pre = ((1/(T_pre-1)) * sum_{t=2}^{t0-1} |nu_t|^p)^{1/p}
 *
 * @param v Input column vector
 * @param p Norm index (p >= 1), or missing/1e10 for L-infinity
 *
 * @return Normalized Lp norm:
 *         - p = 1: mean absolute value
 *         - p = 2: root mean square (RMS)
 *         - p = inf: maximum absolute value
 *         - Empty vector returns 0
 *
 * @see Mikhaeil & Harshaw (2025), Section 3.1
 */
real scalar _pretest_lp_norm(real colvector v, real scalar p) {
    real scalar n, result
    
    n = rows(v)
    if (n == 0) return(0)
    
    // Input validation: Paper requires p >= 1 (Section 3.1)
    // Lp norm is mathematically invalid for p < 1 (violates triangle inequality)
    if (p < 1 & p != .) {
        errprintf("Error: p must be >= 1 for valid Lp norm (got %f)\n", p)
        errprintf("  Reference: Section 3.1, severity measure definition\n")
        return(.)
    }
    
    // L-infinity norm special case
    if (p == . | p >= 1e10) {
        return(max(abs(v)))
    }
    
    // Standard Lp norm computation
    // (1/n * sum(|v_i|^p))^(1/p)
    result = (mean(abs(v):^p))^(1/p)
    
    return(result)
}

// ============================================================================
// TIME INDEX MAPPING
// ============================================================================

/**
 * @function _pretest_make_index_map
 * @brief Create time-to-index mapping for non-consecutive time values
 *
 * Maps original time values to consecutive indices {1, 2, ..., T}.
 * Required when the time variable has gaps (e.g., 2000, 2002, 2004).
 *
 * @param orig_times Sorted unique original time values (T x 1)
 *
 * @return Matrix of dimension T x 2:
 *         - Column 1: original time values
 *         - Column 2: consecutive indices 1, 2, ..., T
 *
 * @example Input: (2000, 2002, 2004, 2006)' -> Output: ((2000,1), (2002,2), ...)
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

/**
 * @function _pretest_find_time_index
 * @brief Find consecutive index for an original time value
 *
 * Looks up the consecutive index corresponding to an original time value
 * using the mapping matrix from _pretest_make_index_map().
 *
 * @param index_map Mapping matrix from _pretest_make_index_map()
 * @param orig_time Original time value to look up
 *
 * @return Corresponding consecutive index, or missing (.) if not found
 */
real scalar _pretest_find_time_index(real matrix index_map, real scalar orig_time) {
    real scalar i, T
    
    T = rows(index_map)
    for (i = 1; i <= T; i++) {
        if (index_map[i, 1] == orig_time) {
            return(index_map[i, 2])
        }
    }
    
    // Not found
    return(.)
}

/**
 * @function _pretest_is_consecutive
 * @brief Check if time values form a consecutive sequence
 *
 * Determines whether time values are consecutive integers without gaps.
 * Used to decide if time-to-index mapping is needed.
 *
 * @param times Sorted time value column vector
 *
 * @return 1 if consecutive (no gaps), 0 if has gaps
 */
real scalar _pretest_is_consecutive(real colvector times) {
    real scalar T, t_min, t_max
    
    T = rows(times)
    if (T <= 1) return(1)
    
    t_min = min(times)
    t_max = max(times)
    
    // If consecutive, then t_max - t_min + 1 == T
    return((t_max - t_min + 1) == T)
}

// ============================================================================
// INPUT VALIDATION
// ============================================================================

/**
 * @function _pretest_validate_positive
 * @brief Validate that a scalar is positive
 *
 * Checks whether the input scalar is strictly positive (> 0).
 * Prints an error message if validation fails.
 *
 * @param x    Scalar to validate
 * @param name Variable name for error message
 *
 * @return 1 if positive, 0 if non-positive or missing
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

/**
 * @function _pretest_conjugate_exponent
 * @brief Compute Holder conjugate exponent q
 *
 * Computes the conjugate exponent q satisfying 1/p + 1/q = 1.
 * Used in kappa constant computation (Proposition 1, Section 3.2):
 *
 *   kappa = ((1/T_post) * sum_{t=1}^{T_post} t^q)^{1/q}
 *
 * @param p Norm index (p >= 1)
 *
 * @return Conjugate exponent q = p/(p-1):
 *         - p = 1: returns missing (.) representing infinity
 *         - p = infinity: returns 1
 *
 * @see Mikhaeil & Harshaw (2025), Section 3.2, Proposition 1
 */
real scalar _pretest_conjugate_exponent(real scalar p) {
    // p = infinity case
    if (p == . | p >= 1e10) {
        return(1)
    }
    
    // p = 1 case
    if (abs(p - 1) < 1e-10) {
        return(.)  // Represents infinity
    }
    
    // General case: q = p/(p-1)
    return(p / (p - 1))
}

end
