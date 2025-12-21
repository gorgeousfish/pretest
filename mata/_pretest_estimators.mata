*! _pretest_estimators.mata v0.1.0
*! Difference-in-Differences Estimators for Conditional Extrapolation Pre-Test
*!
*! Description:
*!   Implements the core DID estimators for the conditional extrapolation
*!   pre-test framework using the potential outcomes approach (Abadie, 2005;
*!   Egami & Yamauchi, 2023).
*!
*! Functions:
*!   _pretest_group_mean()  - Period-group sample mean: Y_bar_{t,d}
*!   _pretest_did_t()       - Single-period DID estimator: delta_hat_t
*!   _pretest_did_vector()  - DID vector: (delta_{t0}, ..., delta_T)'
*!   _pretest_did_avg()     - Average DID estimator: delta_bar
*!
*! Mathematical Framework (Section 2.1):
*!   The DID estimand for post-treatment period t >= t0:
*!
*!     delta_t = E[Y_t - Y_{t0} | D=1] - E[Y_t - Y_{t0} | D=0]
*!
*!   Sample analog (unbiased estimator):
*!
*!     delta_hat_t = (Y_bar_{t,1} - Y_bar_{t0,1}) - (Y_bar_{t,0} - Y_bar_{t0,0})
*!
*!   Average post-treatment DID:
*!
*!     delta_bar = (1/T_post) * sum_{t=t0}^T delta_t
*!
*!   Under parallel trends (nu_t = 0 for all t), delta_bar = tau_bar (ATT).
*!
*! Reference:
*!   Mikhaeil, J. M. and C. Harshaw. 2025. In Defense of the Pre-Test: Valid
*!   Inference when Testing Violations of Parallel Trends for Difference-in-
*!   Differences. arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470

version 17.0
mata:

// ============================================================================
// PERIOD-GROUP SAMPLE MEAN
// ============================================================================

/**
 * @function _pretest_group_mean
 * @brief Compute period-group sample mean
 *
 * Calculates the sample mean of outcomes for a specific treatment group
 * at a specific time period:
 *
 *   Y_bar_{t,d} = (1/n_{td}) * sum_{i: D_i=d, t_i=t} Y_i
 *
 * @param Y    Outcome variable (n x 1 column vector)
 * @param D    Treatment indicator (n x 1), values in {0, 1}
 * @param time Time period index (n x 1 column vector)
 * @param t    Target time period (scalar)
 * @param d    Target treatment group (0 = control, 1 = treated)
 *
 * @return Sample mean Y_bar_{t,d}, or missing (.) if no valid observations
 *
 * @see Mikhaeil & Harshaw (2025), Section 2.1
 */
real scalar _pretest_group_mean(real colvector Y, real colvector D,
                                 real colvector time, real scalar t,
                                 real scalar d) {
    real colvector selection, Y_selected
    
    // Select observations satisfying: time==t & D==d & Y not missing
    selection = (time :== t) :& (D :== d) :& !rowmissing(Y)
    Y_selected = select(Y, selection)
    
    // Check for valid observations
    if (rows(Y_selected) == 0) {
        errprintf("Warning: No observations for time=%f, D=%f\n", t, d)
        return(.)
    }
    
    return(mean(Y_selected))
}

// ============================================================================
// SINGLE-PERIOD DID ESTIMATOR
// ============================================================================

/**
 * @function _pretest_did_t
 * @brief Compute single-period difference-in-differences estimator
 *
 * Implements the DID estimator for post-treatment period t (Section 2.1):
 *
 * @population_estimand
 *   delta_t = E[Y_t - Y_{t0} | D=1] - E[Y_t - Y_{t0} | D=0]
 *
 * @sample_analog
 *   delta_hat_t = (Y_bar_{t,1} - Y_bar_{t0,1}) - (Y_bar_{t,0} - Y_bar_{t0,0})
 *
 * @identification
 *   delta_t identifies tau_t (ATT at time t) if and only if nu_t = 0 (parallel trends).
 *   Under violations: tau_t = delta_t + nubar_t, where nubar_t is cumulative violation.
 *
 * @param Y    Outcome vector (n x 1)
 * @param D    Treatment indicator (n x 1), binary D in {0, 1}
 * @param time Time period index (n x 1), mapped to consecutive integers 1:T
 * @param t    Target period (t >= t0 for post-treatment periods)
 * @param t0   Treatment time (first post-treatment period, reference for DID)
 *
 * @return delta_hat_t (DID estimate), or missing (.) if insufficient observations
 *
 * @note delta_{t0} = 0 by construction since Y_{t0} - Y_{t0} = 0 for all groups
 * @note This estimator is unbiased: E[delta_hat_t] = delta_t under random sampling
 *
 * @see Mikhaeil & Harshaw (2025), Section 2.1
 */
real scalar _pretest_did_t(real colvector Y, real colvector D,
                            real colvector time, real scalar t,
                            real scalar t0) {
    real scalar mean_Y_t_1, mean_Y_t0_1, mean_Y_t_0, mean_Y_t0_0
    
    // Treatment group means
    mean_Y_t_1 = _pretest_group_mean(Y, D, time, t, 1)
    mean_Y_t0_1 = _pretest_group_mean(Y, D, time, t0, 1)
    
    // Control group means
    mean_Y_t_0 = _pretest_group_mean(Y, D, time, t, 0)
    mean_Y_t0_0 = _pretest_group_mean(Y, D, time, t0, 0)
    
    // Check for missing values
    if (missing(mean_Y_t_1) | missing(mean_Y_t0_1) | 
        missing(mean_Y_t_0) | missing(mean_Y_t0_0)) {
        return(.)
    }
    
    // Paper formula: delta_t = (Y_t - Y_{t0})_treat - (Y_t - Y_{t0})_control
    return((mean_Y_t_1 - mean_Y_t0_1) - (mean_Y_t_0 - mean_Y_t0_0))
}

// ============================================================================
// DID VECTOR FOR POST-TREATMENT PERIODS
// ============================================================================

/**
 * @function _pretest_did_vector
 * @brief Compute DID vector for all post-treatment periods
 *
 * Constructs the vector of post-treatment DID estimates:
 *
 *   delta = (delta_{t0}, delta_{t0+1}, ..., delta_T)'
 *
 * This vector forms the last T_post elements of the parameter vector theta
 * as defined in Section 2.2, Equation (205).
 *
 * @param Y    Outcome variable (n x 1 column vector)
 * @param D    Treatment indicator (n x 1), values in {0, 1}
 * @param time Time period index (n x 1 column vector)
 * @param t0   Treatment time (first post-treatment period)
 * @param T    Total number of time periods
 *
 * @return Column vector of dimension T_post x 1, where T_post = T - t0 + 1
 *         - delta[1] = delta_{t0} (equals 0 by construction)
 *         - delta[k] = delta_{t0+k-1} for k = 1, ..., T_post
 *
 * @see Mikhaeil & Harshaw (2025), Section 2.2, Equation (205)
 */
real colvector _pretest_did_vector(real colvector Y, real colvector D,
                                    real colvector time, real scalar t0,
                                    real scalar T) {
    real colvector delta
    real scalar t, idx, T_post
    
    // Compute T_post = T - t0 + 1 (includes t0)
    T_post = T - t0 + 1
    
    // Initialize vector
    delta = J(T_post, 1, .)
    
    // Compute DID for each period
    idx = 1
    for (t = t0; t <= T; t++) {
        delta[idx] = _pretest_did_t(Y, D, time, t, t0)
        idx++
    }
    
    return(delta)
}

// ============================================================================
// AVERAGE POST-TREATMENT DID ESTIMATOR
// ============================================================================

/**
 * @function _pretest_did_avg
 * @brief Compute average post-treatment DID estimator
 *
 * Calculates the average DID across all post-treatment periods (Section 2.1):
 *
 *   delta_bar = (1/T_post) * sum_{t=t0}^T delta_t
 *
 * Under parallel trends, delta_bar is a consistent estimator for the
 * average post-treatment ATT (tau_bar).
 *
 * @param delta DID vector from _pretest_did_vector() (T_post x 1)
 *
 * @return Average DID estimator delta_bar_hat, or missing (.) if
 *         vector is empty or contains all missing values
 *
 * @see Mikhaeil & Harshaw (2025), Section 2.1
 */
real scalar _pretest_did_avg(real colvector delta) {
    if (rows(delta) == 0) return(.)
    
    // Compute mean excluding missing values
    if (missing(delta)) {
        // Check if all missing
        if (sum(!rowmissing(delta)) == 0) return(.)
    }
    
    return(mean(delta))
}

// ============================================================================
// AUXILIARY FUNCTIONS
// ============================================================================

/**
 * @function _pretest_get_group_means
 * @brief Compute all period-group sample means
 *
 * Calculates sample means for all period-treatment group combinations.
 * Useful for diagnostics and event study plots.
 *
 * @param Y    Outcome variable (n x 1 column vector)
 * @param D    Treatment indicator (n x 1), values in {0, 1}
 * @param time Time period index (n x 1 column vector)
 * @param T    Total number of time periods
 *
 * @return Matrix of dimension T x 2:
 *         - Column 1: control group means (Y_bar_{t, D=0})
 *         - Column 2: treatment group means (Y_bar_{t, D=1})
 */
real matrix _pretest_get_group_means(real colvector Y, real colvector D,
                                      real colvector time, real scalar T) {
    real matrix means
    real scalar t
    
    means = J(T, 2, .)
    
    for (t = 1; t <= T; t++) {
        means[t, 1] = _pretest_group_mean(Y, D, time, t, 0)  // Control
        means[t, 2] = _pretest_group_mean(Y, D, time, t, 1)  // Treatment
    }
    
    return(means)
}

/**
 * @function _pretest_count_by_group
 * @brief Count observations by period and treatment group
 *
 * Returns observation counts for each period-treatment group combination.
 * Useful for checking data balance and detecting missing cells.
 *
 * @param D    Treatment indicator (n x 1), values in {0, 1}
 * @param time Time period index (n x 1 column vector)
 * @param T    Total number of time periods
 *
 * @return Matrix of dimension T x 2:
 *         - Column 1: control group counts (n_{t,0})
 *         - Column 2: treatment group counts (n_{t,1})
 */
real matrix _pretest_count_by_group(real colvector D, real colvector time,
                                     real scalar T) {
    real matrix counts
    real colvector selection
    real scalar t, d
    
    counts = J(T, 2, 0)
    
    for (t = 1; t <= T; t++) {
        for (d = 0; d <= 1; d++) {
            selection = (time :== t) :& (D :== d)
            counts[t, d+1] = sum(selection)
        }
    }
    
    return(counts)
}

end
