*! _pretest_estimators.mata
*! Difference-in-Differences Estimators for Pre-Test Framework
*!
*! Description:
*!   Implements non-parametric sample mean estimators for the difference-in-
*!   differences (DID) framework under the potential outcomes approach.
*!   These estimators form the identified component of the parameter vector
*!   theta = (nu_2, ..., nu_{t0-1}, delta_{t0}, ..., delta_T).
*!
*! Key Estimands (Mikhaeil & Harshaw, 2025, Section 2.1):
*!
*!   Single-period DID:
*!     delta_t = E[Y_t - Y_{t0} | D=1] - E[Y_t - Y_{t0} | D=0]
*!
*!   Average post-treatment DID:
*!     delta_bar = (1/T_post) * sum_{t=t0}^T delta_t
*!
*! Identification:
*!   Under parallel trends (nu_t = 0 for all t), delta_bar identifies the
*!   average treatment effect on the treated (ATT), tau_bar.
*!   Under violations: |tau_bar - delta_bar| <= kappa * S_pre (Proposition 1).
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
// PERIOD-GROUP SAMPLE MEAN
// ============================================================================

/**
 * @brief Compute period-group sample mean
 *
 * Calculates the sample mean of outcomes for a specific treatment group
 * at a specific time period. This is the fundamental building block for
 * DID estimation.
 *
 * Formula:
 *   Y_bar_{t,d} = (1/n_{td}) * sum_{i: D_i=d, t_i=t} Y_i
 *
 * where n_{td} is the number of observations with treatment status d at time t.
 *
 * @param Y    Outcome variable (n x 1 column vector)
 * @param D    Treatment indicator (n x 1), binary values in {0, 1}
 * @param time Time period index (n x 1 column vector)
 * @param t    Target time period (scalar)
 * @param d    Target treatment group: 0 = control, 1 = treated
 *
 * @return Scalar: sample mean Y_bar_{t,d}, or missing (.) if no valid observations
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
 * @brief Compute single-period difference-in-differences estimator
 *
 * Implements the DID estimator for period t relative to reference period t0.
 *
 * Population estimand (Mikhaeil & Harshaw, 2025, Section 2.1):
 *   delta_t = E[Y_t - Y_{t0} | D=1] - E[Y_t - Y_{t0} | D=0]
 *
 * Sample analog:
 *   delta_hat_t = (Y_bar_{t,1} - Y_bar_{t0,1}) - (Y_bar_{t,0} - Y_bar_{t0,0})
 *
 * Identification:
 *   Under parallel trends, delta_t = tau_t (ATT at time t).
 *   Under violations: tau_t = delta_t + nu_bar_t, where nu_bar_t is the
 *   overall (cumulative) violation of parallel trends.
 *
 * @param Y    Outcome variable (n x 1 column vector)
 * @param D    Treatment indicator (n x 1), binary values in {0, 1}
 * @param time Time period index (n x 1 column vector)
 * @param t    Target period (typically t >= t0 for post-treatment)
 * @param t0   Treatment time (reference period for DID computation)
 *
 * @return Scalar: delta_hat_t, or missing (.) if insufficient observations
 *
 * @note delta_{t0} = 0 by construction since Y_{t0} - Y_{t0} = 0
 * @note Unbiased under random sampling: E[delta_hat_t] = delta_t
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
 * @brief Compute DID vector for all post-treatment periods
 *
 * Constructs the vector of post-treatment DID estimates:
 *   delta = (delta_{t0}, delta_{t0+1}, ..., delta_T)'
 *
 * This vector forms the last T_post elements of the parameter vector theta
 * defined in Mikhaeil & Harshaw (2025), Section 2.2:
 *   theta = (nu_2, ..., nu_{t0-1}, delta_{t0}, ..., delta_T)
 *
 * @param Y    Outcome variable (n x 1 column vector)
 * @param D    Treatment indicator (n x 1), binary values in {0, 1}
 * @param time Time period index (n x 1 column vector)
 * @param t0   Treatment time (first post-treatment period)
 * @param T    Total number of time periods
 *
 * @return Column vector (T_post x 1), where T_post = T - t0 + 1
 *         Element k corresponds to delta_{t0+k-1} for k = 1, ..., T_post
 *
 * @note delta_{t0} = 0 by construction (first element)
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
 * @brief Compute average post-treatment DID estimator
 *
 * Calculates the average DID across all post-treatment periods.
 *
 * Formula (Mikhaeil & Harshaw, 2025, Section 2.1):
 *   delta_bar = (1/T_post) * sum_{t=t0}^T delta_t
 *
 * Identification:
 *   Under parallel trends, delta_bar is a consistent estimator for the
 *   average post-treatment ATT (tau_bar). Under violations, the bias is
 *   bounded by kappa * S_pre per Proposition 1.
 *
 * @param delta DID vector from _pretest_did_vector() (T_post x 1)
 *
 * @return Scalar: delta_bar_hat, or missing (.) if vector is empty or
 *         contains all missing values
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
 * @brief Compute all period-group sample means
 *
 * Calculates sample means for all period-treatment group combinations.
 * Useful for diagnostics, data exploration, and event study plots.
 *
 * @param Y    Outcome variable (n x 1 column vector)
 * @param D    Treatment indicator (n x 1), binary values in {0, 1}
 * @param time Time period index (n x 1 column vector)
 * @param T    Total number of time periods
 *
 * @return Matrix (T x 2):
 *         Column 1: control group means Y_bar_{t,0} for t = 1, ..., T
 *         Column 2: treatment group means Y_bar_{t,1} for t = 1, ..., T
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
 * @brief Count observations by period and treatment group
 *
 * Returns observation counts for each period-treatment group combination.
 * Useful for checking data balance, detecting missing cells, and
 * computing standard errors.
 *
 * @param D    Treatment indicator (n x 1), binary values in {0, 1}
 * @param time Time period index (n x 1 column vector)
 * @param T    Total number of time periods
 *
 * @return Matrix (T x 2):
 *         Column 1: control group counts n_{t,0} for t = 1, ..., T
 *         Column 2: treatment group counts n_{t,1} for t = 1, ..., T
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
