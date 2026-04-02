*! _pretest_covariance.mata
*! Asymptotic Covariance Matrix Estimation for the Conditional Extrapolation Pre-Test
*!
*! Description:
*!   Estimates the asymptotic covariance matrix Sigma of the parameter vector
*!   theta = (nu_2, ..., nu_{t0-1}, delta_{t0}, ..., delta_T)' using influence
*!   functions derived from sample mean estimators.
*!
*! Asymptotic Framework:
*!   This module implements variance estimation under the asymptotic normality
*!   framework of Section 2.3:
*!
*!   Assumption 1 (Asymptotic Normality):
*!     sqrt(n) * (theta_hat - theta) -->d N(0, Sigma)
*!
*!   Assumption 2 (Consistent Variance Estimation):
*!     Sigma_hat -->p Sigma
*!
*! Parameter Vector (Section 2.2):
*!   theta = (nu_2, ..., nu_{t0-1}, delta_{t0}, ..., delta_T)'
*!   - nu_t: Iterative parallel trends violations (pre-treatment, t = 2,...,t0-1)
*!   - delta_t: DID estimands (post-treatment, t = t0,...,T)
*!   - Dimension: dim(theta) = (t0-2) + (T-t0+1) = T - 1
*!
*! Estimation Methods:
*!   Standard (iid):
*!     Sigma_hat = (1/(n-1)) * sum_{i=1}^n psi_i * psi_i'
*!
*!   Cluster-robust (Cameron & Miller, 2015):
*!     Sigma_hat = (G/(G-1)) * (1/n) * sum_{g=1}^G u_g * u_g'
*!     where u_g = sum_{i in cluster g} psi_i
*!
*! Reference:
*!   Mikhaeil, J. M. and C. Harshaw. 2025. In Defense of the Pre-Test: Valid
*!   Inference when Testing Violations of Parallel Trends for Difference-in-
*!   Differences. arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470
*!   Section 2.2-2.3, Assumptions 1-2.

version 17.0

mata:
mata set matastrict on

// ============================================================================
// CORE COVARIANCE ESTIMATION
// ============================================================================

/**
 * @function _pretest_covariance
 * @brief Asymptotic covariance matrix estimator for the parameter vector theta
 *
 * Computes the sample covariance matrix Sigma_hat satisfying Assumption 2,
 * such that sqrt(n)(theta_hat - theta) -->d N(0, Sigma) under Assumption 1.
 *
 * The estimation uses influence functions derived from the representation of
 * sample means as M-estimators. For the group-period mean Y_bar_{t,d}, the
 * influence function is:
 *
 *   psi_i(Y_bar_{t,d}) = (n / n_{td}) * 1{D_i=d, t_i=t} * (Y_i - Y_bar_{t,d})
 *
 * where n_{td} is the number of observations in group d at time t.
 *
 * @param Y        Outcome vector (n x 1)
 * @param D        Treatment indicator (n x 1), binary in {0, 1}
 * @param time     Time period vector (n x 1), integers 1 to T
 * @param t0       Treatment time (first post-treatment period)
 * @param T        Total number of periods
 * @param is_panel 1 = panel data, 0 = repeated cross-section (reserved)
 * @param cluster  Optional cluster variable for cluster-robust SEs (n x 1)
 *
 * @return Sigma_hat: Asymptotic covariance matrix of dimension (T-1) x (T-1)
 *
 * @note Matrix is regularized if not positive definite (eigenvalue correction)
 * @note Cluster-robust estimation follows Cameron & Miller (2015)
 *
 * @see Mikhaeil & Harshaw (2025), Section 2.2-2.3, Assumptions 1-2
 */
real matrix _pretest_covariance(real colvector Y, real colvector D,
                                 real colvector time, real scalar t0,
                                 real scalar T, real scalar is_panel,
                                 | real colvector cluster)
{
    real matrix Sigma, psi_mat
    real scalar has_cluster, n, dim
    pragma unused is_panel  // Reserved for future effective sample size adjustment
    
    dim = T - 1
    
    // Check if cluster variable provided
    has_cluster = (args() >= 7 && rows(cluster) > 0)
    
    // Get influence function matrix (n x dim)
    psi_mat = _pretest_influence_matrix(Y, D, time, t0, T)
    
    n = rows(psi_mat)
    
    if (has_cluster) {
        // Cluster-robust covariance
        Sigma = _pretest_covariance_cluster(psi_mat, cluster)
    }
    else {
        // Standard covariance: (1/(n-1)) * psi' * psi
        Sigma = (psi_mat' * psi_mat) / (n - 1)
    }
    
    // Ensure positive definiteness
    Sigma = _pretest_ensure_pd(Sigma)
    
    return(Sigma)
}


// ============================================================================
// INFLUENCE FUNCTION CALCULATION
// ============================================================================

/**
 * @function _pretest_influence_matrix
 * @brief Construct the n x (T-1) influence function matrix for theta
 *
 * Computes the influence function matrix Psi where row i contains
 * observation i's contribution to the asymptotic variance of theta_hat.
 * The covariance is then: Sigma_hat = (1/(n-1)) * Psi' * Psi.
 *
 * Influence Function Derivation:
 *   For a sample mean estimator Y_bar_{t,d} = (1/n_{td}) sum_{i in (t,d)} Y_i,
 *   the influence function under standard asymptotic theory is:
 *
 *     psi_i(Y_bar_{t,d}) = (n/n_{td}) * 1{D_i=d, t_i=t} * (Y_i - Y_bar_{t,d})
 *
 *   The scaling factor n/n_{td} = 1/pi_{td} arises because we need
 *   sqrt(n)(Y_bar_{t,d} - mu_{t,d}) to be O_p(1).
 *
 * For Iterative Violations nu_t (t = 2,...,t0-1):
 *   nu_t = (Y_bar_{t,1} - Y_bar_{t-1,1}) - (Y_bar_{t,0} - Y_bar_{t-1,0})
 *
 *   psi_i(nu_t) = psi_i(Y_bar_{t,1}) - psi_i(Y_bar_{t-1,1})
 *               - psi_i(Y_bar_{t,0}) + psi_i(Y_bar_{t-1,0})
 *
 * For DID Estimands delta_t (t = t0,...,T):
 *   delta_t = (Y_bar_{t,1} - Y_bar_{t0,1}) - (Y_bar_{t,0} - Y_bar_{t0,0})
 *
 *   psi_i(delta_t) = psi_i(Y_bar_{t,1}) - psi_i(Y_bar_{t0,1})
 *                  - psi_i(Y_bar_{t,0}) + psi_i(Y_bar_{t0,0})
 *
 * @param Y    Outcome vector (n x 1)
 * @param D    Treatment indicator (n x 1), binary in {0, 1}
 * @param time Time period vector (n x 1), integers 1 to T
 * @param t0   Treatment time (first post-treatment period)
 * @param T    Total number of periods
 *
 * @return Psi: Influence function matrix of dimension n x (T-1)
 *         - Columns 1 to (t0-2): influence for nu_2, ..., nu_{t0-1}
 *         - Columns (t0-1) to (T-1): influence for delta_{t0}, ..., delta_T
 *
 * @note Column sums are approximately zero (mean-centering property)
 * @note Returns matrix of missing values if any period-group cell is empty
 *
 * @see Mikhaeil & Harshaw (2025), Section 2.2
 */
real matrix _pretest_influence_matrix(real colvector Y, real colvector D,
                                       real colvector time, real scalar t0,
                                       real scalar T)
{
    real matrix psi
    real scalar n, dim, t, i, T_pre
    real scalar mean_t_1, mean_t_0, mean_tm1_1, mean_tm1_0, mean_t0_1, mean_t0_0
    real scalar n_t_1, n_t_0, n_tm1_1, n_tm1_0, n_t0_1, n_t0_0
    real scalar col_idx
    real scalar has_missing_obs
    
    n = rows(Y)
    dim = T - 1
    T_pre = t0 - 1
    has_missing_obs = 0
    
    // ========================================================================
    // Validate data completeness: each period-group cell must have observations
    // ========================================================================
    for (t = 1; t <= T; t++) {
        n_t_1 = sum((time :== t) :& (D :== 1))
        n_t_0 = sum((time :== t) :& (D :== 0))
        if (n_t_1 == 0 | n_t_0 == 0) {
            has_missing_obs = 1
            errprintf("Period %f has no observations in %s group\n", 
                      t, (n_t_1 == 0 ? "treatment" : "control"))
        }
    }
    
    if (has_missing_obs) {
        errprintf("Cannot compute valid covariance matrix due to empty cells.\n")
        return(J(n, dim, .))
    }
    
    // Initialize influence matrix (all zeros)
    psi = J(n, dim, 0)
    
    col_idx = 1
    
    // ========================================================================
    // Part 1: Influence functions for iterative violations nu_t (t = 2,...,t0-1)
    // ========================================================================
    //
    // Population estimand (Section 2.1):
    //   nu_t = E[Y_t - Y_{t-1} | D=1] - E[Y_t - Y_{t-1} | D=0]
    //
    // Sample analog:
    //   nu_hat_t = (Y_bar_{t,1} - Y_bar_{t-1,1}) - (Y_bar_{t,0} - Y_bar_{t-1,0})
    //
    // Influence function (by linearity):
    //   psi_i(nu_t) = (n/n_{t,1})   * 1{D=1,t} * (Y_i - Y_bar_{t,1})
    //               - (n/n_{t-1,1}) * 1{D=1,t-1} * (Y_i - Y_bar_{t-1,1})
    //               - (n/n_{t,0})   * 1{D=0,t} * (Y_i - Y_bar_{t,0})
    //               + (n/n_{t-1,0}) * 1{D=0,t-1} * (Y_i - Y_bar_{t-1,0})
    
    for (t = 2; t <= t0 - 1; t++) {
        // Get group means and counts
        mean_t_1 = _pretest_group_mean(Y, D, time, t, 1)
        mean_t_0 = _pretest_group_mean(Y, D, time, t, 0)
        mean_tm1_1 = _pretest_group_mean(Y, D, time, t - 1, 1)
        mean_tm1_0 = _pretest_group_mean(Y, D, time, t - 1, 0)
        
        n_t_1 = sum((time :== t) :& (D :== 1))
        n_t_0 = sum((time :== t) :& (D :== 0))
        n_tm1_1 = sum((time :== (t-1)) :& (D :== 1))
        n_tm1_0 = sum((time :== (t-1)) :& (D :== 0))
        
        // Accumulate influence contributions for each observation
        for (i = 1; i <= n; i++) {
            if (D[i] == 1) {
                // Treatment group contributions
                if (time[i] == t) {
                    psi[i, col_idx] = psi[i, col_idx] + (Y[i] - mean_t_1) * n / n_t_1
                }
                if (time[i] == t - 1) {
                    psi[i, col_idx] = psi[i, col_idx] - (Y[i] - mean_tm1_1) * n / n_tm1_1
                }
            }
            else {
                // Control group contributions (opposite signs)
                if (time[i] == t) {
                    psi[i, col_idx] = psi[i, col_idx] - (Y[i] - mean_t_0) * n / n_t_0
                }
                if (time[i] == t - 1) {
                    psi[i, col_idx] = psi[i, col_idx] + (Y[i] - mean_tm1_0) * n / n_tm1_0
                }
            }
        }
        
        col_idx++
    }
    
    // ========================================================================
    // Part 2: Influence functions for DID estimands delta_t (t = t0,...,T)
    // ========================================================================
    //
    // Population estimand (Section 2.1):
    //   delta_t = E[Y_t - Y_{t0} | D=1] - E[Y_t - Y_{t0} | D=0]
    //
    // Sample analog:
    //   delta_hat_t = (Y_bar_{t,1} - Y_bar_{t0,1}) - (Y_bar_{t,0} - Y_bar_{t0,0})
    //
    // Influence function (by linearity):
    //   psi_i(delta_t) = (n/n_{t,1})  * 1{D=1,t} * (Y_i - Y_bar_{t,1})
    //                  - (n/n_{t0,1}) * 1{D=1,t0} * (Y_i - Y_bar_{t0,1})
    //                  - (n/n_{t,0})  * 1{D=0,t} * (Y_i - Y_bar_{t,0})
    //                  + (n/n_{t0,0}) * 1{D=0,t0} * (Y_i - Y_bar_{t0,0})
    
    // Get t0 reference means and counts
    mean_t0_1 = _pretest_group_mean(Y, D, time, t0, 1)
    mean_t0_0 = _pretest_group_mean(Y, D, time, t0, 0)
    n_t0_1 = sum((time :== t0) :& (D :== 1))
    n_t0_0 = sum((time :== t0) :& (D :== 0))
    
    for (t = t0; t <= T; t++) {
        // Get group means and counts
        mean_t_1 = _pretest_group_mean(Y, D, time, t, 1)
        mean_t_0 = _pretest_group_mean(Y, D, time, t, 0)
        n_t_1 = sum((time :== t) :& (D :== 1))
        n_t_0 = sum((time :== t) :& (D :== 0))
        
        // Accumulate influence contributions for each observation
        for (i = 1; i <= n; i++) {
            if (D[i] == 1) {
                // Treatment group contributions
                if (time[i] == t) {
                    psi[i, col_idx] = psi[i, col_idx] + (Y[i] - mean_t_1) * n / n_t_1
                }
                if (time[i] == t0) {
                    psi[i, col_idx] = psi[i, col_idx] - (Y[i] - mean_t0_1) * n / n_t0_1
                }
            }
            else {
                // Control group contributions (opposite signs)
                if (time[i] == t) {
                    psi[i, col_idx] = psi[i, col_idx] - (Y[i] - mean_t_0) * n / n_t_0
                }
                if (time[i] == t0) {
                    psi[i, col_idx] = psi[i, col_idx] + (Y[i] - mean_t0_0) * n / n_t0_0
                }
            }
        }
        
        col_idx++
    }
    
    return(psi)
}


// ============================================================================
// CLUSTER-ROBUST COVARIANCE ESTIMATION
// ============================================================================

/**
 * @function _pretest_covariance_cluster
 * @brief Cluster-robust covariance matrix estimator
 *
 * Implements the cluster-robust variance estimator for settings where
 * observations are correlated within clusters (e.g., panel data with
 * individual-level clustering, or geographic clustering).
 *
 * Formula (Cameron & Miller, 2015):
 *   Sigma_hat = (G/(G-1)) * (1/n) * sum_{g=1}^G u_g * u_g'
 *
 * where:
 *   - G is the number of clusters
 *   - u_g = sum_{i in cluster g} psi_i is the cluster-aggregated influence
 *   - G/(G-1) is the finite-sample degrees of freedom adjustment
 *   - 1/n normalizes for the influence function scaling
 *
 * Consistency requires G -> infinity; small G leads to downward bias.
 *
 * @param psi_mat Influence function matrix (n x dim)
 * @param cluster Cluster identifier vector (n x 1)
 *
 * @return Sigma_hat: Cluster-robust covariance matrix (dim x dim)
 *
 * @note For balanced clusters with no within-cluster correlation,
 *       this estimator converges to the standard (iid) estimator
 *
 * @see Cameron, A. C. and D. L. Miller. 2015. A Practitioner's Guide to
 *      Cluster-Robust Inference. Journal of Human Resources 50(2): 317-372.
 */
real matrix _pretest_covariance_cluster(real matrix psi_mat, real colvector cluster)
{
    real matrix Sigma, cluster_sums
    real colvector cluster_ids
    real scalar G, n, dim, g, i
    real rowvector u_g
    
    n = rows(psi_mat)
    dim = cols(psi_mat)
    
    // Get unique cluster IDs
    cluster_ids = uniqrows(cluster)
    G = rows(cluster_ids)
    
    // Step 1: Aggregate influence functions by cluster
    // u_g = sum_{i in cluster g} psi_i
    cluster_sums = J(G, dim, 0)
    
    for (i = 1; i <= n; i++) {
        for (g = 1; g <= G; g++) {
            if (cluster[i] == cluster_ids[g]) {
                cluster_sums[g, .] = cluster_sums[g, .] + psi_mat[i, .]
                break
            }
        }
    }
    
    // Step 2: Compute outer product sum
    // sum_{g=1}^G u_g * u_g'
    Sigma = J(dim, dim, 0)
    for (g = 1; g <= G; g++) {
        u_g = cluster_sums[g, .]
        Sigma = Sigma + u_g' * u_g
    }
    
    // Step 3: Apply scaling factors
    // (G/(G-1)): finite-sample degrees of freedom correction
    // (1/n): normalization for influence function scaling
    Sigma = (G / (G - 1)) * Sigma / n   
    return(Sigma)
}

// ============================================================================
// COVARIANCE TRANSFORMATION FOR OVERALL VIOLATIONS MODE
// ============================================================================

/**
 * @function _pretest_transform_sigma
 * @brief Transform covariance from iterative to overall violations mode
 *
 * Applies the linear transformation Sigma^Delta = A * Sigma * A' to convert
 * the covariance matrix from the iterative violations parameterization to
 * the overall (cumulative) violations parameterization.
 *
 * Mathematical Background (Appendix C):
 *   The overall violations are cumulative sums of iterative violations:
 *     nubar_t = sum_{s=2}^t nu_s  for t = 2, ..., t0-1
 *
 *   This is a linear transformation: theta^Delta = A * theta, where
 *   the transformation matrix A has block structure:
 *
 *     A = | L           0        |
 *         | 0       I_{T_post}   |
 *
 *   where L is a lower triangular matrix of 1's implementing cumulative sum,
 *   and I_{T_post} is the identity (DID estimands unchanged).
 *
 *   By the delta method: Sigma^Delta = A * Sigma * A'
 *
 * @param Sigma  Covariance matrix in iterative mode (T-1 x T-1)
 * @param T_pre  Number of pre-treatment periods (t0 - 1)
 * @param T_post Number of post-treatment periods (T - t0 + 1)
 *
 * @return Sigma^Delta: Covariance matrix in overall mode (T-1 x T-1)
 *
 * @see Mikhaeil & Harshaw (2025), Appendix C, Lemma on variance transformation
 */
real matrix _pretest_transform_sigma(real matrix Sigma, real scalar T_pre,
                                      real scalar T_post)
{
    real matrix A, Sigma_delta
    real scalar dim, i, j, n_pre
    
    n_pre = T_pre - 1  // Number of violation parameters
    dim = n_pre + T_post
    
    // Validate input dimensions
    if (rows(Sigma) != dim || cols(Sigma) != dim) {
        errprintf("Sigma dimension mismatch: expected %f x %f, got %f x %f\n", 
                  dim, dim, rows(Sigma), cols(Sigma))
        _error(3200)
    }
    
    // Construct transformation matrix A
    // Start with identity, then modify upper-left block
    A = I(dim)
    
    // Build lower triangular block L for cumulative sum
    // L[i,j] = 1 for j <= i, implementing: nubar_t = sum_{s=2}^t nu_s
    for (i = 1; i <= n_pre; i++) {
        for (j = 1; j < i; j++) {
            A[i, j] = 1
        }
    }
    
    // Apply delta method transformation
    Sigma_delta = A * Sigma * A'
    
    return(Sigma_delta)
}


// ============================================================================
// POSITIVE DEFINITENESS ENFORCEMENT
// ============================================================================

/**
 * @function _pretest_ensure_pd
 * @brief Ensure positive definiteness of the covariance matrix
 *
 * Checks whether the covariance matrix is positive definite by examining
 * its eigenvalues. If any eigenvalue is non-positive (within numerical
 * tolerance), applies Tikhonov regularization to restore positive definiteness.
 *
 * This regularization is necessary because:
 *   1. delta_{t0} = 0 by construction (Y_{t0} - Y_{t0}), causing rank deficiency
 *   2. Numerical errors may produce slightly negative eigenvalues
 *   3. Monte Carlo sampling in f_alpha computation requires Cholesky decomposition
 *
 * Regularization method:
 *   Sigma_reg = Sigma + epsilon * I, where epsilon = 1e-10
 *
 * @param Sigma Covariance matrix (symmetric, should be positive semi-definite)
 *
 * @return Sigma or Sigma_reg: Positive definite covariance matrix
 *
 * @note Regularization magnitude is minimal to preserve statistical properties
 * @note Errors if regularization fails to achieve positive definiteness
 */
real matrix _pretest_ensure_pd(real matrix Sigma)
{
    real matrix Sigma_reg
    real colvector eigenvalues
    real scalar dim, epsilon, min_eig
    
    dim = rows(Sigma)
    epsilon = 1e-10
    
    // Check positive definiteness via eigenvalue decomposition
    eigenvalues = symeigenvalues(Sigma)
    min_eig = min(eigenvalues)
    
    // Matrix is positive definite if all eigenvalues exceed tolerance
    if (min_eig > epsilon) {
        return(Sigma)
    }
    
    // Apply Tikhonov regularization: Sigma_reg = Sigma + epsilon * I
    Sigma_reg = Sigma + epsilon * I(dim)
    
    // Verify regularization succeeded
    eigenvalues = symeigenvalues(Sigma_reg)
    min_eig = min(eigenvalues)
    
    if (min_eig <= 0) {
        errprintf("Covariance regularization failed (min eigenvalue: %e)\n", min_eig)
        _error(3300)
    }
    
    return(Sigma_reg)
}


// ============================================================================
// AUXILIARY FUNCTIONS
// ============================================================================

/**
 * @function _pretest_get_n_eff
 * @brief Compute effective sample size for variance scaling
 *
 * Calculates an effective sample size that accounts for the data structure.
 * This is used for variance normalization and diagnostic purposes.
 *
 * @param Y        Outcome vector (n x 1)
 * @param D        Treatment indicator (n x 1)
 * @param time     Time period vector (n x 1)
 * @param T        Total number of periods
 * @param is_panel 1 = panel data, 0 = repeated cross-section
 *
 * @return Effective sample size:
 *         - Panel: n/T (observations per period, assuming balance)
 *         - RCS: harmonic mean of cell sizes (accounts for imbalance)
 */
real scalar _pretest_get_n_eff(real colvector Y, real colvector D,
                                 real colvector time, real scalar T,
                                 real scalar is_panel)
{
    real scalar n, t, d, n_td, sum_inv
    
    if (is_panel) {
        // Panel data: observations per time period
        // Assume balanced panel: total obs / T
        n = floor(rows(Y) / T)
    }
    else {
        // Repeated cross-section: harmonic mean
        // n_eff = (1/(2T) * sum_{t,d} 1/n_{td})^{-1}
        sum_inv = 0
        for (t = 1; t <= T; t++) {
            for (d = 0; d <= 1; d++) {
                n_td = sum((time :== t) :& (D :== d))
                if (n_td > 0) {
                    sum_inv = sum_inv + 1 / n_td
                }
            }
        }
        n = (sum_inv / (2 * T))^(-1)
    }
    
    return(n)
}


/**
 * @function _pretest_extract_sigma_submatrix
 * @brief Extract principal submatrix from covariance matrix
 *
 * Extracts a square submatrix corresponding to a subset of the theta
 * parameters. Useful for computing variances of estimands that depend
 * only on a subset of parameters (e.g., pre-treatment violations only).
 *
 * @param Sigma     Full covariance matrix (dim x dim)
 * @param row_start Starting index (1-indexed, inclusive)
 * @param row_end   Ending index (1-indexed, inclusive)
 *
 * @return Submatrix: Sigma[row_start..row_end, row_start..row_end]
 */
real matrix _pretest_extract_sigma_submatrix(real matrix Sigma,
                                               real scalar row_start,
                                               real scalar row_end)
{
    if (row_start < 1 || row_end > rows(Sigma) || row_start > row_end) {
        errprintf("Error: Invalid submatrix indices\n")
        _error(3300)
    }
    
    return(Sigma[row_start::row_end, row_start::row_end])
}


end
