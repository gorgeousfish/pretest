*! _pretest_covariance.mata v0.1.0
*! Asymptotic Covariance Matrix Estimation
*!
*! Description:
*!   Implements covariance estimation for the parameter vector theta using
*!   influence function methods. Supports both standard and cluster-robust
*!   variance estimation.
*!
*! Theory:
*!   Under Assumption 1: sqrt(n)(theta_hat - theta) -> N(0, Sigma)
*!   where theta = (nu_2, ..., nu_{t0-1}, delta_{t0}, ..., delta_T)'
*!
*! Contents:
*!   - Influence function computation for DID and violation estimators
*!   - Standard covariance estimation
*!   - Cluster-robust covariance estimation
*!   - Covariance transformation for overall mode
*!   - Positive definiteness enforcement
*!
*! Reference:
*!   Mikhaeil, J.M. and C. Harshaw. 2025. In Defense of the Pre-Test: Valid
*!   Inference when Testing Violations of Parallel Trends for Difference-in-
*!   Differences. arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470

version 17.0

mata:

// ============================================================================
// CORE COVARIANCE ESTIMATION
// ============================================================================

/**
 * @function _pretest_covariance
 * @brief Main asymptotic covariance matrix estimator
 *
 * Estimates Sigma such that sqrt(n)(theta_hat - theta) -> N(0, Sigma).
 * Uses influence function method for sample mean estimators.
 *
 * @param Y        Outcome vector (n x 1)
 * @param D        Treatment indicator (n x 1)
 * @param time     Time period vector (n x 1)
 * @param t0       Treatment time
 * @param T        Total number of periods
 * @param is_panel 1 = panel data, 0 = repeated cross-section
 * @param cluster  Optional cluster variable for robust SEs
 *
 * @return Asymptotic covariance matrix Sigma_hat of dimension (T-1) x (T-1)
 *
 * @note Automatically regularizes if not positive definite
 * @note Supports cluster-robust standard errors
 *
 * @see Mikhaeil & Harshaw (2025), Section 2.2-2.3
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
 * @brief Compute influence function matrix for theta estimator
 *
 * Constructs the n x (T-1) matrix of influence functions where each row
 * is observation i's contribution to the estimator theta_hat.
 *
 * For violation nu_t:
 *   psi_i = (Y_i - mean)*n/n_td terms for each group-time combination
 *
 * For DID delta_t:
 *   Similar linear combination referencing time t0
 *
 * @param Y, D, time, t0, T Same as _pretest_covariance()
 *
 * @return Influence function matrix of dimension n x (T-1)
 *
 * @note Column sums should be approximately zero (centering property)
 *
 * @see Standard asymptotic theory for sample mean estimators
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
    
    // ========================================
    // Check if all periods have observations
    // ========================================
    for (t = 1; t <= T; t++) {
        n_t_1 = sum((time :== t) :& (D :== 1))
        n_t_0 = sum((time :== t) :& (D :== 0))
        if (n_t_1 == 0 | n_t_0 == 0) {
            has_missing_obs = 1
            printf("{err}Error: Period %f has no observations in %s group\n", t, (n_t_1 == 0 ? "treatment" : "control"))
        }
    }
    
    // If missing observations, return matrix with missing values
    // This will trigger f_alpha to be missing later
    if (has_missing_obs) {
        printf("{err}       Cannot compute valid covariance matrix.\n")
        printf("{err}       Setting pretest_pass to 0, confidence interval to missing.\n")
        return(J(n, dim, .))
    }
    
    // Initialize influence matrix (all zeros)
    psi = J(n, dim, 0)
    
    col_idx = 1
    
    // ========================================
    // Part 1: Influence for violations nu_t (t = 2, ..., t0-1)
    // ========================================
    // Formula: nu_t = (mean_Y_t_1 - mean_Y_{t-1}_1) - (mean_Y_t_0 - mean_Y_{t-1}_0)
    // 
    // Influence function derivation:
    //   For mean_Y_t_d = (1/n_{td}) * sum_i Y_i * 1{D_i=d, time_i=t}
    //   Influence function: d(mean_Y_t_d)/dY_i = 1{D_i=d, time_i=t} / n_{td}
    //   
    //   Therefore, influence function for nu_t:
    //   psi_i(nu_t) = 1{D_i=1, t_i=t}/n_t1 - 1{D_i=1, t_i=t-1}/n_{t-1,1}
    //               - 1{D_i=0, t_i=t}/n_t0 + 1{D_i=0, t_i=t-1}/n_{t-1,0}
    
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
        
        // Compute influence for each observation
        // Correct influence function formula: psi_i = (Y_i - mean) * n / n_td
        // Reference: Influence function for sqrt(n)(mean - mu) in asymptotic theory
        for (i = 1; i <= n; i++) {
            if (D[i] == 1) {
                // Treatment group
                if (time[i] == t) {
                    psi[i, col_idx] = psi[i, col_idx] + (Y[i] - mean_t_1) * n / n_t_1
                }
                if (time[i] == t - 1) {
                    psi[i, col_idx] = psi[i, col_idx] - (Y[i] - mean_tm1_1) * n / n_tm1_1
                }
            }
            else {
                // Control group
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
    
    // ========================================
    // Part 2: Influence for DID delta_t (t = t0, ..., T)
    // ========================================
    // Formula: delta_t = (mean_Y_t_1 - mean_Y_t0_1) - (mean_Y_t_0 - mean_Y_t0_0)
    // 
    // Influence function:
    //   psi_i(delta_t) = 1{D_i=1, t_i=t}/n_t1 - 1{D_i=1, t_i=t0}/n_{t0,1}
    //                  - 1{D_i=0, t_i=t}/n_t0 + 1{D_i=0, t_i=t0}/n_{t0,0}
    
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
        
        // Compute influence for each observation
        // Correct influence function formula: psi_i = (Y_i - mean) * n / n_td
        for (i = 1; i <= n; i++) {
            if (D[i] == 1) {
                // Treatment group
                if (time[i] == t) {
                    psi[i, col_idx] = psi[i, col_idx] + (Y[i] - mean_t_1) * n / n_t_1
                }
                if (time[i] == t0) {
                    psi[i, col_idx] = psi[i, col_idx] - (Y[i] - mean_t0_1) * n / n_t0_1
                }
            }
            else {
                // Control group
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
// CLUSTER-ROBUST COVARIANCE
// ============================================================================

/**
 * @function _pretest_covariance_cluster
 * @brief Cluster-robust covariance estimator
 *
 * Implements cluster-robust standard errors following Cameron & Miller (2015):
 *
 *   Sigma = (G/(G-1)) * (1/n) * sum_g u_g * u_g'
 *
 * where u_g = sum_{i in cluster g} psi_i.
 *
 * @param psi_mat Influence function matrix (n x dim)
 * @param cluster Cluster variable (n x 1)
 *
 * @return Cluster-robust covariance matrix (dim x dim)
 *
 * @note G/(G-1) is the degrees of freedom adjustment
 * @note Division by n corrects for influence function scaling
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
    
    // Aggregate influence functions by cluster
    cluster_sums = J(G, dim, 0)
    
    for (i = 1; i <= n; i++) {
        for (g = 1; g <= G; g++) {
            if (cluster[i] == cluster_ids[g]) {
                cluster_sums[g, .] = cluster_sums[g, .] + psi_mat[i, .]
                break
            }
        }
    }
    
    // Compute cluster-robust covariance
    // Formula: Sigma = (G/(G-1)) * (1/n) * sum_g u_g' * u_g
    //
    // This matches the non-clustered formula which uses (1/(n-1)) * psi'*psi
    // Both should give similar results when there's no within-cluster correlation
    //
    // Note: The (G/(G-1)) factor is the standard finite-sample degrees of freedom adjustment
    //       following Cameron & Miller (2015)
    Sigma = J(dim, dim, 0)
    for (g = 1; g <= G; g++) {
        u_g = cluster_sums[g, .]
        Sigma = Sigma + u_g' * u_g
    }
    
    // Apply scaling: (G/(G-1)) for df adjustment, (1/n) for correct magnitude
    Sigma = (G / (G - 1)) * Sigma / n   
    return(Sigma)
}

// ============================================================================
// COVARIANCE TRANSFORMATION (for Overall Mode)
// ============================================================================

/**
 * @function _pretest_transform_sigma
 * @brief Transform covariance matrix for overall violations mode
 *
 * Applies linear transformation Sigma^Delta = A * Sigma * A' where A
 * implements cumulative summation for violations (Appendix C).
 *
 * Transformation structure:
 *   A = | L           0        |
 *       | 0       I_{T_post}   |
 *
 * where L is lower triangular matrix of 1's (cumulative sum operator).
 *
 * @param Sigma  Iterative mode covariance matrix
 * @param T_pre  Number of pre-treatment periods
 * @param T_post Number of post-treatment periods
 *
 * @return Overall mode covariance matrix Sigma^Delta
 *
 * @see Mikhaeil & Harshaw (2025), Appendix C
 */
real matrix _pretest_transform_sigma(real matrix Sigma, real scalar T_pre,
                                      real scalar T_post)
{
    real matrix A, Sigma_delta
    real scalar dim, i, j, n_pre
    
    n_pre = T_pre - 1  // Number of violation parameters
    dim = n_pre + T_post
    
    // Verify dimension
    if (rows(Sigma) != dim || cols(Sigma) != dim) {
        errprintf("Error: Sigma dimension mismatch in transformation\n")
        errprintf("  Expected: %f x %f, Got: %f x %f\n", dim, dim, rows(Sigma), cols(Sigma))
        _error(3200)
    }
    
    // Initialize transformation matrix A as identity
    A = I(dim)
    
    // Modify upper-left block to lower triangular of 1's
    // This implements cumulative sum: nubar_t = sum_{s=2}^t nu_s
    for (i = 1; i <= n_pre; i++) {
        for (j = 1; j < i; j++) {
            A[i, j] = 1
        }
    }
    
    // Compute transformation: Sigma^Delta = A * Sigma * A'
    Sigma_delta = A * Sigma * A'
    
    return(Sigma_delta)
}


// ============================================================================
// POSITIVE DEFINITENESS ENFORCEMENT
// ============================================================================

/**
 * @function _pretest_ensure_pd
 * @brief Ensure positive definiteness of covariance matrix
 *
 * Checks eigenvalues and applies diagonal perturbation if needed to
 * ensure positive definiteness for Cholesky decomposition.
 *
 * @param Sigma Covariance matrix
 *
 * @return Regularized covariance matrix (if needed)
 *
 * @note Regularization: Sigma_reg = Sigma + epsilon * I (epsilon = 1e-10)
 * @note Outputs informational message if regularization applied
 */
real matrix _pretest_ensure_pd(real matrix Sigma)
{
    real matrix Sigma_reg
    real colvector eigenvalues
    real scalar dim, epsilon, min_eig
    
    dim = rows(Sigma)
    epsilon = 1e-10
    
    // Test positive definiteness via eigenvalues
    eigenvalues = symeigenvalues(Sigma)
    min_eig = min(eigenvalues)
    
    // Matrix is PD if all eigenvalues > 0
    // We use epsilon as threshold to account for numerical errors
    if (min_eig > epsilon) {
        // Matrix is already positive definite
        return(Sigma)
    }
    
    // Matrix is not positive definite; apply regularization
    // Note: This is expected when delta_{t0} = 0 (Y_{t0} - Y_{t0} = 0 by definition)
    // Regularization ensures numerical stability for Monte Carlo sampling
    // Users can inspect eigenvalues via the diagnose option
    
    // Add epsilon to diagonal
    Sigma_reg = Sigma + epsilon * I(dim)
    
    // Verify regularization worked
    eigenvalues = symeigenvalues(Sigma_reg)
    min_eig = min(eigenvalues)
    
    if (min_eig <= 0) {
        errprintf("Error: Covariance matrix regularization failed\n")
        errprintf("  Min eigenvalue after regularization: %e\n", min_eig)
        _error(3300)
    }
    
    return(Sigma_reg)
}


// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

/**
 * @function _pretest_get_n_eff
 * @brief Calculate effective sample size
 *
 * @param Y        Outcome vector
 * @param D        Treatment indicator
 * @param time     Time period vector
 * @param T        Total number of periods
 * @param is_panel 1 if panel, 0 if repeated cross-section
 *
 * @return Effective sample size:
 *         - Panel: observations per period (total/T)
 *         - RCS: harmonic mean across group-time cells
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
 * @brief Extract submatrix for subset of theta parameters
 *
 * @param Sigma     Full covariance matrix
 * @param row_start Starting row/col index
 * @param row_end   Ending row/col index
 *
 * @return Extracted submatrix Sigma[row_start..row_end, row_start..row_end]
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




