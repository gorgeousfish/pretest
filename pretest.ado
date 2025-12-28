*! pretest - Conditional Extrapolation Pre-Test for Difference-in-Differences
*! Version 0.1.0
*!
*! Reference:
*!   Mikhaeil, J.M. and C. Harshaw. 2025. "In Defense of the Pre-Test: Valid
*!   Inference when Testing Violations of Parallel Trends for Difference-in-
*!   Differences." arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470
*!
*! Description:
*!   Implements the conditional extrapolation framework (Assumption 3) for DID
*!   designs. This package provides:
*!   (1) A preliminary test to determine whether extrapolation of parallel
*!       trend violations from pre- to post-treatment is justified;
*!   (2) Conditionally valid confidence intervals for the ATT when the test
*!       passes (Theorem 2).
*!
*! Syntax:
*!   pretest depvar, treatment(varname) time(varname) threshold(real) [options]
*!
*! Required arguments:
*!   depvar       - Outcome variable (numeric)
*!   treatment()  - Binary treatment indicator (0=control, 1=treated)
*!   time()       - Time period variable (numeric)
*!   threshold()  - Acceptable violation threshold M > 0 (Section 3.1)
*!
*! Options:
*!   treat_time(#)    - Treatment time t0 (first post-treatment period)
*!   p(#)             - Severity norm exponent p >= 1 (default: 2)
*!   alpha(#)         - Significance level alpha in (0,1) (default: 0.05)
*!   level(#)         - Confidence level (default: 95)
*!   cluster(varname) - Cluster variable for robust standard errors
*!   overall          - Use overall violations (Appendix C); default: iterative
*!   nograph          - Suppress event study plot
*!   simulate(#)      - Monte Carlo draws for critical value (default: 5000)
*!   seed(#)          - Random number seed (default: 12345)
*!
*! Stored results:
*!   e(b)            - Coefficient vector (ATT estimate)
*!   e(V)            - Variance matrix
*!   e(T)            - Total number of time periods T
*!   e(t0)           - Treatment time t0
*!   e(T_pre)        - Pre-treatment periods T_pre = t0 - 1
*!   e(T_post)       - Post-treatment periods T_post = T - t0 + 1
*!   e(S_pre)        - Estimated severity S_pre (Section 4.2)
*!   e(phi)          - Pre-test indicator phi = 1{S_pre > M} (Theorem 1)
*!   e(kappa)        - Constant kappa (Proposition 1)
*!   e(f_alpha)      - Critical value f(alpha, Sigma) (Section 5.1)
*!   e(delta_bar)    - Average DID estimate delta_bar
*!   e(ci_lower)     - CI lower bound (Theorem 2)
*!   e(ci_upper)     - CI upper bound (Theorem 2)
*!   e(pretest_pass) - 1 if pre-test passed (phi = 0), 0 otherwise
*!
*! Example:
*!   . sysuse auto, clear
*!   . pretest price, treatment(foreign) time(rep78) threshold(500) treat_time(4)

program define pretest, eclass
    version 17.0
    
    // =========================================================================
    // STEP 0: Parse command syntax
    // =========================================================================
    syntax varlist(min=1 max=1 numeric) [if] [in], ///
        TREATment(varname numeric) ///
        TIME(varname numeric) ///
        THREshold(real) ///
        [TREat_time(real -999) ///
         P(real 2) ///
         ALpha(real -999) ///
         Level(cilevel) ///
         CLuster(varname) ///
         OVERall ///
         NOGraph ///
         SIMulate(integer 5000) ///
         SEED(integer 12345) ///
         DIAGnose ///
         *]
    
    local outcome `varlist'
    
    // Create sample marker for estimates store compatibility
    marksample touse
    
    // =========================================================================
    // STEP 0.1: Validate input parameters
    // =========================================================================
    
    // Threshold M must be strictly positive (Section 3.1, Assumption 3)
    if `threshold' <= 0 {
        di as error "threshold() must be strictly positive"
        exit 198
    }
    
    // Severity norm p must satisfy p >= 1 (Section 3.1)
    if `p' < 1 {
        di as error "p() must be >= 1"
        exit 198
    }
    
    // Reconcile alpha and level options
    if `alpha' == -999 {
        local alpha = 1 - `level'/100
    }
    else {
        if `alpha' <= 0 | `alpha' >= 1 {
            di as error "alpha() must be in the open interval (0, 1)"
            exit 198
        }
        local level = 100 * (1 - `alpha')
    }
    
    // Monte Carlo simulations must be sufficient for stable critical values
    if `simulate' < 100 {
        di as error "simulate() must be >= 100 for reliable critical value estimation"
        exit 198
    }
    
    // Set violation mode: iterative (default) or overall (Appendix C)
    if "`overall'" != "" {
        local mode "overall"
    }
    else {
        local mode "iterative"
    }
    
    local sims = `simulate'
    
    // =========================================================================
    // STEP 0.2: Load Mata library (if not already in memory)
    // =========================================================================
    // Priority: 1) Pre-compiled .mlib (fast, recommended)
    //           2) Source compilation via _pretest_mata.do (fallback)
    
    capture mata: mata describe _pretest_main_compute()
    local mata_exists = (_rc == 0)
    
    if !`mata_exists' {
        local mata_loaded = 0
        
        // Method 1: Load pre-compiled Mata library (lpretest.mlib)
        // This is the preferred method - fast and reliable
        capture findfile lpretest.mlib
        if _rc == 0 {
            // Library found, update mlib index to load it
            quietly mata: mata mlib index
            capture mata: mata describe _pretest_main_compute()
            if _rc == 0 {
                local mata_loaded = 1
            }
        }
        
        // Method 2: Fallback to source compilation if .mlib not found
        if `mata_loaded' == 0 {
            // Get directory containing ado file
            local ado_dir ""
            local base_dir ""
            capture findfile pretest.ado
            if _rc == 0 {
                local ado_dir = subinstr("`r(fn)'", "pretest.ado", "", .)
                // Check if ado file is in /ado/ subdirectory - if so, base_dir is parent
                if strpos("`ado_dir'", "/ado/") > 0 {
                    local base_dir = subinstr("`ado_dir'", "/ado/", "/", .)
                }
                else {
                    local base_dir "`ado_dir'"
                }
            }
            
            // Try base_dir (parent of /ado/)
            if "`base_dir'" != "" & `mata_loaded' == 0 {
                capture quietly do "`base_dir'_pretest_mata.do"
                if _rc == 0 {
                    local mata_loaded = 1
                }
            }
            
            // Try ado_dir directly
            if "`ado_dir'" != "" & `mata_loaded' == 0 {
                capture quietly do "`ado_dir'_pretest_mata.do"
                if _rc == 0 {
                    local mata_loaded = 1
                }
            }
            
            // Final fallback - use findfile
            if `mata_loaded' == 0 {
                capture findfile _pretest_mata.do
                if _rc == 0 {
                    capture quietly do "`r(fn)'"
                }
            }
        }
    }
    
    // =========================================================================
    // STEP 1: Validate input data
    // =========================================================================
    _pretest_validate `outcome', treatment(`treatment') time(`time')
    local n_total_obs = r(n_total)
    
    // Verify complete cell coverage: the parameter vector
    // theta = (nu_2,...,nu_{t0-1}, delta_{t0},...,delta_T) requires observations
    // in all (time, treatment) cells to satisfy Assumption 1.
    // Note: Coverage check removed - validation handled by Mata computation
    
    // =========================================================================
    // STEP 2: Detect panel structure (panel vs. repeated cross-sections)
    // =========================================================================
    _pretest_detect_panel
    local is_panel = r(is_panel)
    if `is_panel' {
        local panelvar = r(panelvar)
        local timevar_xt = r(timevar)
        di as text "Panel data detected: panelvar=`panelvar', timevar=`timevar_xt'"
    }
    else {
        di as text "Repeated cross-section data detected"
    }
    
    // =========================================================================
    // STEP 3: Parse time structure
    // =========================================================================
    if `treat_time' == -999 {
        _pretest_parse_time `time', treatment(`treatment') is_panel(`is_panel')
    }
    else {
        _pretest_parse_time `time', treatment(`treatment') ///
            treat_time(`treat_time') is_panel(`is_panel')
    }
    
    local T = r(T)
    local t0 = r(t0)
    local T_pre = r(T_pre)
    local T_post = r(T_post)
    
    // Display time structure
    di as text _n "Time Structure:"
    di as text "  T (total periods)     = " as result `T'
    di as text "  t0 (treatment time)   = " as result `t0'
    di as text "  T_pre (pre-treatment) = " as result `T_pre'
    di as text "  T_post (post-treatment) = " as result `T_post'
    
    // =========================================================================
    // STEP 4: Compute DID estimates and conduct pre-test
    // =========================================================================
    // Reference: Mikhaeil & Harshaw (2025), Section 2
    //   delta_t = E[Y_t - Y_{t0} | D=1] - E[Y_t - Y_{t0} | D=0]
    //   delta_bar = (1/T_post) * sum_{t=t0}^T delta_t
    
    // Get time mapping (handles gaps in time variable)
    tempname time_vals_mat
    matrix `time_vals_mat' = r(time_vals)
    
    // Pass parameters to Mata via Stata scalars and matrices
    scalar __pretest_T = `T'
    scalar __pretest_t0 = `t0'
    scalar __pretest_T_pre = `T_pre'
    scalar __pretest_T_post = `T_post'
    scalar __pretest_p = `p'
    scalar __pretest_alpha = `alpha'
    scalar __pretest_threshold = `threshold'
    scalar __pretest_sims = `sims'
    scalar __pretest_seed = `seed'
    scalar __pretest_is_panel = `is_panel'
    if "`mode'" == "overall" {
        scalar __pretest_mode_overall = 1
    }
    else {
        scalar __pretest_mode_overall = 0
    }
    
    // Pass time values matrix
    matrix __pretest_time_vals = `time_vals_mat'
    
    // Call main Mata computation function
    mata: _pretest_main_compute("`outcome'", "`treatment'", "`time'", "`cluster'")
    
    // =========================================================================
    // STEP 5: Retrieve results from Mata
    // =========================================================================
    // Core quantities: S_pre (severity), kappa, f_alpha (critical value),
    // phi (pre-test indicator), delta_bar (ATT estimate), CI bounds
    local delta_bar_result = scalar(__pretest_delta_bar)
    local S_pre_result = scalar(__pretest_S_pre)
    local kappa_result = scalar(__pretest_kappa)
    local f_alpha_result = scalar(__pretest_f_alpha)
    local phi_result = scalar(__pretest_phi)
    local n_result = scalar(__pretest_n)
    local ci_lo = scalar(__pretest_ci_lower)
    local ci_hi = scalar(__pretest_ci_upper)
    
    // Prepare CI arguments for display (handle missing values)
    if missing(`ci_lo') {
        local ci_lo_arg = -999
    }
    else {
        local ci_lo_arg = `ci_lo'
    }
    if missing(`ci_hi') {
        local ci_hi_arg = -999
    }
    else {
        local ci_hi_arg = `ci_hi'
    }
    
    // Prepare phi argument (handle missing values)
    // Definition: phi = 1{S_pre > M}, values: 0 (pass), 1 (fail), or . (data issue)
    if missing(`phi_result') {
        // If phi missing, recompute from S_pre and threshold
        if missing(`S_pre_result') {
            local phi_int = .  // Data issue
        }
        else if `S_pre_result' > `threshold' {
            local phi_int = 1
        }
        else {
            local phi_int = 0
        }
    }
    else {
        // Ensure integer conversion (if not missing)
        if !missing(`phi_result') {
            local phi_int = int(`phi_result')
        }
        else {
            local phi_int = .  // Data issue
        }
    }
    
    // =========================================================================
    // STEP 6: Store estimation results in e()
    // =========================================================================
    // Following Stata convention: e(b) = coefficient vector, e(V) = variance matrix
    // This enables estimates store, esttab, and other post-estimation commands.
    tempname b_mat V_mat delta_mat nu_mat Sigma_mat Sigma_full
    
    // Get full covariance matrix
    capture matrix `Sigma_full' = __pretest_Sigma
    local has_sigma = (_rc == 0)
    
    // Build b vector (1 x 1): ATT = delta_bar
    matrix `b_mat' = J(1, 1, `delta_bar_result')
    matrix colnames `b_mat' = ATT
    matrix rownames `b_mat' = y1
    
    // Build V matrix (1 x 1): Var(ATT) = Var(delta_bar)
    // Sigma is (T-1) x (T-1) covariance matrix of theta = (nu, delta)
    local se_delta_bar_val = scalar(__pretest_se_delta_bar)
    local var_att = 0
    if !missing(`se_delta_bar_val') & `se_delta_bar_val' > 0 {
        local var_att = `se_delta_bar_val' * `se_delta_bar_val'
    }
    else if `has_sigma' {
        // Fallback: compute from Sigma matrix
        local dim = rowsof(`Sigma_full')
        local expected_dim = `T' - 1  // Sigma is (T-1) x (T-1)
        
        // Delta portion starts at T_pre (1-indexed)
        // because ν has T_pre-1 elements
        if `dim' == `expected_dim' & `T_post' > 0 {
            local delta_start = `T_pre'
            local delta_end = `T_pre' + `T_post' - 1
            
            // Extract delta portion of covariance matrix
            tempname Sigma_delta
            capture matrix `Sigma_delta' = `Sigma_full'[`delta_start'..`delta_end', `delta_start'..`delta_end']
            
            if _rc == 0 {
                // Compute ATT variance: Var(ATT) = (1/T_post^2) * 1' * Sigma_delta * 1
                local var_att = 0
                local delta_dim = rowsof(`Sigma_delta')
                forvalues i = 1/`delta_dim' {
                    forvalues j = 1/`delta_dim' {
                        local elem = `Sigma_delta'[`i', `j']
                        if !missing(`elem') {
                            local var_att = `var_att' + `elem'
                        }
                    }
                }
                local var_att = `var_att' / (`T_post' * `T_post')
            }
        }
    }
    
    // Ensure variance is non-negative and non-missing
    if `var_att' <= 0 | missing(`var_att') {
        // Use small positive value as default variance
        local var_att = 1e-10
    }
    
    matrix `V_mat' = J(1, 1, `var_att')
    matrix colnames `V_mat' = ATT
    matrix rownames `V_mat' = ATT
    
    // --- Post results via ereturn (supports estimates store) ---
    // Report valid observations (non-missing Y, D, time) for correct sample size
    ereturn post `b_mat' `V_mat', esample(`touse') obs(`n_result')
    
    // --- Time structure ---
    ereturn scalar T = `T'
    ereturn scalar t0 = `t0'
    ereturn scalar T_pre = `T_pre'
    ereturn scalar T_post = `T_post'
    
    // --- Model parameters ---
    ereturn scalar M = `threshold'           // Threshold (preferred name)
    ereturn scalar threshold = `threshold'   // Threshold (alias)
    ereturn scalar p = `p'                   // Norm parameter
    ereturn scalar alpha = `alpha'           // Significance level
    ereturn scalar level = `level'           // Confidence level
    
    // --- Sample information ---
    ereturn scalar n = `n_result'            // Valid sample size (after time mapping)
    ereturn scalar N = `n_result'            // Valid sample size (Stata convention alias)
    ereturn scalar n_total = `n_total_obs'   // Total rows (for reference)
    ereturn scalar is_panel = `is_panel'     // Panel data indicator
    
    // --- Pre-test results (Theorem 1, Proposition 1) ---
    ereturn scalar S_pre = `S_pre_result'    // Severity S_pre
    ereturn scalar kappa = `kappa_result'    // kappa constant
    ereturn scalar f_alpha = `f_alpha_result' // Critical value f(alpha, Sigma)
    ereturn scalar phi = `phi_int'           // phi = 1{S_pre > M}
    
    // --- Pre-test pass indicator ---
    // pretest_pass = 1 iff phi = 0, i.e., S_pre <= M (extrapolation justified)
    if `phi_int' == 0 {
        ereturn scalar pretest_pass = 1      // Pass
    }
    else {
        ereturn scalar pretest_pass = 0      // Fail (phi=1: violation too large, phi=.: data issue)
    }
    
    // Check for data validity issues
    if missing(`phi_int') | missing(`S_pre_result') | missing(`f_alpha_result') {
        ereturn scalar data_valid = 0    // Data invalid (missing periods or computation failed)
        di as error "Warning: Some periods have no observations in treatment or control group."
        di as error "         Cannot compute valid confidence interval. pretest_pass set to 0."
    }
    else {
        ereturn scalar data_valid = 1    // Data valid
    }
    
    // --- ATT estimate ---
    ereturn scalar delta_bar = `delta_bar_result' // delta_bar
    ereturn scalar ATT = `delta_bar_result'       // ATT (alias)
    
    // --- Conditional confidence interval (Theorem 2) ---
    if !missing(`ci_lo') & !missing(`ci_hi') {
        ereturn scalar ci_lower = `ci_lo'    // CI lower bound
        ereturn scalar ci_upper = `ci_hi'    // CI upper bound
    }
    
    // --- Conventional DID confidence interval (for comparison) ---
    local ci_conv_lo = scalar(__pretest_ci_conv_lower)
    local ci_conv_hi = scalar(__pretest_ci_conv_upper)
    local se_delta_bar_val = scalar(__pretest_se_delta_bar)
    
    if !missing(`ci_conv_lo') & !missing(`ci_conv_hi') {
        ereturn scalar ci_conv_lower = `ci_conv_lo'    // Conventional CI lower
        ereturn scalar ci_conv_upper = `ci_conv_hi'    // Conventional CI upper
    }
    if !missing(`se_delta_bar_val') {
        ereturn scalar se_delta_bar = `se_delta_bar_val'  // ATT standard error
    }
    
    // --- Severity standard error (Delta method) ---
    local S_pre_se_val = scalar(__pretest_S_pre_se)
    if !missing(`S_pre_se_val') {
        ereturn scalar S_pre_se = `S_pre_se_val'   // SE(S_pre) via Delta method
    }
    
    // --- Monte Carlo parameters ---
    ereturn scalar sims = `simulate'         // Monte Carlo simulations
    ereturn scalar seed = `seed'             // Random seed
    
    // --- Matrices ---
    // delta vector: individual DID estimates (T_post x 1)
    capture matrix `delta_mat' = __pretest_delta
    if _rc == 0 {
        matrix rownames `delta_mat' = delta
        matrix colnames `delta_mat' = value
        ereturn matrix delta = `delta_mat'
    }
    
    // nu vector: iterative violations (T_pre-1 x 1)
    capture matrix `nu_mat' = __pretest_nu
    if _rc == 0 {
        ereturn matrix nu = `nu_mat'
    }
    
    // theta vector: full parameter vector ((T-1) x 1)
    // theta = (nu_2, ..., nu_{t0-1}, delta_{t0}, ..., delta_T) per Section 2.2
    tempname theta_mat
    capture matrix `theta_mat' = __pretest_theta
    if _rc == 0 {
        ereturn matrix theta = `theta_mat'
    }
    
    // Sigma: covariance matrix of theta
    if `has_sigma' {
        ereturn matrix Sigma = `Sigma_full'
    }
    
    // Clean up global scalars and matrices
    capture scalar drop __pretest_*
    capture matrix drop __pretest_*
    
    // --- Macros ---
    ereturn local depvar "`outcome'"         // Dependent variable (Stata convention)
    ereturn local outcome "`outcome'"        // Dependent variable (alias)
    ereturn local treatment "`treatment'"    // Treatment variable
    ereturn local time "`time'"              // Time variable
    ereturn local mode "`mode'"              // Mode (iterative/overall)
    if "`cluster'" != "" {
        ereturn local clustvar "`cluster'"   // Cluster variable (Stata convention)
        ereturn local cluster "`cluster'"    // Cluster variable (alias)
    }
    ereturn local cmd "pretest"              // Command name
    ereturn local cmdline "pretest `0'"      // Full command line
    ereturn local title "Conditional Extrapolation Pre-Test"
    ereturn local predict ""                 // Predict command (not supported)
    
    // =========================================================================
    // STEP 7: Display results
    // =========================================================================
    // Note: Display must come after ereturn to preserve estimation results
    
    // Prepare conventional CI arguments
    if missing(`ci_conv_lo') {
        local ci_conv_lo_arg = -999
    }
    else {
        local ci_conv_lo_arg = `ci_conv_lo'
    }
    if missing(`ci_conv_hi') {
        local ci_conv_hi_arg = -999
    }
    else {
        local ci_conv_hi_arg = `ci_conv_hi'
    }
    
    _pretest_output, ///
        threshold(`threshold') ///
        mode("`mode'") ///
        alpha(`alpha') ///
        level(`level') ///
        t(`T') ///
        t0(`t0') ///
        t_pre(`T_pre') ///
        t_post(`T_post') ///
        n(`n_result') ///
        p(`p') ///
        s_pre(`S_pre_result') ///
        kappa(`kappa_result') ///
        f_alpha(`f_alpha_result') ///
        phi(`phi_int') ///
        delta_bar(`delta_bar_result') ///
        ci_lower(`ci_lo_arg') ///
        ci_upper(`ci_hi_arg') ///
        ci_conv_lower(`ci_conv_lo_arg') ///
        ci_conv_upper(`ci_conv_hi_arg') ///
        sims(`sims') ///
        seed(`seed')
    
    // =========================================================================
    // STEP 8: Reasonability checks and warnings
    // =========================================================================
    local show_warnings = 0
    
    // Warning: T_post = 1 yields degenerate conventional CI
    // (delta_{t0} = 0 by definition when t = t0)
    if `T_post' == 1 {
        if `show_warnings' == 0 {
            di as text ""
            di as text "{bf:Warnings:}"
            local show_warnings = 1
        }
        di as text "  {it:- T_post=1: Only one post-treatment period (t=t0).}"
        di as text "  {it:  By definition delta_{t0}=0, so Conventional CI is degenerate.}"
        di as text "  {it:  Recommend T_post >= 2 for meaningful ATT estimation.}"
    }
    
    // Note: Very small S_pre suggests near-perfect parallel trends
    if !missing(`S_pre_result') & `S_pre_result' < 0.001 & `S_pre_result' >= 0 {
        if `show_warnings' == 0 {
            di as text ""
            di as text "{bf:Warnings:}"
            local show_warnings = 1
        }
        di as text "  {it:- Very small S_pre suggests near-perfect parallel trends.}"
    }
    
    // =========================================================================
    // STEP 9: Diagnostic output (optional)
    // =========================================================================
    if "`diagnose'" != "" {
        di as text ""
        di as text "{hline 72}"
        di as text "{bf:Diagnostic Information}"
        di as text "{hline 72}"
        
        // theta vector components
        di as text ""
        di as text "{bf:Parameter Vector theta (T-1 components):}"
        tempname theta_mat
        capture matrix `theta_mat' = e(theta)
        if _rc == 0 {
            matrix list `theta_mat', noheader format(%9.4f)
        }
        
        // Pre-treatment violations nu_t
        di as text ""
        di as text "{bf:Pre-treatment Violations nu_t (t=2,...,t0-1):}"
        tempname nu_mat
        capture matrix `nu_mat' = e(nu)
        if _rc == 0 {
            matrix list `nu_mat', noheader format(%9.4f)
        }
        
        // Post-treatment DID estimates delta_t
        di as text ""
        di as text "{bf:Post-treatment DID Estimates delta_t (t=t0,...,T):}"
        tempname delta_mat
        capture matrix `delta_mat' = e(delta)
        if _rc == 0 {
            matrix list `delta_mat', noheader format(%9.4f)
        }
        
        // Sigma matrix eigenvalues (for positive definiteness check)
        di as text ""
        di as text "{bf:Covariance Matrix Sigma Eigenvalues:}"
        tempname Sigma_mat eigenvals
        capture matrix `Sigma_mat' = e(Sigma)
        if _rc == 0 {
            capture matrix symeigen eigenvecs `eigenvals' = `Sigma_mat'
            if _rc == 0 {
                matrix list `eigenvals', noheader format(%9.4e)
            }
        }
        
        // Monte Carlo simulation parameters
        di as text ""
        di as text "{bf:Monte Carlo Simulation:}"
        di as text "  Simulations: " e(sims)
        di as text "  Seed: " e(seed)
        di as text "  f(alpha, Sigma): " %9.4f e(f_alpha)
        
        di as text "{hline 72}"
    }
    
    // =========================================================================
    // STEP 10: Generate event study plot (unless suppressed)
    // =========================================================================
    if "`nograph'" == "" {
        capture _pretest_graph, `options'
        if _rc != 0 {
            di as text "(Note: Graph generation failed with error `=_rc'. Results still valid.)"
        }
    }
end
