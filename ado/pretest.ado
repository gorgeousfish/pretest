*! pretest v0.1.0 - Conditional Extrapolation Pre-Test for DID
*! Based on Mikhaeil & Harshaw (2025), "In Defense of the Pre-Test"
*! arXiv:2510.26470
*! 
*! Description:
*!   Implements the conditional extrapolation pre-test framework for 
*!   difference-in-differences designs. Provides a preliminary test to 
*!   determine whether extrapolation from pre-treatment parallel trend 
*!   violations to post-treatment is justified, along with valid confidence 
*!   intervals conditional on passing the test.
*!
*! Syntax:
*!   pretest outcome, treatment(varname) time(varname) threshold(#) [options]
*!
*! Required:
*!   outcome      - Outcome variable (numeric)
*!   treatment()  - Treatment indicator (0/1 binary)
*!   time()       - Time variable
*!   threshold()  - Acceptable violation threshold M > 0
*!
*! Options:
*!   treat_time(#)    - Treatment time point (first post-treatment period)
*!   p(#)             - Severity norm p >= 1 (default: 2)
*!   alpha(#)         - Significance level (default: 0.05)
*!   level(#)         - Confidence level (default: 95, alternative to alpha)
*!   cluster(varname) - Cluster variable for standard errors
*!   overall          - Use overall violations mode (default: iterative)
*!   nograph          - Suppress event study graph
*!   simulate(#)      - Monte Carlo simulations (default: 5000)
*!   seed(#)          - Random seed (default: 12345)
*!
*! Returns:
*!   e(T)            - Total number of time periods
*!   e(t0)           - Treatment time point
*!   e(T_pre)        - Number of pre-treatment periods
*!   e(T_post)       - Number of post-treatment periods
*!   e(S_pre)        - Estimated pre-treatment severity
*!   e(phi)          - Pre-test result (1=reject extrapolation, 0=accept)
*!   e(delta_bar)    - Average post-treatment DID estimate
*!   e(kappa)        - κ constant
*!   e(f_alpha)      - Critical value
*!   e(ci_lower)     - Lower bound of confidence interval
*!   e(ci_upper)     - Upper bound of confidence interval
*!   e(pretest_pass) - 1 if pretest passed (0 if failed)
*!
*! Example:
*!   . webuse nlswork, clear
*!   . xtset idcode year
*!   . pretest ln_wage, treatment(union) time(year) treat_time(80)
*!
*! References:
*!   Mikhaeil, J.M. and Harshaw, C. (2025). "In Defense of the Pre-Test: 
*!   Valid Inference when Testing Violations of Parallel Trends for 
*!   Difference-in-Differences." arXiv:2510.26470.
*!
*! Authors:
*!   Stata implementation: Xuanyu Cai, Wenli Xu

program define pretest, eclass
    version 17.0
    
    // ========================================
    // Step 0: Parse syntax
    // ========================================
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
         DIAGnose]
    
    local outcome `varlist'
    
    // Create sample marker for estimates store compatibility
    marksample touse
    
    // ========================================
    // Step 0.1: Parameter Validation
    // ========================================
    
    // Validate threshold > 0 (error code 105)
    if `threshold' <= 0 {
        di as error "Error 105: threshold() must be positive"
        di as error "  threshold = `threshold' is invalid"
        exit 105
    }
    
    // Validate p >= 1 (error code 106)
    if `p' < 1 {
        di as error "Error 106: p() must be >= 1"
        di as error "  p = `p' is invalid (use p >= 1, default p = 2)"
        exit 106
    }
    
    // Handle alpha/level compatibility
    if `alpha' == -999 {
        // User did not specify alpha, compute from level
        local alpha = 1 - `level'/100
    }
    else {
        // User specified alpha, validate range (error code 107)
        if `alpha' <= 0 | `alpha' >= 1 {
            di as error "Error 107: alpha() must be in (0, 1)"
            di as error "  alpha = `alpha' is invalid"
            exit 107
        }
    }
    
    // Validate simulate >= 100 (error code 110)
    if `simulate' < 100 {
        di as error "Error 110: simulate() must be >= 100"
        di as error "  simulate = `simulate' is too small for reliable critical values"
        exit 110
    }
    
    // Set mode (iterative vs overall)
    if "`overall'" != "" {
        local mode "overall"
    }
    else {
        local mode "iterative"
    }
    
    // Store simulate value for downstream code
    local sims = `simulate'
    
    // Display header delegated to _pretest_output
    
    // Step 0: Load Mata functions
    // Get directory containing ado file
    local ado_dir ""
    capture findfile pretest.ado
    if _rc == 0 {
        local ado_dir = subinstr("`r(fn)'", "pretest.ado", "", .)
    }
    
    // Attempt to load Mata functions
    local mata_loaded = 0
    if "`ado_dir'" != "" {
        capture quietly do "`ado_dir'_pretest_mata.do"
        if _rc == 0 {
            local mata_loaded = 1
        }
    }
    
    // Final fallback: use findfile
    if `mata_loaded' == 0 {
        capture findfile _pretest_mata.do
        if _rc == 0 {
            capture quietly do "`r(fn)'"
        }
    }
    
    // Step 1: Validate input data
    _pretest_validate `outcome', treatment(`treatment') time(`time')
    
    // Step 2: Detect panel structure
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
    
    // Step 3: Parse time structure
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
    
    // ========================================
    // Step 4: Compute DID estimates
    // ========================================
    // Reference: Mikhaeil & Harshaw (2025), Section 2.1
    // δ̂_t = (ĀŸ_{t,D=1} - ĀŸ_{t0,D=1}) - (ĀŸ_{t,D=0} - ĀŸ_{t0,D=0})
    // δ̄̂ = (1/T_post) × Σ_{t=t0}^T δ̂_t
    
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
    
    // ========================================
    // Retrieve Results from Mata
    // ========================================
    
    // Retrieve results from global scalars
    local delta_bar_result = scalar(__pretest_delta_bar)
    local S_pre_result = scalar(__pretest_S_pre)
    local kappa_result = scalar(__pretest_kappa)
    local f_alpha_result = scalar(__pretest_f_alpha)
    local phi_result = scalar(__pretest_phi)
    local n_result = scalar(__pretest_n)
    local ci_lo = scalar(__pretest_ci_lower)
    local ci_hi = scalar(__pretest_ci_upper)
    
    // ========================================
    // Display Results
    // ========================================
    
    // Prepare CI arguments (handle missing values)
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
    // Definition: φ = 1{Ŝ_pre > M}, values: 0, 1, or . (data issue)
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
    
    // ========================================
    // Store Results in e()
    // ========================================
    // Use ereturn post to publish b and V matrices for estimates store
    // b = [ATT] (1x1 coefficient vector)
    // V = [Var(ATT)] (1x1 variance matrix)
    
    // --- Prepare b and V matrices ---
    tempname b_mat V_mat delta_mat nu_mat Sigma_mat Sigma_full
    
    // Get full covariance matrix
    capture matrix `Sigma_full' = __pretest_Sigma
    local has_sigma = (_rc == 0)
    
    // Build b vector (1 x 1): ATT
    matrix `b_mat' = J(1, 1, `delta_bar_result')
    matrix colnames `b_mat' = ATT
    matrix rownames `b_mat' = y1
    
    // Build V matrix (1 x 1): Var(ATT) = SE(δ̄)²
    // Use se_delta_bar from Mata (most accurate variance estimate)
    // Note: Sigma is (T-1) x (T-1) matrix where:
    //   - First T_pre-1 elements correspond to ν (t=2 to t=t0-1)
    //   - Last T_post elements correspond to δ (t=t0 to t=T)
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
    ereturn post `b_mat' `V_mat', esample(`touse') obs(`n_result')
    
    // --- Scalars: Time structure ---
    ereturn scalar T = `T'
    ereturn scalar t0 = `t0'
    ereturn scalar T_pre = `T_pre'
    ereturn scalar T_post = `T_post'
    
    // --- Scalars: Model parameters ---
    ereturn scalar M = `threshold'           // Threshold (preferred name)
    ereturn scalar threshold = `threshold'   // Threshold (alias)
    ereturn scalar p = `p'                   // Norm parameter
    ereturn scalar alpha = `alpha'           // Significance level
    ereturn scalar level = `level'           // Confidence level
    
    // --- Scalars: Sample information (N set by ereturn post) ---
    ereturn scalar n = `n_result'            // Sample size (lowercase alias)
    ereturn scalar is_panel = `is_panel'     // Panel data indicator
    
    // --- Scalars: Pre-test core results ---
    ereturn scalar S_pre = `S_pre_result'    // Severity Ŝ_pre
    ereturn scalar kappa = `kappa_result'    // κ constant
    ereturn scalar f_alpha = `f_alpha_result' // Critical value f(α, Σ̂)
    ereturn scalar phi = `phi_int'           // φ = 1{Ŝ_pre > M}
    
    // --- Scalars: Pre-test result (key) ---
    // pretest_pass = (S_pre <= M) = (phi == 0)
    // Note: phi=0 means pass, phi=1 means fail, phi=. means data issue
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
    
    // --- Scalars: ATT estimates ---
    ereturn scalar delta_bar = `delta_bar_result' // δ̄̂
    ereturn scalar ATT = `delta_bar_result'       // ATT (alias)
    
    // --- Scalars: Confidence interval ---
    if !missing(`ci_lo') & !missing(`ci_hi') {
        ereturn scalar ci_lower = `ci_lo'    // CI lower bound
        ereturn scalar ci_upper = `ci_hi'    // CI upper bound
    }
    
    // --- Scalars: Conventional confidence interval ---
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
    
    // --- Scalars: Severity standard error ---
    local S_pre_se_val = scalar(__pretest_S_pre_se)
    if !missing(`S_pre_se_val') {
        ereturn scalar S_pre_se = `S_pre_se_val'   // SE(Ŝ_pre) via Delta Method
    }
    
    // --- Scalars: Simulation parameters ---
    ereturn scalar sims = `simulate'         // Monte Carlo simulations
    ereturn scalar seed = `seed'             // Random seed
    
    // --- Matrices ---
    // δ vector (T_post x 1)
    capture matrix `delta_mat' = __pretest_delta
    if _rc == 0 {
        matrix rownames `delta_mat' = delta
        matrix colnames `delta_mat' = value
        ereturn matrix delta = `delta_mat'
    }
    
    // ν vector (T_pre x 1) - iterative violations
    capture matrix `nu_mat' = __pretest_nu
    if _rc == 0 {
        ereturn matrix nu = `nu_mat'
    }
    
    // θ̂ vector ((T-1) x 1) - full parameter vector
    // Reference: Section 2.2, θ̂ = (ν̂_2, ..., ν̂_{t0-1}, δ̂_{t0}, ..., δ̂_T)
    tempname theta_mat
    capture matrix `theta_mat' = __pretest_theta
    if _rc == 0 {
        ereturn matrix theta = `theta_mat'
    }
    
    // Σ̂ covariance matrix
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
    
    // ========================================
    // Display Results
    // ========================================
    // Note: Output must come after ereturn to avoid clearing estimation results
    
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
    
    // Conventional CI is now displayed by _pretest_output via ci_conv_lower/ci_conv_upper parameters
    // (Format number helper is defined within _pretest_output.ado)
    
    // ========================================
    // Result Reasonability Checks and Warnings
    // ========================================
    local show_warnings = 0
    
    // Check 0: T_post=1 boundary condition
    // When T_post=1, delta_bar = delta_{t0} = 0 by definition (Y_{t0} - Y_{t0} = 0)
    // This renders the conventional CI degenerate (se_delta_bar → 0)
    if `T_post' == 1 {
        if `show_warnings' == 0 {
            di as text ""
            di as text "{bf:Warnings:}"
            local show_warnings = 1
        }
        di as text "  {it:- T_post=1: Only one post-treatment period (t=t₀).}"
        di as text "  {it:  By definition δ̂_{t₀}=0, so Conventional CI is degenerate.}"
        di as text "  {it:  Recommend T_post ≥ 2 for meaningful ATT estimation.}"
    }
    
    // Check 1: Small sample size (n < 100)
    if `n_result' < 100 {
        if `show_warnings' == 0 {
            di as text ""
            di as text "{bf:Warnings:}"
            local show_warnings = 1
        }
        di as text "  {it:- Small sample size (n=`n_result' < 100). Results may be imprecise.}"
    }
    
    // Check 2: Very wide CI (width > 10)
    if !missing(`ci_lo') & !missing(`ci_hi') {
        local ci_width = `ci_hi' - `ci_lo'
        if `ci_width' > 10 {
            if `show_warnings' == 0 {
                di as text ""
                di as text "{bf:Warnings:}"
                local show_warnings = 1
            }
            local ci_width_fmt : di %9.1f `ci_width'
            local ci_width_fmt = strtrim("`ci_width_fmt'")
            di as text "  {it:- Wide confidence interval (width=`ci_width_fmt'). Consider data quality.}"
        }
    }
    
    // Check 3: Very small S_pre (near-perfect parallel trends)
    if !missing(`S_pre_result') & `S_pre_result' < 0.001 & `S_pre_result' >= 0 {
        if `show_warnings' == 0 {
            di as text ""
            di as text "{bf:Warnings:}"
            local show_warnings = 1
        }
        di as text "  {it:- Very small Ŝ_pre suggests near-perfect parallel trends.}"
    }
    
    // ========================================
    // Diagnostic Output
    // ========================================
    if "`diagnose'" != "" {
        di as text ""
        di as text "{hline 72}"
        di as text "{bf:Diagnostic Information}"
        di as text "{hline 72}"
        
        // θ̂ vector components
        di as text ""
        di as text "{bf:Parameter Vector θ̂ (T-1 components):}"
        tempname theta_mat
        capture matrix `theta_mat' = e(theta)
        if _rc == 0 {
            matrix list `theta_mat', noheader format(%9.4f)
        }
        
        // Pre-treatment violations ν_t
        di as text ""
        di as text "{bf:Pre-treatment Violations ν̂_t (t=2,...,t0-1):}"
        tempname nu_mat
        capture matrix `nu_mat' = e(nu)
        if _rc == 0 {
            matrix list `nu_mat', noheader format(%9.4f)
        }
        
        // Post-treatment DID estimates δ_t
        di as text ""
        di as text "{bf:Post-treatment DID Estimates δ̂_t (t=t0,...,T):}"
        tempname delta_mat
        capture matrix `delta_mat' = e(delta)
        if _rc == 0 {
            matrix list `delta_mat', noheader format(%9.4f)
        }
        
        // Σ̂ matrix eigenvalues
        di as text ""
        di as text "{bf:Covariance Matrix Σ̂ Eigenvalues:}"
        tempname Sigma_mat eigenvals
        capture matrix `Sigma_mat' = e(Sigma)
        if _rc == 0 {
            capture matrix symeigen eigenvecs `eigenvals' = `Sigma_mat'
            if _rc == 0 {
                matrix list `eigenvals', noheader format(%9.4e)
            }
        }
        
        // Monte Carlo convergence diagnostics
        di as text ""
        di as text "{bf:Monte Carlo Simulation:}"
        di as text "  Simulations: " e(sims)
        di as text "  Seed: " e(seed)
        di as text "  f(α,Σ̂): " %9.4f e(f_alpha)
        
        di as text "{hline 72}"
    }
    
    // ========================================
    // Generate Event Study Graph
    // ========================================
    // Only call graph if nograph option not specified
    if "`nograph'" == "" {
        capture _pretest_graph
        if _rc != 0 {
            di as text "(Note: Graph generation failed with error `=_rc'. Results still valid.)"
        }
    }
end
