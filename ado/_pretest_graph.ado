*! _pretest_graph.ado v0.1.0
*! Event Study Visualization for Conditional Pre-Test Inference
*!
*! Description:
*!   Generates event study plot displaying pre-treatment parallel trend 
*!   violations and post-treatment DID estimates with confidence intervals.
*!   Implements the graphical representation described in Section 1.1 of 
*!   Mikhaeil & Harshaw (2025).
*!
*! Syntax:
*!   _pretest_graph [, TItle(string) SAVing(string) NAME(string) REPLACE]
*!
*! Options:
*!   title(string) - Custom graph title (default: includes author citation)
*!   saving(string)- Export graph to file (.gph, .png, .pdf, .eps)
*!   name(string)  - Stata graph window name (for managing multiple graphs)
*!   replace       - Overwrite existing file when saving
*!
*! Required e() Returns (from pretest command):
*!   e(nu)           - Iterative violations nu_t, t=2,...,t0-1 (T_pre-1 x 1)
*!   e(delta)        - DID estimates delta_t, t=t0,...,T (T_post x 1)
*!   e(Sigma)        - Asymptotic covariance matrix of theta (T-1 x T-1)
*!   e(threshold)    - Acceptable violation threshold M (Section 3.1)
*!   e(t0)           - Treatment time (first post-treatment period)
*!   e(T_pre)        - Number of pre-treatment periods
*!   e(T_post)       - Number of post-treatment periods
*!   e(N)            - Sample size n
*!   e(alpha)        - Significance level (default: 0.05)
*!   e(mode)         - "iterative" (nu_t) or "overall" (nubar_t) mode
*!   e(pretest_pass) - Preliminary test result (1=pass, 0=fail)
*!   e(delta_bar)    - Average post-treatment ATT estimator
*!   e(ci_lower)     - Conditional CI lower bound (Theorem 2)
*!   e(ci_upper)     - Conditional CI upper bound (Theorem 2)
*!   e(ci_conv_lower)- Conventional CI lower bound (for comparison)
*!   e(ci_conv_upper)- Conventional CI upper bound (for comparison)
*!
*! Reference:
*!   Mikhaeil, J. M. and C. Harshaw. 2025. In Defense of the Pre-Test: Valid
*!   Inference when Testing Violations of Parallel Trends for Difference-in-
*!   Differences. arXiv preprint arXiv:2510.26470.
*!   https://arxiv.org/abs/2510.26470

program define _pretest_graph
    version 17.0
    
    // Parse options (follows Stata post-estimation command conventions)
    // Element-specific options allow granular control over visual elements
    // User-provided twoway_options (*) are passed through to the final graph
    syntax [, TItle(string) SAVing(string) NAME(string) REPLACE ///
        CI_opt_pass(string asis) ///
        CI_opt_fail(string asis) ///
        LINE_opt_m(string asis) ///
        MARKER_opt_pre(string asis) ///
        MARKER_opt_post(string asis) ///
        *]
    
    // =========================================================================
    // STEP 1: Validate e() Returns
    // =========================================================================
    // This command requires prior execution of -pretest- which stores
    // estimation results in e(). Exit with error 301 if results not found.
    
    // Validate required matrices: nu (violations), delta (DID), Sigma (covariance)
    capture confirm matrix e(nu)
    if _rc != 0 {
        di as error "e(nu) not found. Run {bf:pretest} command first."
        exit 301
    }
    
    capture confirm matrix e(delta)
    if _rc != 0 {
        di as error "e(delta) not found. Run {bf:pretest} command first."
        exit 301
    }
    
    capture confirm matrix e(Sigma)
    if _rc != 0 {
        di as error "e(Sigma) not found. Run {bf:pretest} command first."
        exit 301
    }
    
    // Validate required scalars for inference
    local required_scalars "t0 T_pre T_post N alpha threshold pretest_pass"
    foreach s of local required_scalars {
        capture confirm scalar e(`s')
        if _rc != 0 {
            di as error "e(`s') not found. Run {bf:pretest} command first."
            exit 301
        }
    }
    
    // =========================================================================
    // STEP 2: Extract Estimation Results from e()
    // =========================================================================
    
    // Time structure parameters
    local t0 = e(t0)              // Treatment time (first post-treatment period)
    local T_pre = e(T_pre)        // Number of pre-treatment periods
    local T_post = e(T_post)      // Number of post-treatment periods
    local n = e(N)                // Sample size
    
    // Inference parameters
    local alpha = e(alpha)        // Significance level (typically 0.05)
    local threshold = e(threshold) // Acceptable violation threshold M (Section 3.1)
    local pretest_pass = e(pretest_pass)  // Preliminary test result (Theorem 1)
    local delta_bar = e(delta_bar)        // Average DID: delta_bar (Section 2.1)
    
    // Violation mode: "iterative" uses nu_t, "overall" uses nubar_t (Appendix C)
    local mode = e(mode)
    if "`mode'" == "" {
        local mode "iterative"
    }
    
    // Confidence interval bounds (missing if pretest failed)
    // Conditional CI from Theorem 2; conventional CI for comparison
    local ci_lower = e(ci_lower)
    local ci_upper = e(ci_upper)
    local ci_conv_lower = e(ci_conv_lower)
    local ci_conv_upper = e(ci_conv_upper)
    
    // Parameter vector components (Section 2.2, Equation 205):
    // theta = (nu_2, ..., nu_{t0-1}, delta_{t0}, ..., delta_T)'
    tempname nu_mat delta_mat Sigma_mat
    matrix `nu_mat' = e(nu)       // Iterative violations (T_pre-1 x 1)
    matrix `delta_mat' = e(delta) // DID estimates (T_post x 1)
    matrix `Sigma_mat' = e(Sigma) // Asymptotic covariance (T-1 x T-1)
    
    // Dimension verification
    local T_pre_minus_1 = rowsof(`nu_mat')  // Should equal T_pre - 1
    local T_post_actual = rowsof(`delta_mat')
    
    // Critical value for Wald-type confidence intervals
    local z_crit = invnormal(1 - `alpha'/2)
    
    // =========================================================================
    // STEP 2.5: Define Default Styles (Mikhaeil-Harshaw 2025 Figure 1)
    // =========================================================================
    // These defaults implement the original paper's visual design.
    // User options are appended AFTER defaults to enable overrides via Stata's
    // "last option wins" precedence rule.
    
    // Pre-treatment confidence intervals (navy)
    local def_ci_pre "lcolor(navy) lwidth(medthick)"
    
    // Pre-treatment point markers (navy circles)
    local def_marker_pre "mcolor(navy) msymbol(O) msize(medlarge)"
    
    // Post-treatment confidence intervals (styling depends on test result)
    local def_ci_pass "lcolor(maroon) lwidth(medthick)"
    local def_ci_fail "lcolor(maroon%50) lwidth(medthick) lpattern(dash)"
    
    // Post-treatment point markers (maroon)
    local def_marker_post_pass "mcolor(maroon) msymbol(O) msize(medlarge)"
    local def_marker_post_fail "mcolor(maroon%50) msymbol(Oh) msize(medlarge)"
    
    // Threshold M horizontal lines (orange dashed)
    local def_line_m "lcolor(orange%70) lpattern(shortdash)"
    
    // Apply user overrides: append user options to defaults
    // Pre-treatment CI (uses pass styling as base for pre-treatment)
    local ci_pre_opts "`def_ci_pre'"
    if `"`ci_opt_pass'"' != "" {
        local ci_pre_opts `"`ci_pre_opts' `ci_opt_pass'"'
    }
    
    // Pre-treatment markers
    local marker_pre_opts "`def_marker_pre'"
    if `"`marker_opt_pre'"' != "" {
        local marker_pre_opts `"`marker_pre_opts' `marker_opt_pre'"'
    }
    
    // Post-treatment CI (pass case)
    local ci_pass_opts "`def_ci_pass'"
    if `"`ci_opt_pass'"' != "" {
        local ci_pass_opts `"`ci_pass_opts' `ci_opt_pass'"'
    }
    
    // Post-treatment CI (fail case)
    local ci_fail_opts "`def_ci_fail'"
    if `"`ci_opt_fail'"' != "" {
        local ci_fail_opts `"`ci_fail_opts' `ci_opt_fail'"'
    }
    
    // Post-treatment markers
    if `pretest_pass' == 1 {
        local marker_post_opts "`def_marker_post_pass'"
    }
    else {
        local marker_post_opts "`def_marker_post_fail'"
    }
    if `"`marker_opt_post'"' != "" {
        local marker_post_opts `"`marker_post_opts' `marker_opt_post'"'
    }
    
    // Threshold M line
    local line_m_opts "`def_line_m'"
    if `"`line_opt_m'"' != "" {
        local line_m_opts `"`line_m_opts' `line_opt_m'"'
    }
    
    // =========================================================================
    // STEP 3: Construct Plot Data
    // =========================================================================
    // Create dataset with event study estimates and confidence intervals.
    // Time axis uses relative time: negative for pre-treatment, 0+ for post.
    
    preserve
    clear
    
    // Total observations: one per violation + one per DID estimate
    local total_points = `T_pre_minus_1' + `T_post_actual'
    quietly set obs `total_points'
    
    quietly {
        // Initialize plot variables
        gen period = .      // Relative time (treatment at 0)
        gen estimate = .    // Point estimate (nu_t or delta_t)
        gen ci_lo = .       // Lower CI bound
        gen ci_hi = .       // Upper CI bound
        gen is_pre = .      // Indicator: 1=pre-treatment, 0=post-treatment
        
        // -----------------------------------------------------------------
        // Pre-treatment: Parallel trend violations
        // -----------------------------------------------------------------
        // Map absolute time t in {2,...,t0-1} to relative time {-(T_pre-1),...,-1}
        
        if "`mode'" == "overall" {
            // Overall mode: nubar_t = sum_{s=2}^t nu_s (Appendix C, Eq. 700)
            // Requires linear transformation of covariance: Sigma^Delta = A*Sigma*A'
            tempname nu_bar A_mat Sigma_nu Sigma_nu_overall
            
            // Construct cumulative sum operator A (lower triangular of 1s)
            matrix `A_mat' = J(`T_pre_minus_1', `T_pre_minus_1', 0)
            forvalues i = 1/`T_pre_minus_1' {
                forvalues j = 1/`i' {
                    matrix `A_mat'[`i', `j'] = 1
                }
            }
            
            // Transform violations: nubar = A * nu
            matrix `nu_bar' = `A_mat' * `nu_mat'
            
            // Extract and transform covariance for violations block
            matrix `Sigma_nu' = `Sigma_mat'[1..`T_pre_minus_1', 1..`T_pre_minus_1']
            matrix `Sigma_nu_overall' = `A_mat' * `Sigma_nu' * `A_mat''
            
            // Populate pre-treatment observations (overall violations)
            forvalues i = 1/`T_pre_minus_1' {
                local rel_time = `i' - `T_pre'
                replace period = `rel_time' in `i'
                
                local est = `nu_bar'[`i', 1]
                replace estimate = `est' in `i'
                
                // Standard error derivation (Assumption 1, Section 2.3):
                // sqrt(n)(theta_hat - theta) -> N(0, Sigma)
                // => Var(theta_hat) = Sigma/n, SE = sqrt(Sigma_ii/n)
                local var_i = `Sigma_nu_overall'[`i', `i']
                local se_i = sqrt(`var_i' / `n')
                replace ci_lo = `est' - `z_crit' * `se_i' in `i'
                replace ci_hi = `est' + `z_crit' * `se_i' in `i'
                replace is_pre = 1 in `i'
            }
        }
        else {
            // Iterative mode: use nu_t directly (default, Section 2.1)
            forvalues i = 1/`T_pre_minus_1' {
                local rel_time = `i' - `T_pre'
                replace period = `rel_time' in `i'
                
                local est = `nu_mat'[`i', 1]
                replace estimate = `est' in `i'
                
                // Standard error from Sigma diagonal (Assumption 1)
                local var_i = `Sigma_mat'[`i', `i']
                local se_i = sqrt(`var_i' / `n')
                replace ci_lo = `est' - `z_crit' * `se_i' in `i'
                replace ci_hi = `est' + `z_crit' * `se_i' in `i'
                replace is_pre = 1 in `i'
            }
        }
        
        // -----------------------------------------------------------------
        // Post-treatment: DID estimates delta_t
        // -----------------------------------------------------------------
        // Map absolute time t in {t0,...,T} to relative time {0,...,T_post-1}
        
        forvalues i = 1/`T_post_actual' {
            local obs = `T_pre_minus_1' + `i'
            local rel_time = `i' - 1  // t0 maps to 0
            replace period = `rel_time' in `obs'
            
            local est = `delta_mat'[`i', 1]
            replace estimate = `est' in `obs'
            
            // Delta estimates are in theta[T_pre:T-1] (1-indexed: T_pre_minus_1+1)
            local idx = `T_pre_minus_1' + `i'
            local var_i = `Sigma_mat'[`idx', `idx']
            local se_i = sqrt(`var_i' / `n')
            replace ci_lo = `est' - `z_crit' * `se_i' in `obs'
            replace ci_hi = `est' + `z_crit' * `se_i' in `obs'
            replace is_pre = 0 in `obs'
        }
    }
    
    // =========================================================================
    // STEP 4: Construct Event Study Graph
    // =========================================================================
    // Visual design follows Figure 1 of Mikhaeil & Harshaw (2025):
    // - Pre-treatment: violations nu_t or nubar_t with pointwise CIs
    // - Post-treatment: DID estimates delta_t with pointwise CIs
    // - Threshold M shown as horizontal reference lines
    // - Conditional CI (orange) vs conventional CI (gray) for ATT
    
    // Separate pre- and post-treatment series for distinct styling
    tempvar pre_est pre_ci_lo pre_ci_hi post_est post_ci_lo post_ci_hi
    quietly {
        gen `pre_est' = estimate if is_pre == 1
        gen `pre_ci_lo' = ci_lo if is_pre == 1
        gen `pre_ci_hi' = ci_hi if is_pre == 1
        gen `post_est' = estimate if is_pre == 0
        gen `post_ci_lo' = ci_lo if is_pre == 0
        gen `post_ci_hi' = ci_hi if is_pre == 0
    }
    
    // Graph title (default includes author citation)
    if "`title'" == "" {
        local title "Event Study with Pre-test (Mikhaeil-Harshaw 2025)"
    }
    
    // Mode label for legend
    if "`mode'" == "overall" {
        local mode_label "Overall violations"
    }
    else {
        local mode_label "Iterative violations"
    }
    
    // Visual encoding of preliminary test result (Theorem 1)
    if `pretest_pass' == 1 {
        local pretest_label "Pretest: PASS"
        local post_marker "O"
        local post_color "maroon"
    }
    else {
        local pretest_label "Pretest: FAIL"
        local post_marker "Oh"
        local post_color "maroon%50"  // Faded to indicate invalid inference
    }
    
    // X-axis: relative time with integer ticks
    quietly summarize period
    local x_min_data = r(min)
    local x_max_data = r(max)
    local x_min = floor(`x_min_data') - 1
    local x_max = ceil(`x_max_data') + 1
    
    // Y-axis: ensure threshold M is visible
    quietly summarize estimate
    local y_min = min(r(min), -`threshold') - 0.1 * abs(r(max) - r(min))
    local y_max = max(r(max), `threshold') + 0.1 * abs(r(max) - r(min))
    
    // Handle degenerate case
    if `y_min' == `y_max' {
        local y_min = `y_min' - 1
        local y_max = `y_max' + 1
    }
    
    // -----------------------------------------------------------------
    // Build graph layers
    // -----------------------------------------------------------------
    local graph_cmd ""
    
    // Layer 1: Pre-treatment CIs (default: navy, user-overridable)
    local graph_cmd `"`graph_cmd' (rcap `pre_ci_lo' `pre_ci_hi' period, `ci_pre_opts')"'
    
    // Layer 2: Pre-treatment point estimates (default: navy circles, user-overridable)
    local graph_cmd `"`graph_cmd' (scatter `pre_est' period, `marker_pre_opts')"'
    
    // Layer 3: Post-treatment CIs (styling depends on test result, user-overridable)
    if `pretest_pass' == 1 {
        local graph_cmd `"`graph_cmd' (rcap `post_ci_lo' `post_ci_hi' period, `ci_pass_opts')"'
    }
    else {
        // Dashed lines indicate inference may be invalid
        local graph_cmd `"`graph_cmd' (rcap `post_ci_lo' `post_ci_hi' period, `ci_fail_opts')"'
    }
    
    // Layer 4: Post-treatment point estimates (user-overridable)
    local graph_cmd `"`graph_cmd' (scatter `post_est' period, `marker_post_opts')"'
    
    // -----------------------------------------------------------------
    // Average ATT with CI comparison (Figure 1 bottom panel)
    // -----------------------------------------------------------------
    local att_pos = `x_max' - 0.3
    
    tempvar att_period att_est att_ci_lo att_ci_hi
    quietly {
        gen `att_period' = `att_pos' in 1
        gen `att_est' = `delta_bar' in 1
    }
    
    // Conventional CI (gray, always shown for comparison)
    if !missing(`ci_conv_lower') & !missing(`ci_conv_upper') {
        local att_pos_conv = `att_pos' - 0.15
        tempvar att_period_conv att_ci_conv_lo att_ci_conv_hi
        quietly {
            gen `att_period_conv' = `att_pos_conv' in 1
            gen `att_ci_conv_lo' = `ci_conv_lower' in 1
            gen `att_ci_conv_hi' = `ci_conv_upper' in 1
        }
        local graph_cmd `"`graph_cmd' (rcap `att_ci_conv_lo' `att_ci_conv_hi' `att_period_conv', lcolor(gs6) lwidth(medthick))"'
        local graph_cmd `"`graph_cmd' (scatter `att_est' `att_period_conv', mcolor(gs6) msymbol(Oh) msize(medium))"'
    }
    
    // Conditional CI (orange, Theorem 2; only if pretest passed)
    if `pretest_pass' == 1 & !missing(`ci_lower') & !missing(`ci_upper') {
        local att_pos_cond = `att_pos' + 0.15
        tempvar att_period_cond att_ci_cond_lo att_ci_cond_hi
        quietly {
            gen `att_period_cond' = `att_pos_cond' in 1
            gen `att_ci_cond_lo' = `ci_lower' in 1
            gen `att_ci_cond_hi' = `ci_upper' in 1
        }
        local graph_cmd `"`graph_cmd' (rcap `att_ci_cond_lo' `att_ci_cond_hi' `att_period_cond', lcolor(orange) lwidth(thick))"'
        local graph_cmd `"`graph_cmd' (scatter `att_est' `att_period_cond', mcolor(orange) msymbol(D) msize(large))"'
    }
    
    // -----------------------------------------------------------------
    // Annotation and legend
    // -----------------------------------------------------------------
    local note_text "M = `: di %6.3f `threshold''  |  `pretest_label'"
    
    // Legend labels
    if "`mode'" == "overall" {
        local pre_legend "Overall viol."
    }
    else {
        local pre_legend "Iter. viol."
    }
    
    // Build legend dynamically based on available CI types
    local legend_order "2 4"
    local legend_labels `"label(2 "`pre_legend'") label(4 "DID est.")"'
    
    if !missing(`ci_conv_lower') {
        local legend_order "2 4 6"
        local legend_labels `"label(2 "`pre_legend'") label(4 "DID est.") label(6 "Conv. CI")"'
    }
    
    if `pretest_pass' == 1 & !missing(`ci_lower') {
        if !missing(`ci_conv_lower') {
            local legend_order "2 4 6 8"
            local legend_labels `"label(2 "`pre_legend'") label(4 "DID est.") label(6 "Conv. CI") label(8 "Cond. CI")"'
        }
        else {
            local legend_order "2 4 6"
            local legend_labels `"label(2 "`pre_legend'") label(4 "DID est.") label(6 "Cond. CI")"'
        }
    }
    
    // Graph window name option
    local name_opt ""
    if "`name'" != "" {
        local name_opt "name(`name', replace)"
    }
    
    // -----------------------------------------------------------------
    // Render graph
    // -----------------------------------------------------------------
    // Append user-provided twoway_options at the end to enable full customization
    // (Stata's "last option wins" rule ensures user options override defaults)
    twoway `graph_cmd', ///
        xline(-0.5, lcolor(gs10) lpattern(dash) lwidth(medium)) ///
        yline(0, lcolor(gs6) lwidth(thin)) ///
        yline(`threshold', `line_m_opts') ///
        yline(-`threshold', `line_m_opts') ///
        xlabel(`x_min'(1)`x_max', labsize(small)) ///
        ylabel(, angle(horizontal) labsize(small) format(%9.2f)) ///
        xtitle("Time relative to treatment", size(small)) ///
        ytitle("Estimate", size(small)) ///
        title("`title'", size(medium)) ///
        note("`note_text'", size(vsmall)) ///
        legend(order(`legend_order') `legend_labels' pos(6) rows(1) size(vsmall) symxsize(small)) ///
        graphregion(color(white)) plotregion(color(white)) ///
        `name_opt' ///
        `options'
    
    // =========================================================================
    // STEP 5: Export Graph (Optional)
    // =========================================================================
    // Supports Stata native (.gph) and publication formats (.png, .pdf, .eps)
    
    if "`saving'" != "" {
        local ext = substr("`saving'", -4, 4)
        
        if "`ext'" == ".gph" {
            graph save "`saving'", `replace'
        }
        else if "`ext'" == ".png" | "`ext'" == ".pdf" | "`ext'" == ".eps" {
            graph export "`saving'", `replace'
        }
        else {
            // Default to Stata native format
            graph save "`saving'.gph", `replace'
        }
        
        di as text "Graph saved to: `saving'"
    }
    
    restore
end
