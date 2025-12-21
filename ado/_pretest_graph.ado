*! _pretest_graph v0.1.0 - Event Study Plot for Pre-Test
*!
*! Description:
*!   Generates event study plot showing pre-treatment violations and 
*!   post-treatment DID estimates with confidence intervals.
*!   Follows Stata post-estimation command conventions (reads from e()).
*!
*! Syntax:
*!   _pretest_graph [, TItle(string) SAVing(string) NAME(string) REPLACE CONVentional]
*!
*! Options:
*!   title()      - Custom graph title
*!   saving()     - Save graph to file (.gph, .png, .pdf, .eps)
*!   name()       - Graph window name (for multiple graphs)
*!   replace      - Replace existing file
*!   conventional - Also show conventional CI (black) for comparison
*!
*! Data source: All from e() returns
*!   e(nu)           - Iterative violations vector
*!   e(delta)        - DID estimates vector  
*!   e(Sigma)        - Asymptotic covariance matrix
*!   e(threshold)    - Threshold M
*!   e(t0)           - Treatment time point
*!   e(T_pre)        - Pre-treatment periods
*!   e(T_post)       - Post-treatment periods
*!   e(N)            - Sample size
*!   e(alpha)        - Significance level
*!   e(mode)         - Mode ("iterative"/"overall")
*!   e(pretest_pass) - Pretest result (1=PASS, 0=FAIL)
*!   e(delta_bar)    - Average ATT
*!   e(ci_lower/upper) - Conditional CI bounds
*!   e(ci_conv_lower/upper) - Conventional CI bounds
*!
*! Reference:
*!   Mikhaeil & Harshaw (2025), "In Defense of the Pre-Test"
*!   arXiv:2510.26470

program define _pretest_graph
    version 17.0
    
    // Parse syntax - minimal parameters, read from e()
    syntax [, TItle(string) SAVing(string) NAME(string) REPLACE CONVentional]
    
    // ========================================
    // Step 1: Validate e() returns exist
    // ========================================
    
    // Check required e() returns
    capture confirm matrix e(nu)
    if _rc != 0 {
        di as error "Error: e(nu) not found. Run pretest command first."
        exit 301
    }
    
    capture confirm matrix e(delta)
    if _rc != 0 {
        di as error "Error: e(delta) not found. Run pretest command first."
        exit 301
    }
    
    capture confirm matrix e(Sigma)
    if _rc != 0 {
        di as error "Error: e(Sigma) not found. Run pretest command first."
        exit 301
    }
    
    // Check required scalars
    local required_scalars "t0 T_pre T_post N alpha threshold pretest_pass"
    foreach s of local required_scalars {
        capture confirm scalar e(`s')
        if _rc != 0 {
            di as error "Error: e(`s') not found. Run pretest command first."
            exit 301
        }
    }
    
    // ========================================
    // Step 2: Extract data from e()
    // ========================================
    
    // Scalars
    local t0 = e(t0)
    local T_pre = e(T_pre)
    local T_post = e(T_post)
    local n = e(N)
    local alpha = e(alpha)
    local threshold = e(threshold)
    local pretest_pass = e(pretest_pass)
    local delta_bar = e(delta_bar)
    
    // Mode
    local mode = e(mode)
    if "`mode'" == "" {
        local mode "iterative"
    }
    
    // CI bounds (may be missing if pretest failed)
    local ci_lower = e(ci_lower)
    local ci_upper = e(ci_upper)
    local ci_conv_lower = e(ci_conv_lower)
    local ci_conv_upper = e(ci_conv_upper)
    
    // Matrices
    tempname nu_mat delta_mat Sigma_mat
    matrix `nu_mat' = e(nu)
    matrix `delta_mat' = e(delta)
    matrix `Sigma_mat' = e(Sigma)
    
    // Dimensions
    local T_pre_minus_1 = rowsof(`nu_mat')
    local T_post_actual = rowsof(`delta_mat')
    
    // Critical value for CI
    local z_crit = invnormal(1 - `alpha'/2)
    
    // ========================================
    // Step 3: Prepare plot data
    // ========================================
    
    preserve
    clear
    
    // Total points: T_pre-1 (violations) + T_post (DID estimates)
    local total_points = `T_pre_minus_1' + `T_post_actual'
    quietly set obs `total_points'
    
    // Generate variables
    quietly {
        gen period = .
        gen estimate = .
        gen ci_lo = .
        gen ci_hi = .
        gen is_pre = .
        
        // --- Pre-treatment violations ---
        // Period mapping: t=2 to t=t0-1 -> relative time = 2-t0 to -1
        // In relative time: violations at -T_pre+1, ..., -1
        
        // For overall mode: need cumulative sum transformation
        if "`mode'" == "overall" {
            // Cumulative sum of violations: nu_bar_t = sum_{s=2}^{t} nu_s
            // Also need covariance transformation: Sigma^Delta = A * Sigma * A'
            tempname nu_bar A_mat Sigma_nu Sigma_nu_overall
            
            // Build cumulative sum matrix A (lower triangular of 1s)
            matrix `A_mat' = J(`T_pre_minus_1', `T_pre_minus_1', 0)
            forvalues i = 1/`T_pre_minus_1' {
                forvalues j = 1/`i' {
                    matrix `A_mat'[`i', `j'] = 1
                }
            }
            
            // Compute nu_bar = A * nu
            matrix `nu_bar' = `A_mat' * `nu_mat'
            
            // Extract nu part of Sigma (upper-left T_pre_minus_1 x T_pre_minus_1)
            matrix `Sigma_nu' = `Sigma_mat'[1..`T_pre_minus_1', 1..`T_pre_minus_1']
            
            // Transform: Sigma^Delta = A * Sigma_nu * A'
            matrix `Sigma_nu_overall' = `A_mat' * `Sigma_nu' * `A_mat''
            
            // Fill pre-treatment data (overall violations)
            forvalues i = 1/`T_pre_minus_1' {
                local rel_time = `i' - `T_pre'  // e.g., T_pre=4, i=1 -> -3
                replace period = `rel_time' in `i'
                
                local est = `nu_bar'[`i', 1]
                replace estimate = `est' in `i'
                
                // SE from transformed Sigma (Sigma is asymptotic cov, SE = sqrt(Sigma_ii/n))
                // Paper Assumption 1: sqrt(n)(θ̂-θ) → N(0,Σ), so Var(θ̂) = Σ/n
                local var_i = `Sigma_nu_overall'[`i', `i']
                local se_i = sqrt(`var_i' / `n')
                replace ci_lo = `est' - `z_crit' * `se_i' in `i'
                replace ci_hi = `est' + `z_crit' * `se_i' in `i'
                replace is_pre = 1 in `i'
            }
        }
        else {
            // Iterative mode: use nu directly
            forvalues i = 1/`T_pre_minus_1' {
                local rel_time = `i' - `T_pre'  // e.g., T_pre=4, i=1 -> -3
                replace period = `rel_time' in `i'
                
                local est = `nu_mat'[`i', 1]
                replace estimate = `est' in `i'
                
                // SE from Sigma diagonal (Sigma is asymptotic cov, SE = sqrt(Sigma_ii/n))
                // Paper Assumption 1: sqrt(n)(θ̂-θ) → N(0,Σ), so Var(θ̂) = Σ/n
                local var_i = `Sigma_mat'[`i', `i']
                local se_i = sqrt(`var_i' / `n')
                replace ci_lo = `est' - `z_crit' * `se_i' in `i'
                replace ci_hi = `est' + `z_crit' * `se_i' in `i'
                replace is_pre = 1 in `i'
            }
        }
        
        // --- Post-treatment DID estimates ---
        // Period mapping: t=t0 to t=T -> relative time = 0 to T_post-1
        forvalues i = 1/`T_post_actual' {
            local obs = `T_pre_minus_1' + `i'
            local rel_time = `i' - 1  // t0 -> 0, t0+1 -> 1, etc.
            replace period = `rel_time' in `obs'
            
            local est = `delta_mat'[`i', 1]
            replace estimate = `est' in `obs'
            
            // SE from Sigma (delta part starts at T_pre_minus_1 + 1)
            // Paper Assumption 1: sqrt(n)(θ̂-θ) → N(0,Σ), so Var(θ̂) = Σ/n
            local idx = `T_pre_minus_1' + `i'
            local var_i = `Sigma_mat'[`idx', `idx']
            local se_i = sqrt(`var_i' / `n')
            replace ci_lo = `est' - `z_crit' * `se_i' in `obs'
            replace ci_hi = `est' + `z_crit' * `se_i' in `obs'
            replace is_pre = 0 in `obs'
        }
    }
    
    // ========================================
    // Step 4: Build twoway graph
    // ========================================
    
    // Separate pre and post data
    tempvar pre_est pre_ci_lo pre_ci_hi post_est post_ci_lo post_ci_hi
    quietly {
        gen `pre_est' = estimate if is_pre == 1
        gen `pre_ci_lo' = ci_lo if is_pre == 1
        gen `pre_ci_hi' = ci_hi if is_pre == 1
        gen `post_est' = estimate if is_pre == 0
        gen `post_ci_lo' = ci_lo if is_pre == 0
        gen `post_ci_hi' = ci_hi if is_pre == 0
    }
    
    // --- Graph title ---
    if "`title'" == "" {
        local title "Event Study with Pre-test (Mikhaeil-Harshaw 2025)"
    }
    
    // --- Mode label ---
    if "`mode'" == "overall" {
        local mode_label "Overall violations"
    }
    else {
        local mode_label "Iterative violations"
    }
    
    // --- Pretest result label ---
    if `pretest_pass' == 1 {
        local pretest_label "Pretest: PASS"
        local post_marker "O"
        local post_color "maroon"
    }
    else {
        local pretest_label "Pretest: FAIL"
        local post_marker "Oh"
        local post_color "maroon%50"  // Faded color for FAIL
    }
    
    // --- X-axis range (integer ticks) ---
    quietly summarize period
    local x_min_data = r(min)
    local x_max_data = r(max)
    local x_min = floor(`x_min_data') - 1
    local x_max = ceil(`x_max_data') + 1
    
    // --- Y-axis: include threshold lines ---
    quietly summarize estimate
    local y_min = min(r(min), -`threshold') - 0.1 * abs(r(max) - r(min))
    local y_max = max(r(max), `threshold') + 0.1 * abs(r(max) - r(min))
    
    // Handle edge case where range is zero
    if `y_min' == `y_max' {
        local y_min = `y_min' - 1
        local y_max = `y_max' + 1
    }
    
    // --- Build graph command ---
    local graph_cmd ""
    
    // Pre-treatment CI (navy) - thicker lines for visibility
    local graph_cmd `"`graph_cmd' (rcap `pre_ci_lo' `pre_ci_hi' period, lcolor(navy) lwidth(medthick))"'
    
    // Pre-treatment points (navy, solid circle)
    local graph_cmd `"`graph_cmd' (scatter `pre_est' period, mcolor(navy) msymbol(O) msize(medlarge))"'
    
    // Post-treatment CI
    if `pretest_pass' == 1 {
        // PASS: normal maroon
        local graph_cmd `"`graph_cmd' (rcap `post_ci_lo' `post_ci_hi' period, lcolor(maroon) lwidth(medthick))"'
    }
    else {
        // FAIL: faded with dashed lines
        local graph_cmd `"`graph_cmd' (rcap `post_ci_lo' `post_ci_hi' period, lcolor(maroon%50) lwidth(medthick) lpattern(dash))"'
    }
    
    // Post-treatment points (maroon)
    local graph_cmd `"`graph_cmd' (scatter `post_est' period, mcolor(`post_color') msymbol(`post_marker') msize(medlarge))"'
    
    // --- ATT with CI comparison ---
    // Always show conventional CI for comparison, conditional CI only if PASS
    local att_pos = `x_max' - 0.3
    
    // Create ATT marker data
    tempvar att_period att_est att_ci_lo att_ci_hi
    quietly {
        gen `att_period' = `att_pos' in 1
        gen `att_est' = `delta_bar' in 1
    }
    
    // Conventional CI (always shown, black, slightly offset)
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
    
    // Conditional CI (only if PASS, orange, prominent)
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
    
    // --- Reference lines and formatting ---
    local note_text "M = `: di %6.3f `threshold''  |  `pretest_label'"
    
    // Legend entries - simplified labels, multi-row layout
    if "`mode'" == "overall" {
        local pre_legend "Overall viol."
    }
    else {
        local pre_legend "Iter. viol."
    }
    
    // Legend setup for CI comparison
    local legend_order "2 4"
    local legend_labels `"label(2 "`pre_legend'") label(4 "DID est.")"'
    local legend_rows = 1
    
    // Always add conventional CI to legend if available
    if !missing(`ci_conv_lower') {
        local legend_order "2 4 6"
        local legend_labels `"label(2 "`pre_legend'") label(4 "DID est.") label(6 "Conv. CI")"'
    }
    
    // Add conditional CI to legend if PASS
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
    
    // --- Build name option ---
    local name_opt ""
    if "`name'" != "" {
        local name_opt "name(`name', replace)"
    }
    
    // --- Execute graph ---
    twoway `graph_cmd', ///
        xline(-0.5, lcolor(gs10) lpattern(dash) lwidth(medium)) ///
        yline(0, lcolor(gs6) lwidth(thin)) ///
        yline(`threshold', lcolor(orange%70) lpattern(shortdash)) ///
        yline(-`threshold', lcolor(orange%70) lpattern(shortdash)) ///
        xlabel(`x_min'(1)`x_max', labsize(small)) ///
        ylabel(, angle(0) labsize(small) format(%9.2f)) ///
        xtitle("Time relative to treatment", size(small)) ///
        ytitle("Estimate", size(small)) ///
        title("`title'", size(medium)) ///
        note("`note_text'", size(vsmall)) ///
        legend(order(`legend_order') `legend_labels' pos(6) rows(1) size(vsmall) symxsize(*.5)) ///
        graphregion(color(white)) plotregion(color(white)) ///
        `name_opt'
    
    // ========================================
    // Step 5: Handle saving
    // ========================================
    
    if "`saving'" != "" {
        // Determine file extension
        local ext = substr("`saving'", -4, 4)
        
        if "`ext'" == ".gph" {
            graph save "`saving'", `replace'
        }
        else if "`ext'" == ".png" | "`ext'" == ".pdf" | "`ext'" == ".eps" {
            graph export "`saving'", `replace'
        }
        else {
            // Default to .gph if no recognized extension
            graph save "`saving'.gph", `replace'
        }
        
        di as text "Graph saved to: `saving'"
    }
    
    restore
end
